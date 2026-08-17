#!/usr/bin/env python3

import argparse
import copy
import os
import re
import sys
from collections import defaultdict

import gemmi
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Geometry import Point3D


WATER_NAMES = {
    "HOH", "WAT", "H2O", "DOD",
}


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Read a PDB/mmCIF structure, optionally merge ligand(s) from "
            ".sdf/.mol/.mol2, canonicalize atom names for all non-polymer, "
            "non-water residues (including ions), discard hydrogens, optionally "
            "reassign ligand residue name/chain/sequence numbering, preserve "
            "connectivity for externally added ligands in mmCIF via _chem_comp_bond, "
            "and write mmCIF. If the ligand file contains multiple molecules, one "
            "output mmCIF is generated per ligand model."
        )
    )
    p.add_argument("input", help="Input structure file (.pdb, .ent, .cif, .mmcif)")
    p.add_argument("output", help="Output mmCIF path")
    p.add_argument(
        "--ligand-file",
        help="Optional external ligand file (.sdf, .mol, .mol2) to merge into the structure",
    )
    p.add_argument(
        "--ligand-resname",
        default=None,
        help="Optional residue name for ligands",
    )
    p.add_argument(
        "--ligand-chain",
        default=None,
        help="Optional chain name for ligands",
    )
    p.add_argument(
        "--ligand-seq",
        type=int,
        default=None,
        help="Optional starting residue sequence number for ligands",
    )
    p.add_argument(
        "--model-index",
        type=int,
        default=0,
        help="Model index where external ligand will be inserted (default: 0)",
    )
    p.add_argument(
        "--remove-ligands-from-receptor",
        action="store_true",
        help="Remove existing ligand residues from the input receptor structure before processing",
    )
    p.add_argument(
        "--verbose",
        action="store_true",
        help="Print processing details",
    )
    return p.parse_args()


def read_structure(path: str) -> gemmi.Structure:
    return gemmi.read_structure(path)


def residue_is_polymer(res: gemmi.Residue) -> bool:
    info = gemmi.find_tabulated_residue(res.name)
    if info is None:
        return False
    return info.is_amino_acid() or info.is_nucleic_acid()


def residue_is_water(res: gemmi.Residue) -> bool:
    return res.name.strip().upper() in WATER_NAMES


def is_ligand_residue(res: gemmi.Residue) -> bool:
    if residue_is_polymer(res):
        return False
    if residue_is_water(res):
        return False

    atoms = list(res)
    if not atoms:
        return False

    heavy_atoms = [a for a in atoms if a.element.name != "H"]
    return len(heavy_atoms) > 0


def discard_hydrogens_from_residue(res: gemmi.Residue) -> int:
    atom_indices_to_remove = [i for i, atom in enumerate(res) if atom.element.name == "H"]
    if not atom_indices_to_remove:
        return 0
    for i in reversed(atom_indices_to_remove):
        del res[i]
    return len(atom_indices_to_remove)


def discard_hydrogens_from_structure(st: gemmi.Structure) -> int:
    total_removed = 0
    for model in st:
        for chain in model:
            for res in chain:
                total_removed += discard_hydrogens_from_residue(res)
    return total_removed


def remove_ligands_from_structure(st: gemmi.Structure) -> int:
    removed = 0
    for model in st:
        rebuilt_chains = []
        for chain in model:
            new_chain = gemmi.Chain(chain.name)
            for res in chain:
                if is_ligand_residue(res):
                    removed += 1
                    continue
                new_chain.add_residue(copy.deepcopy(res))
            rebuilt_chains.append(new_chain)

        for i in reversed(range(len(model))):
            del model[i]
        for chain in rebuilt_chains:
            if len(chain) > 0:
                model.add_chain(chain)
    return removed


def gemmi_atom_charge(atom: gemmi.Atom) -> int:
    try:
        return int(atom.charge)
    except Exception:
        return 0


def build_rdkit_mol_from_residue(res: gemmi.Residue):
    rw = Chem.RWMol()
    conf = Chem.Conformer()
    selected_atoms = []
    total_charge = 0

    for atom in res:
        if atom.element.name == "H":
            continue

        rd_atom = Chem.Atom(atom.element.atomic_number)
        chg = gemmi_atom_charge(atom)
        if chg != 0:
            rd_atom.SetFormalCharge(chg)

        idx = rw.AddAtom(rd_atom)
        pos = atom.pos
        conf.SetAtomPosition(idx, Point3D(float(pos.x), float(pos.y), float(pos.z)))

        selected_atoms.append(atom)
        total_charge += chg

    mol = rw.GetMol()
    conf.SetId(0)
    mol.AddConformer(conf, assignId=True)

    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    if mol.GetNumAtoms() == 0:
        return mol, selected_atoms

    if mol.GetNumAtoms() == 1:
        return mol, selected_atoms

    try:
        rdDetermineBonds.DetermineBonds(mol, charge=total_charge)
    except Exception:
        try:
            rdDetermineBonds.DetermineConnectivity(mol)
        except Exception:
            pass

    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    return mol, selected_atoms


def canonical_names_for_mol(mol: Chem.Mol):
    n = mol.GetNumAtoms()
    if n == 0:
        return []

    try:
        mol.UpdatePropertyCache(strict=False)

        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass

        ranks = list(
            Chem.CanonicalRankAtoms(
                mol,
                breakTies=True,
                includeChirality=True,
                includeIsotopes=False,
                includeAtomMaps=False,
                includeChiralPresence=True,
            )
        )
        order = sorted(range(n), key=lambda i: ranks[i])

    except Exception:
        conf = mol.GetConformer()
        order = sorted(
            range(n),
            key=lambda i: (
                mol.GetAtomWithIdx(i).GetAtomicNum(),
                round(conf.GetAtomPosition(i).x, 3),
                round(conf.GetAtomPosition(i).y, 3),
                round(conf.GetAtomPosition(i).z, 3),
                i,
            ),
        )

    counters = defaultdict(int)
    names = [None] * n

    for i in order:
        atom = mol.GetAtomWithIdx(i)
        symbol = atom.GetSymbol().upper()
        counters[symbol] += 1
        names[i] = f"{symbol}{counters[symbol]}"

    return names


def rename_residue_atoms(res: gemmi.Residue):
    mol, selected_atoms = build_rdkit_mol_from_residue(res)
    if mol.GetNumAtoms() == 0:
        return 0

    new_names = canonical_names_for_mol(mol)

    renamed = 0
    for atom, new_name in zip(selected_atoms, new_names):
        if atom.name != new_name:
            atom.name = new_name
            renamed += 1

    return renamed


def get_mol_name(mol: Chem.Mol, fallback: str) -> str:
    for key in ("_Name", "NAME", "name", "Title", "title"):
        if mol.HasProp(key):
            value = mol.GetProp(key).strip()
            if value:
                return value
    return fallback


def sanitize_suffix(text: str) -> str:
    text = text.strip()
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9._-]", "_", text)
    text = re.sub(r"_+", "_", text)
    text = text.strip("._")
    return text or "model"


def derive_output_path(base_output_path: str, suffix: str) -> str:
    root, ext = os.path.splitext(base_output_path)
    if not ext:
        ext = ".cif"
    return f"{root}_{suffix}{ext}"


def _safe_int(value, default=0):
    try:
        return int(value)
    except Exception:
        return default


def _safe_float(value, default=0.0):
    try:
        return float(value)
    except Exception:
        return default


def _element_from_mol2_type(atom_type: str, atom_name: str = "") -> str:
    pt = Chem.GetPeriodicTable()

    candidates = []

    if atom_type:
        t = atom_type.strip()
        part = t.split(".")[0]
        candidates.append(part)
        m = re.match(r"([A-Za-z]{1,2})", t)
        if m:
            candidates.append(m.group(1))

    if atom_name:
        m = re.match(r"([A-Za-z]{1,2})", atom_name.strip())
        if m:
            candidates.append(m.group(1))

    for c in candidates:
        if not c:
            continue
        symbol = c[0].upper() + c[1:].lower()
        try:
            z = pt.GetAtomicNumber(symbol)
            if z > 0:
                return symbol
        except Exception:
            pass

    raise RuntimeError(
        f"Could not infer element from MOL2 atom type '{atom_type}' / atom name '{atom_name}'"
    )


def _bond_type_from_mol2(bond_type: str):
    bt = (bond_type or "").strip().lower()
    if bt == "1":
        return Chem.BondType.SINGLE, False
    if bt == "2":
        return Chem.BondType.DOUBLE, False
    if bt == "3":
        return Chem.BondType.TRIPLE, False
    if bt == "ar":
        return Chem.BondType.AROMATIC, True
    if bt == "am":
        return Chem.BondType.SINGLE, False
    if bt in ("du", "un", "nc"):
        return Chem.BondType.SINGLE, False
    return Chem.BondType.SINGLE, False


def _parse_mol2_block(block_text: str):
    lines = block_text.splitlines()
    sections = {}
    current = None

    for line in lines:
        if line.startswith("@<TRIPOS>"):
            current = line.strip()
            sections[current] = []
        elif current is not None:
            sections[current].append(line.rstrip("\n"))

    if "@<TRIPOS>MOLECULE" not in sections or "@<TRIPOS>ATOM" not in sections:
        return None

    mol_name = "model"
    mol_lines = sections["@<TRIPOS>MOLECULE"]
    for line in mol_lines:
        stripped = line.strip()
        if stripped:
            mol_name = stripped
            break

    atom_lines = [x for x in sections.get("@<TRIPOS>ATOM", []) if x.strip()]
    bond_lines = [x for x in sections.get("@<TRIPOS>BOND", []) if x.strip()]

    rw = Chem.RWMol()
    conf = Chem.Conformer()
    id_to_idx = {}

    for line in atom_lines:
        parts = line.split()
        if len(parts) < 6:
            continue

        atom_id = _safe_int(parts[0], None)
        atom_name = parts[1]
        x = _safe_float(parts[2])
        y = _safe_float(parts[3])
        z = _safe_float(parts[4])
        atom_type = parts[5]

        if atom_id is None:
            continue

        symbol = _element_from_mol2_type(atom_type, atom_name)
        atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(symbol)

        rd_atom = Chem.Atom(atomic_num)

        idx = rw.AddAtom(rd_atom)
        conf.SetAtomPosition(idx, Point3D(float(x), float(y), float(z)))
        id_to_idx[atom_id] = idx

    for line in bond_lines:
        parts = line.split()
        if len(parts) < 4:
            continue

        a1 = _safe_int(parts[1], None)
        a2 = _safe_int(parts[2], None)
        btype_str = parts[3]

        if a1 not in id_to_idx or a2 not in id_to_idx:
            continue

        btype, aromatic = _bond_type_from_mol2(btype_str)
        i1 = id_to_idx[a1]
        i2 = id_to_idx[a2]

        try:
            rw.AddBond(i1, i2, btype)
            bond = rw.GetBondBetweenAtoms(i1, i2)
            if bond is not None and aromatic:
                bond.SetIsAromatic(True)
                rw.GetAtomWithIdx(i1).SetIsAromatic(True)
                rw.GetAtomWithIdx(i2).SetIsAromatic(True)
        except Exception:
            pass

    mol = rw.GetMol()
    conf.SetId(0)
    mol.AddConformer(conf, assignId=True)
    try:
        mol.SetProp("_Name", mol_name)
    except Exception:
        pass
    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    return mol


def _standardize_external_mol(mol: Chem.Mol) -> Chem.Mol:
    """
    Convert any imported external molecule (SDF/MOL/MOL2) into the same internal
    simplified representation:
      - hydrogens removed
      - explicit heavy atoms
      - explicit bonds copied when available
      - coordinates preserved
      - imported formal charges ignored for consistency across formats
    """
    rw = Chem.RWMol()
    conf_out = Chem.Conformer()
    idx_map = {}

    conf_in = mol.GetConformer()

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue

        new_atom = Chem.Atom(atom.GetAtomicNum())
        try:
            new_atom.SetFormalCharge(0)
        except Exception:
            pass
        try:
            new_atom.SetIsAromatic(atom.GetIsAromatic())
        except Exception:
            pass

        new_idx = rw.AddAtom(new_atom)
        pos = conf_in.GetAtomPosition(atom.GetIdx())
        conf_out.SetAtomPosition(new_idx, Point3D(float(pos.x), float(pos.y), float(pos.z)))
        idx_map[atom.GetIdx()] = new_idx

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if a1 not in idx_map or a2 not in idx_map:
            continue

        btype = bond.GetBondType()
        if btype not in (
            Chem.BondType.SINGLE,
            Chem.BondType.DOUBLE,
            Chem.BondType.TRIPLE,
            Chem.BondType.AROMATIC,
        ):
            btype = Chem.BondType.SINGLE

        try:
            rw.AddBond(idx_map[a1], idx_map[a2], btype)
            new_bond = rw.GetBondBetweenAtoms(idx_map[a1], idx_map[a2])
            if new_bond is not None and btype == Chem.BondType.AROMATIC:
                new_bond.SetIsAromatic(True)
                rw.GetAtomWithIdx(idx_map[a1]).SetIsAromatic(True)
                rw.GetAtomWithIdx(idx_map[a2]).SetIsAromatic(True)
        except Exception:
            pass

    out = rw.GetMol()
    conf_out.SetId(0)
    out.AddConformer(conf_out, assignId=True)

    try:
        if mol.HasProp("_Name"):
            out.SetProp("_Name", mol.GetProp("_Name"))
    except Exception:
        pass

    try:
        out.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    return out


def read_external_ligands(path: str):
    lower = path.lower()
    molecules = []

    if lower.endswith(".sdf"):
        suppl = Chem.SDMolSupplier(path, removeHs=False, sanitize=False)
        valid_idx = 0
        for mol in suppl:
            if mol is None:
                continue
            if mol.GetNumAtoms() == 0 or mol.GetNumConformers() == 0:
                continue
            mol = _standardize_external_mol(mol)
            valid_idx += 1
            name = get_mol_name(mol, f"model{valid_idx}")
            molecules.append((name, mol))

    elif lower.endswith(".mol"):
        mol = Chem.MolFromMolFile(path, removeHs=False, sanitize=False)
        if mol is None:
            raise RuntimeError(f"Could not read MOL file: {path}")
        if mol.GetNumAtoms() == 0 or mol.GetNumConformers() == 0:
            raise RuntimeError(f"MOL file does not contain usable atoms/coordinates: {path}")
        mol = _standardize_external_mol(mol)
        name = get_mol_name(mol, "model1")
        molecules.append((name, mol))

    elif lower.endswith(".mol2"):
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read()

        parts = text.split("@<TRIPOS>MOLECULE")
        blocks = ["@<TRIPOS>MOLECULE" + part for part in parts[1:]]

        if not blocks:
            raise RuntimeError(f"Could not read any molecule from MOL2 file: {path}")

        valid_idx = 0
        for block in blocks:
            mol = _parse_mol2_block(block)
            if mol is None:
                continue
            if mol.GetNumAtoms() == 0 or mol.GetNumConformers() == 0:
                continue
            mol = _standardize_external_mol(mol)
            valid_idx += 1
            name = get_mol_name(mol, f"model{valid_idx}")
            molecules.append((name, mol))

        if not molecules:
            raise RuntimeError(f"Could not parse any valid molecule from MOL2 file: {path}")

    else:
        raise RuntimeError("Ligand file must be .sdf, .mol, or .mol2")

    for name, mol in molecules:
        if mol.GetNumAtoms() == 0:
            raise RuntimeError(f"Ligand model '{name}' contains no heavy atoms")
        if mol.GetNumConformers() == 0:
            raise RuntimeError(f"Ligand model '{name}' does not contain coordinates")

    return molecules


def find_or_create_chain(model: gemmi.Model, chain_name: str) -> gemmi.Chain:
    for chain in model:
        if chain.name == chain_name:
            return chain

    model.add_chain(gemmi.Chain(chain_name))

    for chain in model:
        if chain.name == chain_name:
            return chain

    raise RuntimeError(f"Failed to create chain '{chain_name}'")


def element_from_atomic_num(atomic_num: int) -> gemmi.Element:
    symbol = Chem.GetPeriodicTable().GetElementSymbol(atomic_num)
    return gemmi.Element(symbol)


def relabel_existing_ligands(
    st: gemmi.Structure,
    ligand_resname: str | None,
    ligand_chain: str | None,
    ligand_seq_start: int | None,
    verbose: bool = False,
):
    total_count = 0

    if ligand_resname is None and ligand_chain is None and ligand_seq_start is None:
        return 0

    for model_index, model in enumerate(st):
        next_seq = ligand_seq_start

        rebuilt_chains = {}
        chain_order = []

        def get_or_create_chain(chain_name: str) -> gemmi.Chain:
            if chain_name not in rebuilt_chains:
                rebuilt_chains[chain_name] = gemmi.Chain(chain_name)
                chain_order.append(chain_name)
            return rebuilt_chains[chain_name]

        for chain in model:
            get_or_create_chain(chain.name)

            for res in chain:
                if not is_ligand_residue(res):
                    get_or_create_chain(chain.name).add_residue(copy.deepcopy(res))
                    continue

                old_chain_name = chain.name
                old_resname = res.name
                old_seqid = str(res.seqid)

                new_chain_name = ligand_chain if ligand_chain is not None else chain.name
                new_resname = ligand_resname if ligand_resname is not None else res.name

                new_res = gemmi.Residue()
                new_res.name = new_resname
                new_res.het_flag = "H"
                try:
                    new_res.subchain = new_chain_name
                except Exception:
                    pass

                if ligand_seq_start is not None:
                    new_res.seqid = gemmi.SeqId(str(next_seq))
                    new_seqid = str(next_seq)
                    next_seq += 1
                else:
                    try:
                        new_res.seqid = copy.deepcopy(res.seqid)
                    except Exception:
                        new_res.seqid = gemmi.SeqId(str(res.seqid))
                    new_seqid = old_seqid

                for atom in res:
                    new_res.add_atom(copy.deepcopy(atom))

                get_or_create_chain(new_chain_name).add_residue(new_res)
                total_count += 1

                if verbose:
                    print(
                        f"[RELABEL] model{model_index}: "
                        f"{old_chain_name}/{old_resname}/{old_seqid} -> "
                        f"{new_chain_name}/{new_resname}/{new_seqid}",
                        file=sys.stderr,
                    )

        for i in reversed(range(len(model))):
            del model[i]

        for chain_name in chain_order:
            chain = rebuilt_chains[chain_name]
            if len(chain) > 0:
                model.add_chain(chain)

    return total_count


def next_available_ligand_seqnum(
    st: gemmi.Structure,
    model_index: int,
    chain_name: str,
    requested_start: int,
) -> int:
    if model_index < 0 or model_index >= len(st):
        return requested_start

    used = set()
    model = st[model_index]
    for chain in model:
        if chain.name != chain_name:
            continue
        for res in chain:
            if is_ligand_residue(res):
                try:
                    used.add(int(res.seqid.num))
                except Exception:
                    try:
                        used.add(int(str(res.seqid)))
                    except Exception:
                        pass

    seq = requested_start
    while seq in used:
        seq += 1
    return seq


def add_external_ligand_to_structure(
    st: gemmi.Structure,
    ligand_mol: Chem.Mol,
    resname: str = "LIG",
    chain_name: str = "L",
    seqnum: int = 1,
    model_index: int = 0,
):
    if model_index < 0 or model_index >= len(st):
        raise RuntimeError(f"Model index {model_index} is out of range")

    model = st[model_index]
    chain = find_or_create_chain(model, chain_name)

    res = gemmi.Residue()
    res.name = resname
    res.seqid = gemmi.SeqId(str(seqnum))
    res.het_flag = "H"
    try:
        res.subchain = chain_name
    except Exception:
        pass

    canonical_names = canonical_names_for_mol(ligand_mol)
    conf = ligand_mol.GetConformer()

    heavy_atom_indices = [a.GetIdx() for a in ligand_mol.GetAtoms() if a.GetAtomicNum() != 1]

    if len(heavy_atom_indices) != len(canonical_names):
        raise RuntimeError(
            f"Internal error: heavy atom count ({len(heavy_atom_indices)}) does not "
            f"match canonical name count ({len(canonical_names)})"
        )

    for rd_idx, atom_name in zip(heavy_atom_indices, canonical_names):
        atom = ligand_mol.GetAtomWithIdx(rd_idx)
        atomic_num = atom.GetAtomicNum()

        gatom = gemmi.Atom()
        gatom.element = element_from_atomic_num(atomic_num)

        pos = conf.GetAtomPosition(rd_idx)
        gatom.pos = gemmi.Position(float(pos.x), float(pos.y), float(pos.z))
        gatom.occ = 1.0
        gatom.b_iso = 20.0

        try:
            gatom.charge = 0
        except Exception:
            pass

        gatom.name = atom_name
        res.add_atom(gatom)

    chain.add_residue(res)
    return model_index, chain.name, res.name, str(res.seqid)


def process_structure(
    st: gemmi.Structure,
    verbose: bool = False,
    skip_residues=None,
):
    if skip_residues is None:
        skip_residues = set()

    total_ligands = 0
    total_atoms_renamed = 0
    failures = []

    for model_index, model in enumerate(st):
        for chain in model:
            for res in chain:
                if not is_ligand_residue(res):
                    continue

                total_ligands += 1
                res_id = f"model{model_index}/{chain.name}/{res.name} {res.seqid}"
                res_key = (model_index, chain.name, res.name, str(res.seqid))

                if res_key in skip_residues:
                    if verbose:
                        print(
                            f"[SKIP] {res_id}: external ligand already canonically named",
                            file=sys.stderr,
                        )
                    continue

                try:
                    removed_h = discard_hydrogens_from_residue(res)
                    n = rename_residue_atoms(res)
                    total_atoms_renamed += n
                    if verbose:
                        if removed_h > 0:
                            print(
                                f"[OK] {res_id}: removed {removed_h} H atoms, renamed {n} atoms",
                                file=sys.stderr,
                            )
                        else:
                            print(f"[OK] {res_id}: renamed {n} atoms", file=sys.stderr)
                except Exception as e:
                    failures.append((res_id, str(e)))
                    if verbose:
                        print(f"[FAIL] {res_id}: {e}", file=sys.stderr)

    return total_ligands, total_atoms_renamed, failures


def _seqid_num_as_str(seqid) -> str:
    try:
        return str(int(seqid.num))
    except Exception:
        m = re.match(r"^\s*(-?\d+)", str(seqid).strip())
        if m:
            return m.group(1)
    return "?"


def _seqid_ins_code_as_str(seqid) -> str:
    try:
        icode = str(seqid.icode).strip()
        if icode and icode not in ("?", ".", "\x00"):
            return icode
    except Exception:
        pass

    s = str(seqid).strip()
    m = re.match(r"^\s*-?\d+([A-Za-z])\s*$", s)
    if m:
        return m.group(1)

    return "?"


def rewrite_atom_site_identifiers_from_structure(
    doc: gemmi.cif.Document,
    st: gemmi.Structure,
) -> None:
    block = doc.sole_block()

    tags = [
        "label_comp_id",
        "label_asym_id",
        "label_seq_id",
        "auth_seq_id",
        "auth_asym_id",
        "pdbx_PDB_ins_code",
    ]

    table = block.find("_atom_site.", tags)
    if not table:
        return

    atom_records = []
    for model in st:
        for chain in model:
            for res in chain:
                seq_num = _seqid_num_as_str(res.seqid)
                ins_code = _seqid_ins_code_as_str(res.seqid)
                for atom in res:
                    atom_records.append(
                        (
                            res.name,
                            chain.name,
                            seq_num,
                            seq_num,
                            chain.name,
                            ins_code,
                        )
                    )

    rows = list(table)
    if len(rows) != len(atom_records):
        raise RuntimeError(
            f"_atom_site row count ({len(rows)}) does not match structure atom count ({len(atom_records)})"
        )

    label_comp_idx = tags.index("label_comp_id")
    label_asym_idx = tags.index("label_asym_id")
    label_seq_idx = tags.index("label_seq_id")
    auth_seq_idx = tags.index("auth_seq_id")
    auth_asym_idx = tags.index("auth_asym_id")
    ins_code_idx = tags.index("pdbx_PDB_ins_code")

    for row, (
        resname,
        chain_name,
        label_seq_id,
        auth_seq_id,
        auth_asym_id,
        ins_code,
    ) in zip(rows, atom_records):
        row[label_comp_idx] = resname
        row[label_asym_idx] = chain_name
        row[label_seq_idx] = label_seq_id
        row[auth_seq_idx] = auth_seq_id
        row[auth_asym_idx] = auth_asym_id
        row[ins_code_idx] = ins_code


def _bond_order_as_cif_value(bond: Chem.Bond) -> str:
    bt = bond.GetBondType()
    if bt == Chem.BondType.SINGLE:
        return "sing"
    if bt == Chem.BondType.DOUBLE:
        return "doub"
    if bt == Chem.BondType.TRIPLE:
        return "trip"
    if bt == Chem.BondType.AROMATIC:
        return "arom"
    return "sing"


def add_chem_comp_bond_loop(block: gemmi.cif.Block, comp_id: str, mol: Chem.Mol) -> int:
    """
    Add _chem_comp_bond records for one ligand component using the same canonical
    atom naming scheme that is used when atoms are inserted into the Gemmi residue.
    Only heavy-atom bonds are written, because hydrogens are removed from imported
    ligands earlier in the pipeline.
    """
    canonical_names = canonical_names_for_mol(mol)
    heavy_atom_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]

    if len(heavy_atom_indices) != len(canonical_names):
        raise RuntimeError(
            f"Internal error: heavy atom count ({len(heavy_atom_indices)}) does not "
            f"match canonical name count ({len(canonical_names)}) for _chem_comp_bond"
        )

    idx_to_name = dict(zip(heavy_atom_indices, canonical_names))

    loop = block.init_loop(
        "_chem_comp_bond.",
        [
            "comp_id",
            "atom_id_1",
            "atom_id_2",
            "value_order",
            "pdbx_aromatic_flag",
        ],
    )

    n_rows = 0
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()

        if a1 not in idx_to_name or a2 not in idx_to_name:
            continue

        loop.add_row(
            [
                comp_id,
                idx_to_name[a1],
                idx_to_name[a2],
                _bond_order_as_cif_value(bond),
                "Y" if bond.GetIsAromatic() else "N",
            ]
        )
        n_rows += 1

    return n_rows


def write_mmcif(
    st: gemmi.Structure,
    output_path: str,
    ligand_mol: Chem.Mol | None = None,
    ligand_resname: str | None = None,
):
    doc = st.make_mmcif_document()
    rewrite_atom_site_identifiers_from_structure(doc, st)
    block = doc.sole_block()

    added_bonds = 0
    if ligand_mol is not None and ligand_resname is not None:
        added_bonds = add_chem_comp_bond_loop(block, ligand_resname, ligand_mol)

    doc.write_file(output_path)
    return added_bonds


def process_one_output(
    base_structure: gemmi.Structure,
    output_path: str,
    args,
    ligand_name=None,
    ligand_mol=None,
):
    st = copy.deepcopy(base_structure)

    removed_h_global = discard_hydrogens_from_structure(st)

    removed_existing_ligands = 0
    if args.remove_ligands_from_receptor:
        removed_existing_ligands = remove_ligands_from_structure(st)

    relabeled_existing = relabel_existing_ligands(
        st,
        ligand_resname=args.ligand_resname,
        ligand_chain=args.ligand_chain,
        ligand_seq_start=args.ligand_seq,
        verbose=args.verbose,
    )

    merged_info = None
    ext_resname = None

    if ligand_mol is not None:
        ext_resname = args.ligand_resname if args.ligand_resname is not None else "LIG"
        ext_chain = args.ligand_chain if args.ligand_chain is not None else "L"
        ext_seq_start = args.ligand_seq if args.ligand_seq is not None else 1
        ext_seq = next_available_ligand_seqnum(
            st=st,
            model_index=args.model_index,
            chain_name=ext_chain,
            requested_start=ext_seq_start,
        )

        merged_info = add_external_ligand_to_structure(
            st,
            ligand_mol,
            resname=ext_resname,
            chain_name=ext_chain,
            seqnum=ext_seq,
            model_index=args.model_index,
        )
        if args.verbose:
            print(
                f"[MERGED] ligand model '{ligand_name}' -> "
                f"model_index={merged_info[0]} chain={merged_info[1]} "
                f"residue={merged_info[2]} seqid={merged_info[3]}",
                file=sys.stderr,
            )

    skip_residues = set()
    if merged_info is not None:
        skip_residues.add(merged_info)

    total_ligands, total_atoms_renamed, failures = process_structure(
        st,
        verbose=args.verbose,
        skip_residues=skip_residues,
    )

    added_bond_rows = write_mmcif(
        st,
        output_path,
        ligand_mol=ligand_mol,
        ligand_resname=ext_resname,
    )

    print(f"Input structure:   {args.input}", file=sys.stderr)
    if args.ligand_file:
        print(f"Input ligand file: {args.ligand_file}", file=sys.stderr)
    if ligand_name is not None:
        print(f"Ligand model:      {ligand_name}", file=sys.stderr)
    print(f"Output:            {output_path}", file=sys.stderr)
    print(f"Ligands processed: {total_ligands}", file=sys.stderr)
    print(f"Atoms renamed:     {total_atoms_renamed}", file=sys.stderr)
    print(f"Hydrogens removed: {removed_h_global}", file=sys.stderr)
    print(f"Existing ligands removed: {removed_existing_ligands}", file=sys.stderr)
    print(f"Existing ligands relabeled: {relabeled_existing}", file=sys.stderr)
    print(f"External ligand bond rows written: {added_bond_rows}", file=sys.stderr)

    if failures:
        print(f"Residues failed:   {len(failures)}", file=sys.stderr)
        for res_id, msg in failures:
            print(f"  - {res_id}: {msg}", file=sys.stderr)


def main():
    args = parse_args()

    base_structure = read_structure(args.input)

    if not args.ligand_file:
        process_one_output(
            base_structure=base_structure,
            output_path=args.output,
            args=args,
            ligand_name=None,
            ligand_mol=None,
        )
        return

    ligand_models = read_external_ligands(args.ligand_file)

    if len(ligand_models) == 1:
        ligand_name, ligand_mol = ligand_models[0]

        input_ext = os.path.splitext(args.ligand_file)[1].lower()
        if input_ext in (".sdf", ".mol2"):
            suffix = sanitize_suffix(ligand_name)
            output_path = derive_output_path(args.output, suffix)
        else:
            output_path = args.output

        process_one_output(
            base_structure=base_structure,
            output_path=output_path,
            args=args,
            ligand_name=ligand_name,
            ligand_mol=ligand_mol,
        )
        return

    used_suffixes = set()
    for i, (ligand_name, ligand_mol) in enumerate(ligand_models, start=1):
        suffix = sanitize_suffix(ligand_name)
        if suffix in used_suffixes:
            suffix = f"{suffix}_{i}"
        used_suffixes.add(suffix)

        output_path = derive_output_path(args.output, suffix)
        process_one_output(
            base_structure=base_structure,
            output_path=output_path,
            args=args,
            ligand_name=ligand_name,
            ligand_mol=ligand_mol,
        )


if __name__ == "__main__":
    main()