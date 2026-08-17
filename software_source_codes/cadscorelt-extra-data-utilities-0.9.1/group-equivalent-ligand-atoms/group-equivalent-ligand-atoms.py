#!/usr/bin/env python3

import sys
import os
import shlex
from collections import defaultdict, Counter

import gemmi
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
}
NT = {"A", "C", "G", "U", "T", "DA", "DC", "DG", "DT"}


def is_polymer_resname(name: str) -> bool:
    name = (name or "").strip().upper()
    return name in AA or name in NT


def is_water(name: str) -> bool:
    return (name or "").strip().upper() in {"HOH", "WAT", "DOD"}


def normalize_atom_name(name: str) -> str:
    return (name or "").strip()


def cif_unquote(value: str) -> str:
    if value is None:
        return ""
    v = value.strip()
    if len(v) >= 2 and ((v[0] == "'" and v[-1] == "'") or (v[0] == '"' and v[-1] == '"')):
        return v[1:-1]
    return v


def tokenize_cif_line(line: str):
    """
    Lightweight tokenizer for ordinary mmCIF loop rows.
    Good enough for _chem_comp_bond / _chem_comp_atom rows.
    """
    line = line.rstrip("\n")
    if not line.strip():
        return []
    if line.lstrip().startswith("#"):
        return []
    lex = shlex.shlex(line, posix=True)
    lex.whitespace_split = True
    lex.commenters = ""
    return list(lex)


def parse_mmcif_loops(path: str):
    """
    Minimal loop parser for mmCIF files.
    Returns a list of loops, each as (tags, rows), where:
      tags = [tag1, tag2, ...]
      rows = [[v11, v12, ...], [v21, v22, ...], ...]
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    loops = []
    i = 0
    n = len(lines)

    while i < n:
        line = lines[i].strip()
        if line != "loop_":
            i += 1
            continue

        i += 1
        tags = []
        while i < n:
            s = lines[i].strip()
            if s.startswith("_"):
                tags.append(s.split()[0])
                i += 1
            else:
                break

        if not tags:
            continue

        rows = []
        pending = []
        width = len(tags)

        while i < n:
            raw = lines[i]
            s = raw.strip()

            # end of loop/data block
            if not s:
                i += 1
                continue
            if s == "loop_" or s.startswith("data_"):
                break
            if s.startswith("_"):
                break
            if s.startswith("#"):
                i += 1
                break

            tokens = tokenize_cif_line(raw)
            if tokens:
                pending.extend(tokens)
                while len(pending) >= width:
                    row = pending[:width]
                    pending = pending[width:]
                    rows.append(row)
            i += 1

        loops.append((tags, rows))

    return loops


def load_chemcomp_tables_from_mmcif(path: str):
    """
    Extract:
      chem_comp_bonds[comp_id] = list of bond dicts
      chem_comp_atoms[comp_id][atom_id] = atom dict
    """
    chem_comp_bonds = defaultdict(list)
    chem_comp_atoms = defaultdict(dict)

    if not os.path.isfile(path):
        return chem_comp_bonds, chem_comp_atoms

    lower = path.lower()
    if not (lower.endswith(".cif") or lower.endswith(".mmcif")):
        return chem_comp_bonds, chem_comp_atoms

    try:
        loops = parse_mmcif_loops(path)
    except Exception:
        return chem_comp_bonds, chem_comp_atoms

    for tags, rows in loops:
        tagset = set(tags)

        # _chem_comp_bond loop
        needed_bond = {
            "_chem_comp_bond.comp_id",
            "_chem_comp_bond.atom_id_1",
            "_chem_comp_bond.atom_id_2",
        }
        if needed_bond.issubset(tagset):
            idx = {tag: j for j, tag in enumerate(tags)}
            for row in rows:
                comp_id = cif_unquote(row[idx["_chem_comp_bond.comp_id"]]).upper()
                atom1 = normalize_atom_name(cif_unquote(row[idx["_chem_comp_bond.atom_id_1"]]))
                atom2 = normalize_atom_name(cif_unquote(row[idx["_chem_comp_bond.atom_id_2"]]))
                order = cif_unquote(row[idx["_chem_comp_bond.value_order"]]) if "_chem_comp_bond.value_order" in idx else ""
                aromatic_flag = cif_unquote(row[idx["_chem_comp_bond.pdbx_aromatic_flag"]]) if "_chem_comp_bond.pdbx_aromatic_flag" in idx else ""
                chem_comp_bonds[comp_id].append({
                    "atom_id_1": atom1,
                    "atom_id_2": atom2,
                    "value_order": order.upper(),
                    "pdbx_aromatic_flag": aromatic_flag.upper(),
                })
            continue

        # _chem_comp_atom loop
        needed_atom = {
            "_chem_comp_atom.comp_id",
            "_chem_comp_atom.atom_id",
        }
        if needed_atom.issubset(tagset):
            idx = {tag: j for j, tag in enumerate(tags)}
            for row in rows:
                comp_id = cif_unquote(row[idx["_chem_comp_atom.comp_id"]]).upper()
                atom_id = normalize_atom_name(cif_unquote(row[idx["_chem_comp_atom.atom_id"]]))
                type_symbol = cif_unquote(row[idx["_chem_comp_atom.type_symbol"]]) if "_chem_comp_atom.type_symbol" in idx else ""
                charge = cif_unquote(row[idx["_chem_comp_atom.charge"]]) if "_chem_comp_atom.charge" in idx else ""
                chem_comp_atoms[comp_id][atom_id] = {
                    "type_symbol": type_symbol,
                    "charge": charge,
                }

    return chem_comp_bonds, chem_comp_atoms


def residue_to_pdb_block(res: gemmi.Residue, chain_id: str) -> str:
    st = gemmi.Structure()
    st.name = "one_res"
    st.cell = gemmi.UnitCell(1, 1, 1, 90, 90, 90)
    st.spacegroup_hm = "P 1"

    model = gemmi.Model("1")
    chain = gemmi.Chain(chain_id)

    r2 = gemmi.Residue()
    r2.name = res.name
    r2.seqid = gemmi.SeqId(res.seqid.num, res.seqid.icode)

    for a in res:
        a2 = gemmi.Atom()
        a2.name = a.name
        a2.element = a.element
        a2.pos = a.pos
        a2.occ = a.occ
        a2.b_iso = a.b_iso
        r2.add_atom(a2)

    chain.add_residue(r2)
    model.add_chain(chain)
    st.add_model(model)
    return st.make_pdb_string()


def bond_type_from_cif(value_order: str, aromatic_flag: str):
    vo = (value_order or "").strip().upper()
    af = (aromatic_flag or "").strip().upper()

    if vo == "AROM" or af == "Y":
        return Chem.BondType.AROMATIC, True
    if vo in {"SING", "SINGLE"}:
        return Chem.BondType.SINGLE, False
    if vo in {"DOUB", "DOUBLE"}:
        return Chem.BondType.DOUBLE, False
    if vo in {"TRIP", "TRIPLE"}:
        return Chem.BondType.TRIPLE, False
    if vo in {"QUAD", "QUADRUPLE"}:
        return Chem.BondType.QUADRUPLE, False
    # fallback: unknown order, at least keep connectivity
    return Chem.BondType.SINGLE, False


def safe_int(x, default=0):
    try:
        return int(x)
    except Exception:
        return default


def build_rdkit_mol_from_template(res: gemmi.Residue,
                                  chain_id: str,
                                  chem_comp_bonds,
                                  chem_comp_atoms) -> Chem.Mol:
    """
    Build an RDKit molecule for a residue using mmCIF _chem_comp_bond
    and _chem_comp_atom template data.
    """
    comp_id = (res.name or "").strip().upper()
    bond_rows = chem_comp_bonds.get(comp_id, [])
    atom_rows = chem_comp_atoms.get(comp_id, {})

    if not bond_rows:
        raise RuntimeError(f"No _chem_comp_bond template found for component {comp_id}")

    rw = Chem.RWMol()
    atom_name_to_idx = {}

    # Add atoms in residue order, preserving PDB residue info.
    for atom in res:
        atom_name = normalize_atom_name(atom.name)
        template = atom_rows.get(atom_name, {})

        # Prefer _chem_comp_atom.type_symbol when present, else coordinates file atom element.
        type_symbol = (template.get("type_symbol") or "").strip()
        if type_symbol:
            # mmCIF may use "CL", "BR", etc.; RDKit expects normal symbol case
            symbol = type_symbol.capitalize()
            if len(type_symbol) > 1:
                symbol = type_symbol[0].upper() + type_symbol[1:].lower()
        else:
            symbol = atom.element.name.capitalize()

        rd_atom = Chem.Atom(symbol)
        rd_atom.SetNoImplicit(True)

        if "charge" in template and str(template["charge"]).strip() not in {"", "?", "."}:
            rd_atom.SetFormalCharge(safe_int(template["charge"], 0))

        info = Chem.AtomPDBResidueInfo()
        info.SetName(atom.name)
        info.SetResidueName(res.name)
        info.SetChainId(chain_id)
        info.SetResidueNumber(res.seqid.num)
        try:
            info.SetInsertionCode(res.seqid.icode or "")
        except Exception:
            pass
        rd_atom.SetMonomerInfo(info)

        idx = rw.AddAtom(rd_atom)
        atom_name_to_idx[atom_name] = idx

    aromatic_atom_indices = set()

    # Add bonds from _chem_comp_bond.
    for b in bond_rows:
        a1 = normalize_atom_name(b["atom_id_1"])
        a2 = normalize_atom_name(b["atom_id_2"])

        # Skip bonds to atoms absent in this residue instance.
        if a1 not in atom_name_to_idx or a2 not in atom_name_to_idx:
            continue

        i = atom_name_to_idx[a1]
        j = atom_name_to_idx[a2]

        if rw.GetBondBetweenAtoms(i, j) is not None:
            continue

        bond_type, aromatic = bond_type_from_cif(b.get("value_order", ""), b.get("pdbx_aromatic_flag", ""))
        rw.AddBond(i, j, bond_type)

        bond = rw.GetBondBetweenAtoms(i, j)
        if aromatic:
            bond.SetIsAromatic(True)
            aromatic_atom_indices.add(i)
            aromatic_atom_indices.add(j)

    mol = rw.GetMol()

    for idx in aromatic_atom_indices:
        mol.GetAtomWithIdx(idx).SetIsAromatic(True)

    # Try to compute caches/symmetry-friendly properties.
    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    try:
        Chem.GetSymmSSSR(mol)
    except Exception:
        pass

    # Avoid forcing full sanitization when template chemistry is slightly nonstandard.
    try:
        Chem.SanitizeMol(
            mol,
            sanitizeOps=(
                Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                Chem.SanitizeFlags.SANITIZE_SETAROMATICITY |
                Chem.SanitizeFlags.SANITIZE_SETCONJUGATION |
                Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION |
                Chem.SanitizeFlags.SANITIZE_SYMMRINGS
            ),
        )
    except Exception:
        pass

    return Chem.RemoveHs(mol)


def rdkit_mol_from_residue_by_coordinates(res: gemmi.Residue, chain_id: str) -> Chem.Mol:
    """
    Fallback path: infer bonding from coordinates only.
    """
    pdb_block = residue_to_pdb_block(res, chain_id)
    mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    if mol is None:
        raise RuntimeError("RDKit failed to parse residue")

    try:
        rdDetermineBonds.DetermineBonds(mol)
    except Exception:
        rdDetermineBonds.DetermineConnectivity(mol)

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    return Chem.RemoveHs(mol)


def rdkit_mol_from_residue(res: gemmi.Residue,
                           chain_id: str,
                           chem_comp_bonds,
                           chem_comp_atoms) -> Chem.Mol:
    """
    Prefer mmCIF component-template bonding; fall back to coordinate-based inference.
    """
    try:
        return build_rdkit_mol_from_template(res, chain_id, chem_comp_bonds, chem_comp_atoms)
    except Exception:
        return rdkit_mol_from_residue_by_coordinates(res, chain_id)


def equivalence_classes(mol: Chem.Mol):
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    uniq = sorted(set(ranks))
    remap = {r: i for i, r in enumerate(uniq)}
    return [remap[r] for r in ranks]


def iter_ligand_residue_instances(structure: gemmi.Structure):
    for model in structure:
        for chain in model:
            chain_id = chain.name.strip() or "A"
            for res in chain:
                resname = (res.name or "").strip()
                if not resname:
                    continue
                if is_water(resname):
                    continue
                if is_polymer_resname(resname):
                    continue
                natoms = len(res)
                if natoms < 2:
                    continue
                yield resname, natoms, chain_id, res


def select_ligands_max_atoms_per_resname(structure: gemmi.Structure):
    best = {}
    for resname, natoms, chain_id, res in iter_ligand_residue_instances(structure):
        if resname not in best or natoms > best[resname][0]:
            best[resname] = (natoms, chain_id, res)
    return [(resname, v[1], v[2]) for resname, v in best.items()]


def main():
    if len(sys.argv) != 2:
        print("Usage: python3 group-equivalent-ligand-atoms.py input.cif", file=sys.stderr)
        sys.exit(2)

    in_path = sys.argv[1]

    structure = gemmi.read_structure(in_path)
    chem_comp_bonds, chem_comp_atoms = load_chemcomp_tables_from_mmcif(in_path)
    selected = select_ligands_max_atoms_per_resname(structure)

    if not selected:
        print("No ligand-like residues found (non-polymer, non-water).", file=sys.stderr)
        return

    for resname, chain_id, res in selected:
        mol = rdkit_mol_from_residue(res, chain_id, chem_comp_bonds, chem_comp_atoms)
        eq = equivalence_classes(mol)

        counts = Counter(eq)
        keep_classes = {c for c, n in counts.items() if n >= 2}
        if not keep_classes:
            continue

        for atom in mol.GetAtoms():
            cls = eq[atom.GetIdx()]
            if cls not in keep_classes:
                continue

            info = atom.GetPDBResidueInfo()
            if info is None:
                continue

            element = atom.GetSymbol()
            eq_label = f"{element}X{cls}"
            print(f"{info.GetResidueName().strip()}\t{info.GetName().strip()}\t{eq_label}")


if __name__ == "__main__":
    main()