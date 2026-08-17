#!/bin/bash

set -u

################################################################################

mkdir -p "$OUTPUT_SUMMARY_DIR"

################################################################################

R --vanilla --args "$OUTPUT_DATA_DIR" "$OUTPUT_SUMMARY_DIR" << 'EOF'
args=commandArgs(TRUE);
dir_data=args[1];
dir_summary=args[2];

v1=scan(paste0(dir_data, "/extremely_compact/raw/results/residue_residue_cadscore_global.tsv"));
v2=scan(paste0(dir_data, "/extremely_compact/renumbered/results/residue_residue_cadscore_global.tsv"));
v3=scan(paste0(dir_data, "/extremely_compact/renumbered_remapped_clustered/results/residue_residue_cadscore_global.tsv"));

options(scipen=8);

png(paste0(dir_summary, "/histograms_of_scores_by_preparation_level.png"), width=6, height=12, units="in", res=200);
par(mfrow=c(3, 1));
hist(v1/100, breaks=seq(0, 1, 0.05), xlim=c(0, 1), xlab="residue-residue CAD-score",
  main="Residue-residue CAD-scores of models without preprocessing");
hist(v2/100, breaks=seq(0, 1, 0.05), xlim=c(0, 1), xlab="residue-residue CAD-score",
  main="Residue-residue CAD-scores of models after renumbering residues \n and assigning chain names based on refrence sequences");
hist(v3/100, breaks=seq(0, 1, 0.05), xlim=c(0, 1), xlab="residue-residue CAD-score",
  main="Residue-residue CAD-scores of models after renumbering residues \n and assigning chain names, and remapping chain names");
dev.off();

EOF

################################################################################

R --vanilla --args "$OUTPUT_DATA_DIR" "$OUTPUT_SUMMARY_DIR" << 'EOF'
args=commandArgs(TRUE);
dir_data=args[1];
dir_summary=args[2];

df1=read.table(paste0(dir_data, "/normal/reference_based_renumbered/results/global_scores.tsv"), header=TRUE, stringsAsFactors=FALSE);
df2=read.table(paste0(dir_data, "/normal/reference_based_renumbered_remapped/results/global_scores.tsv"), header=TRUE, stringsAsFactors=FALSE);
df=merge(df1, df2, by=c("target", "model"));
nrow(df);

options(scipen=8);

png(paste0(dir_summary, "/plot_reference_based_remapping_effect.png"), height=5, width=5, units="in", res=200);

plot(df$residue_residue_cadscore.x, df$residue_residue_cadscore.y,
  xlab="CAD-score BEFORE chain remapping", ylab="CAD-score AFTER chain remapping", main="Reference-based residue-residue CAD-score \n before and after automatic chain remapping", cex=0.6);
abline(0, 1, col="red");

dev.off();

EOF

################################################################################

exit 0

