#!/bin/bash

set -u

cd "$WORKROOTDIR"

###########################################################################

mkdir -p ./output/analysis

###########################################################################

R --vanilla << 'EOF'
df_pp1pn1rrc_nrm=read.table("./output/run_rb1_pp1_pn1_rc0/results/global_scores.tsv", header=TRUE, stringsAsFactors=FALSE);
df_pp1pn1rrc_nrm$pp1pn1_rrc_nrm_score=df_pp1pn1rrc_nrm$residue_residue_cadscore;
sdf_pp1pn1rrc_nrm=df_pp1pn1rrc_nrm[, c("target", "model", "pp1pn1_rrc_nrm_score")];

df_pp1pn1rrc=read.table("./output/run_rb1_pp1_pn1_rc1e/results/global_scores.tsv", header=TRUE, stringsAsFactors=FALSE);
df_pp1pn1rs=read.table("./output/run_sites_rb1_pp1_pn1_rc1e/results/global_scores.tsv", header=TRUE, stringsAsFactors=FALSE);
df_pp1pn1rrc$pp1pn1_rrc_score=df_pp1pn1rrc$residue_residue_cadscore;
df_pp1pn1rrc$pp1pn1_ccc_score=df_pp1pn1rrc$chain_chain_cadscore;
df_pp1pn1rs$pp1pn1_rs_score=df_pp1pn1rs$residue_sites_cadscore;
sdf_pp1pn1rrc=df_pp1pn1rrc[, c("target", "model", "pp1pn1_rrc_score", "pp1pn1_ccc_score")];
sdf_pp1pn1rs=df_pp1pn1rs[, c("target", "model", "pp1pn1_rs_score")];

df_pp0pn1rrc=read.table("./output/run_rb1_pp0_pn1_rc1e/results/global_scores.tsv", header=TRUE, stringsAsFactors=FALSE);
df_pp0pn1rs=read.table("./output/run_sites_rb1_pp0_pn1_rc1e/results/global_scores.tsv", header=TRUE, stringsAsFactors=FALSE);
df_pp0pn1rrc$pp0pn1_rrc_score=df_pp0pn1rrc$residue_residue_cadscore;
df_pp0pn1rrc$pp0pn1_ccc_score=df_pp0pn1rrc$chain_chain_cadscore;
df_pp0pn1rs$pp0pn1_rs_score=df_pp0pn1rs$residue_sites_cadscore;
sdf_pp0pn1rrc=df_pp0pn1rrc[, c("target", "model", "pp0pn1_rrc_score", "pp0pn1_ccc_score")];
sdf_pp0pn1rs=df_pp0pn1rs[, c("target", "model", "pp0pn1_rs_score")];

df_pp1pn0rrc=read.table("./output/run_rb1_pp1_pn0_rc1e/results/global_scores.tsv", header=TRUE, stringsAsFactors=FALSE);
df_pp1pn0rs=read.table("./output/run_sites_rb1_pp1_pn0_rc1e/results/global_scores.tsv", header=TRUE, stringsAsFactors=FALSE);
df_pp1pn0rrc$pp1pn0_rrc_score=df_pp1pn0rrc$residue_residue_cadscore;
df_pp1pn0rrc$pp1pn0_ccc_score=df_pp1pn0rrc$chain_chain_cadscore;
df_pp1pn0rs$pp1pn0_rs_score=df_pp1pn0rs$residue_sites_cadscore;
sdf_pp1pn0rrc=df_pp1pn0rrc[, c("target", "model", "pp1pn0_rrc_score", "pp1pn0_ccc_score")];
sdf_pp1pn0rs=df_pp1pn0rs[, c("target", "model", "pp1pn0_rs_score")];

df=sdf_pp1pn1rrc_nrm;
df=merge(df, sdf_pp1pn1rrc);
df=merge(df, sdf_pp1pn1rs);
df=merge(df, sdf_pp0pn1rrc);
df=merge(df, sdf_pp0pn1rs);
df=merge(df, sdf_pp1pn0rrc);
df=merge(df, sdf_pp1pn0rs);

#df=df[which(df$target=="9DXD-assembly1.cif"),]

options(scipen=8);

###########################################################################

png("./output/analysis/plot_remapping_effect.png", height=5, width=5, units="in", res=200);
plot(df$pp1pn1_rrc_nrm_score, df$pp1pn1_rrc_score, xlim=c(0, 1), ylim=c(0, 1),
  xlab="CAD-score BEFORE chain remapping", ylab="CAD-score AFTER chain remapping", main="Reference-based residue-residue CAD-score \n before and after automatic chain remapping", cex=0.6);
abline(0, 1, col="red");
dev.off();

###########################################################################

png("./output/analysis/plot_rr_vs_cc_scores.png", height=10, width=10, units="in", res=200);

par(mfrow=c(2, 2));

plot(df$pp1pn1_rrc_score, df$pp1pn1_ccc_score, xlim=c(0, 1), ylim=c(0, 1),
  xlab="residue-residue CAD-score of all interfaces", ylab="chain-chain CAD-score of all interfaces", main="residue-residue CAD-score vs \n chain-chain CAD-score (all interfaces)", cex=0.6);
abline(0, 1, col="red");
mtext("a", side=3, adj=-0.16, line=1.8, cex=2.0, font=2);

plot(df$pp0pn1_rrc_score, df$pp0pn1_ccc_score, xlim=c(0, 1), ylim=c(0, 1),
  xlab="residue-residue CAD-score of protein-DNA interfaces", ylab="chain-chain CAD-score of protein-DNA interfaces", main="residue-residue CAD-score vs \n chain-chain CAD-score (protein-DNA interfaces)", cex=0.6);
abline(0, 1, col="red");
mtext("b", side=3, adj=-0.16, line=1.8, cex=2.0, font=2);

plot(df$pp1pn0_rrc_score, df$pp1pn0_ccc_score, xlim=c(0, 1), ylim=c(0, 1),
  xlab="residue-residue CAD-score of protein-protein interfaces", ylab="chain-chain CAD-score of protein-protein interfaces", main="residue-residue CAD-score vs \n chain-chain CAD-score (protein-protein interfaces)", cex=0.6);
abline(0, 1, col="red");
mtext("c", side=3, adj=-0.16, line=1.8, cex=2.0, font=2);

plot(df$pp1pn0_ccc_score, df$pp0pn1_ccc_score, xlim=c(0, 1), ylim=c(0, 1),
  xlab="chain-chain CAD-score of protein-protein interfaces", ylab="chain-chain CAD-score of protein-DNA interfaces", main="chain-chain CAD-score of protein-protein interfaces vs \n chain-chain CAD-score of protein-DNA interfaces", cex=0.6);
abline(0, 1, col="red");
mtext("d", side=3, adj=-0.16, line=1.8, cex=2.0, font=2);

dev.off();

###########################################################################

png("./output/analysis/plot_different_interface_types_scores.png", height=5, width=10, units="in", res=200);

par(mfrow=c(1, 2));

plot(df$pp1pn0_rrc_score, df$pp0pn1_rrc_score, xlim=c(0, 1), ylim=c(0, 1),
  xlab="CAD-score of protein-protein interfaces", ylab="CAD-score of protein-DNA interfaces", main="CAD-score of protein-protein interfaces \n vs CAD-score of protein-DNA interfaces", cex=0.6);
abline(0, 1, col="red");
mtext("a", side=3, adj=-0.16, line=1.8, cex=2.0, font=2);

plot(df$pp1pn1_rrc_score, df$pp0pn1_rrc_score, xlim=c(0, 1), ylim=c(0, 1),
  xlab="CAD-score of all interfaces", ylab="CAD-score of protein-DNA interfaces", main="CAD-score of all interfaces \n vs CAD-score of protein-DNA interfaces", cex=0.6);
abline(0, 1, col="red");
mtext("b", side=3, adj=-0.16, line=1.8, cex=2.0, font=2);

dev.off();

###########################################################################

#df=df[which(df$target=="9DXD-assembly1.cif"),]

png("./output/analysis/plot_binding_site_scores.png", height=5, width=10, units="in", res=200);

par(mfrow=c(1, 2));

plot(df$pp1pn0_rrc_score, df$pp1pn0_rs_score, xlim=c(0, 1), ylim=c(0, 1),
  xlab="CAD-score of protein-protein interfaces", ylab="CAD-score of protein binding sites", main="CAD-score of protein-protein interfaces \n vs CAD-score of protein binding sites", cex=0.6);
abline(0, 1, col="red");
mtext("a", side=3, adj=-0.16, line=1.8, cex=2.0, font=2);

plot(df$pp0pn1_rrc_score, df$pp0pn1_rs_score, xlim=c(0, 1), ylim=c(0, 1),
  xlab="CAD-score of protein-DNA interfaces", ylab="CAD-score of DNA binding sites", main="CAD-score of protein-DNA interfaces \n vs CAD-score of DNA binding sites", cex=0.6);
abline(0, 1, col="red");
mtext("b", side=3, adj=-0.16, line=1.8, cex=2.0, font=2);

dev.off();

###########################################################################

EOF

exit 0

