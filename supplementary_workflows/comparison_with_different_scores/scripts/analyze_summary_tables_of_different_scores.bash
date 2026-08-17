#!/bin/bash

set -u

################################################################################

cd "$OUTPUT_SUMMARY_DIR"

R --vanilla > ./log_of_analysis_of_different_scores.txt << 'EOF'
df1=read.table("comparison_of_cadscore_versions.tsv", header=TRUE, stringsAsFactors=FALSE);
df1=df1[,c("input", "target_atoms", "model_atoms", "total_input_atoms", "CADscoreLT_global_score", "CADscoreLT_remapped_global_score", "CADscoreLT_refseq_global_score", "CADscoreLT_iface_score", "CADscoreLT_remapped_iface_score", "CADscoreLT_refseq_iface_score", "CADscoreLT_remapped_global_time")];

df2_auto=read.table("summary_table_of_openstructure_auto_scores.tsv", header=TRUE, stringsAsFactors=FALSE);

#######################################

df=merge(df1, df2_auto);
write.table(df, file="comparison_of_cadscore_versions_with_openstructure_scores.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t");

cor(df$ost_ilddt, df$CADscoreLT_remapped_iface_score);
cor(df$ost_qs_global, df$CADscoreLT_remapped_iface_score);
cor(df$ost_qs_best, df$CADscoreLT_remapped_iface_score);
cor(df$ost_ics, df$CADscoreLT_remapped_iface_score);
cor(df$ost_ips, df$CADscoreLT_remapped_iface_score);
cor(df$ost_dockq_ave, df$CADscoreLT_remapped_iface_score);
cor(df$ost_dockq_wave, df$CADscoreLT_remapped_iface_score);
cor(df$ost_dockq_ave_full, df$CADscoreLT_remapped_iface_score);
cor(df$ost_dockq_wave_full, df$CADscoreLT_remapped_iface_score);
cor(df$ost_tm_score, df$CADscoreLT_remapped_iface_score);

png("./plot_comparing_CADscoreLT_vs_OpenStructure_scores_part1.png", height=9.0, width=9.0, units="in", res=200);
par(mfrow=c(2, 2));

vx=df$ost_ilddt;
vy=df$CADscoreLT_remapped_iface_score;
vcol=densCols(vx, vy);
plot(x=vx, y=vy, col=vcol, cex=0.6,
  xlab="ilDDT", ylab="CAD-score-LT inter-chain interface score", main=paste0("CAD-score-LT inter-chain vs ilDDT\nPearson cor. coeff. = ", sprintf("%.2f", cor(vx, vy))));
abline(0, 1, col="red");
mtext("a", side=3, adj=-0.11, line=2.0, cex=2.0, font=2);

vx=df$ost_tm_score;
vy=df$CADscoreLT_remapped_iface_score;
vcol=densCols(vx, vy);
plot(x=vx, y=vy, col=vcol, cex=0.6,
  xlab="TM-score", ylab="CAD-score-LT inter-chain interface score", main=paste0("CAD-score-LT inter-chain vs TM-score\nPearson cor. coeff. = ", sprintf("%.2f", cor(vx, vy))));
abline(0, 1, col="red");
mtext("b", side=3, adj=-0.11, line=2.0, cex=2.0, font=2);

vx=df$ost_ics;
vy=df$CADscoreLT_remapped_iface_score;
vcol=densCols(vx, vy);
plot(x=vx, y=vy, col=vcol, cex=0.6,
  xlab="ICS", ylab="CAD-score-LT inter-chain interface score", main=paste0("CAD-score-LT inter-chain vs ICS\nPearson cor. coeff. = ", sprintf("%.2f", cor(vx, vy))));
abline(0, 1, col="red");
mtext("c", side=3, adj=-0.11, line=2.0, cex=2.0, font=2);

vx=df$ost_ips;
vy=df$CADscoreLT_remapped_iface_score;
vcol=densCols(vx, vy);
plot(x=vx, y=vy, col=vcol, cex=0.6,
  xlab="IPS", ylab="CAD-score-LT inter-chain interface score", main=paste0("CAD-score-LT inter-chain vs IPS\nPearson cor. coeff. = ", sprintf("%.2f", cor(vx, vy))));
abline(0, 1, col="red");
mtext("d", side=3, adj=-0.11, line=2.0, cex=2.0, font=2);

dev.off();

png("./plot_comparing_CADscoreLT_vs_OpenStructure_scores_part2.png", height=9.0, width=9.0, units="in", res=200);
par(mfrow=c(2, 2));

vx=df$ost_qs_global;
vy=df$CADscoreLT_remapped_iface_score;
vcol=densCols(vx, vy);
plot(x=vx, y=vy, col=vcol, cex=0.6,
  xlab="QS-global", ylab="CAD-score-LT inter-chain interface score", main=paste0("CAD-score-LT inter-chain vs QS-global\nPearson cor. coeff. = ", sprintf("%.2f", cor(vx, vy))));
abline(0, 1, col="red");
mtext("e", side=3, adj=-0.11, line=2.0, cex=2.0, font=2);

vx=df$ost_qs_best
vy=df$CADscoreLT_remapped_iface_score;
vcol=densCols(vx, vy);
plot(x=vx, y=vy, col=vcol, cex=0.6,
  xlab="QS-best", ylab="CAD-score-LT inter-chain interface score", main=paste0("CAD-score-LT inter-chain vs QS-best\nPearson cor. coeff. = ", sprintf("%.2f", cor(vx, vy))));
abline(0, 1, col="red");
mtext("f", side=3, adj=-0.11, line=2.0, cex=2.0, font=2);

vx=df$ost_dockq_wave;
vy=df$CADscoreLT_remapped_iface_score;
vcol=densCols(vx, vy);
plot(x=vx, y=vy, col=vcol, cex=0.6,
  xlab="DockQ-wave", ylab="CAD-score-LT inter-chain interface score", main=paste0("CAD-score-LT inter-chain vs DockQ-wave\nPearson cor. coeff. = ", sprintf("%.2f", cor(vx, vy))));
abline(0, 1, col="red");
mtext("g", side=3, adj=-0.11, line=2.0, cex=2.0, font=2);

vx=df$ost_dockq_wave_full;
vy=df$CADscoreLT_remapped_iface_score;
vcol=densCols(vx, vy);
plot(x=vx, y=vy, col=vcol, cex=0.6,
  xlab="DockQ-wave-full", ylab="CAD-score-LT inter-chain interface score", main=paste0("CAD-score-LT inter-chain vs DockQ-wave-full\nPearson cor. coeff. = ", sprintf("%.2f", cor(vx, vy))));
abline(0, 1, col="red");
mtext("h", side=3, adj=-0.11, line=2.0, cex=2.0, font=2);

dev.off();

#######################################

df2_noremap=read.table("summary_table_of_openstructure_scores.tsv", header=TRUE, stringsAsFactors=FALSE);
df2_noremap$ost_ilddt_no_remapping=df2_noremap$ost_ilddt;

df2_mix=merge(df2_auto[, c("input", "ost_ilddt")], df2_noremap[, c("input", "ost_ilddt_no_remapping")]);

df=merge(df1[, c("input", "CADscoreLT_refseq_iface_score", "CADscoreLT_remapped_iface_score")], df2_mix);

write.table(df, file="comparison_of_interchain_cadscore_and_openstructure_ilddt_with_and_without_chain_remapping.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t");

cor(df$ost_ilddt_no_remapping, df$CADscoreLT_refseq_iface_score);
cor(df$ost_ilddt, df$CADscoreLT_refseq_iface_score);
cor(df$ost_ilddt_no_remapping, df$CADscoreLT_remapped_iface_score);
cor(df$ost_ilddt, df$CADscoreLT_remapped_iface_score);

png("./plot_comparing_interchain_CADscoreLT_vs_OpenStructure_ilDDT_with_and_without_chain_remapping.png", height=9.0, width=9.0, units="in", res=200);
par(mfrow=c(2, 2));

vx=df$ost_ilddt_no_remapping;
vy=df$CADscoreLT_refseq_iface_score;
vcol=densCols(vx, vy);
plot(x=vx, y=vy, col=vcol, cex=0.6,
  xlab="ilDDT", ylab="CAD-score-LT inter-chain interface score", main=paste0("CAD-score-LT without chain remapping \n vs ilDDT without chain remapping", "\nPearson cor. coeff. = ", sprintf("%.2f", cor(vx, vy))));
abline(0, 1, col="red");
mtext("a", side=3, adj=-0.14, line=2.0, cex=2.0, font=2);

vx=df$ost_ilddt;
vy=df$CADscoreLT_refseq_iface_score;
vcol=densCols(vx, vy);
plot(x=vx, y=vy, col=vcol, cex=0.6,
  xlab="ilDDT", ylab="CAD-score-LT inter-chain interface score", main=paste0("CAD-score-LT without chain remapping \n vs ilDDT with chain remapping", "\nPearson cor. coeff. = ", sprintf("%.2f", cor(vx, vy))));
abline(0, 1, col="red");
mtext("b", side=3, adj=-0.14, line=2.0, cex=2.0, font=2);

vx=df$ost_ilddt_no_remapping;
vy=df$CADscoreLT_remapped_iface_score;
vcol=densCols(vx, vy);
plot(x=vx, y=vy, col=vcol, cex=0.6,
  xlab="ilDDT", ylab="CAD-score-LT inter-chain interface score", main=paste0("CAD-score-LT with chain remapping \n vs ilDDT without chain remapping", "\nPearson cor. coeff. = ", sprintf("%.2f", cor(vx, vy))));
abline(0, 1, col="red");
mtext("c", side=3, adj=-0.14, line=2.0, cex=2.0, font=2);

vx=df$ost_ilddt;
vy=df$CADscoreLT_remapped_iface_score;
vcol=densCols(vx, vy);
plot(x=vx, y=vy, col=vcol, cex=0.6,
  xlab="ilDDT", ylab="CAD-score-LT inter-chain interface score", main=paste0("CAD-score-LT with chain remapping \n vs ilDDT with chain remapping", "\nPearson cor. coeff. = ", sprintf("%.2f", cor(vx, vy))));
abline(0, 1, col="red");
mtext("d", side=3, adj=-0.14, line=2.0, cex=2.0, font=2);

dev.off();

#######################################

EOF

################################################################################

exit 0

