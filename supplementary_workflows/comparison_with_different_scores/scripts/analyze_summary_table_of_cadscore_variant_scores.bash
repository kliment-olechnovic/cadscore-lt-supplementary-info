#!/bin/bash

set -u

################################################################################

cd "$OUTPUT_SUMMARY_DIR"

R --vanilla > ./log_of_analysis_of_cadscore_variant_scores.txt << 'EOF'
df=read.table("summary_table_of_cadscore_variant_scores.tsv", header=TRUE, stringsAsFactors=FALSE);
nrow(df);

# RNA models from the 208 and 304 groups
# often had extreme chain overlaps (e.g. R1283v3TS208_1o and R1283v3TS304_1o)
# therefore they were excluded from the analysis
df=df[grep("^R.*TS208", df$input, invert=TRUE),];
df=df[grep("^R.*TS304", df$input, invert=TRUE),];
nrow(df);

df=df[which(df$cv1_model_atoms==df$cv3_model_atoms & df$cv1_model_atoms==df$cv4_model_atoms & df$cv1_target_atoms==df$cv3_target_atoms & df$cv1_target_atoms==df$cv4_target_atoms),];
nrow(df);

df=df[which((df$cv1_model_atoms/df$cv1_model_residues)>6 & (df$cv1_target_atoms/df$cv1_target_residues)>6),];
nrow(df);

#######################################

df$target_atoms=df$cv1_target_atoms;
df$model_atoms=df$cv1_model_atoms;

df$total_input_atoms=(df$target_atoms+df$model_atoms);

df$CADscore2012_global_score=df$old_global_score;
df$CADscore2021_global_score=df$voronota_js_global_score;
df$CADscoreLT_global_score=df$cadscore_lt_global_score;
df$CADscoreLT_remapped_global_score=df$cadscore_lt_rm_global_score;
df$CADscoreLT_refseq_global_score=df$cadscore_lt_r_global_score;

df$CADscore2012_iface_score=df$old_iface_score;
df$CADscore2021_iface_score=df$voronota_js_iface_score;
df$CADscoreLT_iface_score=df$cadscore_lt_iface_score;
df$CADscoreLT_remapped_iface_score=df$cadscore_lt_rm_iface_score;
df$CADscoreLT_refseq_iface_score=df$cadscore_lt_r_iface_score;

df$CADscore2012_global_time=df$old_global_time;
df$CADscore2021_global_time=df$voronota_js_global_time;
df$CADscoreLT_global_time=df$cadscore_lt_global_time;
df$CADscoreLT_remapped_global_time=df$cadscore_lt_rm_global_time;
df$CADscoreLT_refseq_global_time=df$cadscore_lt_r_global_time;
df$CADscore2021_LT_global_time=df$voronota_js_lt_global_time;

df$CADscore2012_iface_time=df$old_iface_time;
df$CADscore2021_iface_time=df$voronota_js_iface_time;
df$CADscoreLT_iface_time=df$cadscore_lt_iface_time;
df$CADscoreLT_remapped_iface_time=df$cadscore_lt_rm_iface_time;
df$CADscoreLT_refseq_iface_time=df$cadscore_lt_r_iface_time;
df$CADscore2021_LT_iface_time=df$voronota_js_lt_iface_time;

main_colnames=grep("^CADscore", colnames(df), value=TRUE);

df=df[,c("input", "target_atoms", "model_atoms", "total_input_atoms", main_colnames)];

write.table(df, file="comparison_of_cadscore_versions.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t");

name_CADscore2012="CAD-score-2012";
name_CADscore2021="CAD-score-2021";
name_CADscoreLT="CAD-score-LT";

name_CADscore2012_global=paste0(name_CADscore2012, " whole structure score");
name_CADscore2021_global=paste0(name_CADscore2021, " whole structure score");
name_CADscoreLT_global=paste0(name_CADscoreLT, " whole structure score");
name_CADscoreLT_remapped_global=paste0(name_CADscoreLT, " whole structure score after remapping");

name_CADscore2012_iface=paste0(name_CADscore2012, " inter-chain interface score");
name_CADscore2021_iface=paste0(name_CADscore2021, " inter-chain interface score");
name_CADscoreLT_iface=paste0(name_CADscoreLT, " inter-chain interface score");
name_CADscoreLT_remapped_iface=paste0(name_CADscoreLT, " inter-chain interface score after remapping");

#######################################

options(scipen=8);

#######################################

png("./plot_comparing_CADscoreLT_vs_CADscore2012_scores.png", height=4.5, width=9.0, units="in", res=200);
par(mfrow=c(1, 2));

plot(df$CADscore2012_global_score, df$CADscoreLT_global_score, cex=0.6,
  xlab=name_CADscore2012_global, ylab=name_CADscoreLT_global, main=paste0("Whole structure scores\n", name_CADscoreLT, " vs ", name_CADscore2012));
abline(0, 1, col="red");
mtext("a", side=3, adj=-0.17, line=2.5, cex=2.0, font=2);

plot(df$CADscore2012_iface_score, df$CADscoreLT_iface_score, cex=0.6,
  xlab=name_CADscore2012_iface, ylab=name_CADscoreLT_iface, main=paste0("Inter chain interface scores\n", name_CADscoreLT, " vs ", name_CADscore2012));
abline(0, 1, col="red");
mtext("b", side=3, adj=-0.17, line=2.5, cex=2.0, font=2);

dev.off();

#######################################

png("./plot_comparing_CADscoreLT_vs_CADscore2021_scores.png", height=4.5, width=9.0, units="in", res=200);
par(mfrow=c(1, 2));

plot(df$CADscore2021_global_score, df$CADscoreLT_global_score, cex=0.6,
  xlab=name_CADscore2021_global, ylab=name_CADscoreLT_global, main=paste0("Whole structure scores\n", name_CADscoreLT, " vs ", name_CADscore2021));
abline(0, 1, col="red");
mtext("a", side=3, adj=-0.17, line=2.5, cex=2.0, font=2);

plot(df$CADscore2021_iface_score, df$CADscoreLT_iface_score, cex=0.6,
  xlab=name_CADscore2021_iface, ylab=name_CADscoreLT_iface, main=paste0("Inter chain interface scores\n", name_CADscoreLT, " vs ", name_CADscore2021));
abline(0, 1, col="red");
mtext("b", side=3, adj=-0.17, line=2.5, cex=2.0, font=2);

dev.off();

#######################################

col_CADscore2012=palette()[3];
col_CADscore2021=palette()[2];
col_CADscoreLT=palette()[4];
col_CADscore2021modLT=palette()[8];

#######################################

png("./plot_comparing_times_of_CADscore_versions.png", height=7.0, width=7.0, units="in", res=200);
par(mfrow=c(2, 2));

plot(x=c(min(df$total_input_atoms), max(df$total_input_atoms)), y=c(0, max(c(df$CADscore2012_global_time, df$CADscore2021_global_time, df$CADscoreLT_global_time))), type="n",
  xlab="Total input atoms", ylab="Seconds (CPU user time)", main="Whole structure scoring times");
points(x=df$total_input_atoms, y=df$CADscore2012_global_time, col=col_CADscore2012, cex=0.6, pch=2);
points(x=df$total_input_atoms, y=df$CADscore2021_global_time, col=col_CADscore2021, cex=0.6, pch=6);
points(x=df$total_input_atoms, y=df$CADscoreLT_global_time, col=col_CADscoreLT, cex=0.5, pch=1);
legend("topleft", legend=c("CAD-score-2012", "CAD-score-2021", "CAD-score-LT"), col=c(col_CADscore2012, col_CADscore2021, col_CADscoreLT), bty="n", pch=c(2, 6, 1), pt.cex=1.2);
mtext("a", side=3, adj=-0.16, line=1.8, cex=1.5, font=2);

plot(x=c(min(df$total_input_atoms), max(df$total_input_atoms)), y=c(0, 7), type="n",
  xlab="Total input atoms", ylab="Seconds (CPU user time)", main="Whole structure scoring times\n(showing times below 7 seconds)");
points(x=df$total_input_atoms, y=df$CADscore2012_global_time, col=col_CADscore2012, cex=0.6, pch=2);
points(x=df$total_input_atoms, y=df$CADscore2021_global_time, col=col_CADscore2021, cex=0.6, pch=6);
points(x=df$total_input_atoms, y=df$CADscoreLT_global_time, col=col_CADscoreLT, cex=0.5, pch=1);
mtext("b", side=3, adj=-0.16, line=1.8, cex=1.5, font=2);

plot(x=c(min(df$total_input_atoms), max(df$total_input_atoms)), y=c(0, max(c(df$CADscore2012_iface_time, df$CADscore2021_iface_time, df$CADscoreLT_iface_time))), type="n",
  xlab="Total input atoms", ylab="Seconds (CPU user time)", main="Inter-chain interface scoring times");
points(x=df$total_input_atoms, y=df$CADscore2012_iface_time, col=col_CADscore2012, cex=0.6, pch=2);
points(x=df$total_input_atoms, y=df$CADscore2021_iface_time, col=col_CADscore2021, cex=0.6, pch=6);
points(x=df$total_input_atoms, y=df$CADscoreLT_iface_time, col=col_CADscoreLT, cex=0.5, pch=1);
legend("topleft", legend=c("CAD-score-2012", "CAD-score-2021", "CAD-score-LT"), col=c(col_CADscore2012, col_CADscore2021, col_CADscoreLT), bty="n", pch=c(2, 6, 1), pt.cex=1.2);
mtext("c", side=3, adj=-0.16, line=1.8, cex=1.5, font=2);

plot(x=c(min(df$total_input_atoms), max(df$total_input_atoms)), y=c(0, 7), type="n",
  xlab="Total input atoms", ylab="Seconds (CPU user time)", main="Inter-chain interface scoring times\n(showing times below 7 seconds)");
points(x=df$total_input_atoms, y=df$CADscore2012_iface_time, col=col_CADscore2012, cex=0.6);
points(x=df$total_input_atoms, y=df$CADscore2021_iface_time, col=col_CADscore2021, cex=0.6);
points(x=df$total_input_atoms, y=df$CADscoreLT_iface_time, col=col_CADscoreLT, cex=0.5);
mtext("d", side=3, adj=-0.16, line=1.8, cex=1.5, font=2);

dev.off();

#######################################

png("./plot_comparing_times_of_CADscore_LT_and_JS_LT_versions.png", height=5.0, width=10.0, units="in", res=200);
par(mfrow=c(1, 2));

plot(x=c(min(df$total_input_atoms), max(df$total_input_atoms)), y=c(0, max(c(df$CADscore2021_LT_global_time, df$CADscoreLT_global_time))), type="n",
  xlab="Total input atoms", ylab="Seconds (CPU user time)", main="Whole structure scoring times:\n CAD-score-2021-LT-mod vs CAD-score-LT");
points(x=df$total_input_atoms, y=df$CADscore2021_LT_global_time, col=col_CADscore2021modLT, cex=0.6, pch=6);
points(x=df$total_input_atoms, y=df$CADscoreLT_global_time, col=col_CADscoreLT, cex=0.5, pch=1);
legend("topleft", legend=c("CAD-score-2021-LT-mod", "CAD-score-LT"), col=c(col_CADscore2021modLT, col_CADscoreLT), bty="n", pch=c(6, 1), pt.cex=1.2);
mtext("a", side=3, adj=-0.16, line=1.8, cex=1.5, font=2);

plot(x=c(min(df$total_input_atoms), max(df$total_input_atoms)), y=c(0, max(c(df$CADscore2021_LT_iface_time, df$CADscoreLT_iface_time))), type="n",
  xlab="Total input atoms", ylab="Seconds (CPU user time)", main="Inter-chain interface scoring times:\n CAD-score-2021-LT-mod vs CAD-score-LT");
points(x=df$total_input_atoms, y=df$CADscore2021_LT_iface_time, col=col_CADscore2021modLT, cex=0.6, pch=6);
points(x=df$total_input_atoms, y=df$CADscoreLT_iface_time, col=col_CADscoreLT, cex=0.5, pch=1);
legend("topleft", legend=c("CAD-score-2021-LT-mod", "CAD-score-LT"), col=c(col_CADscore2021modLT, col_CADscoreLT), bty="n", pch=c(6, 1), pt.cex=1.2);
mtext("b", side=3, adj=-0.16, line=1.8, cex=1.5, font=2);

dev.off();

#######################################

png("./plot_comparing_CADscore_LT_with_and_without_chain_remapping.png", height=6.0, width=12.0, units="in", res=200);
par(mfrow=c(1, 2));

plot(df$CADscoreLT_refseq_global_score, df$CADscoreLT_remapped_global_score, cex=0.6,
  xlab=name_CADscoreLT_global, ylab=name_CADscoreLT_remapped_global, main=paste0(name_CADscoreLT_global, " vs \n ", name_CADscoreLT_remapped_global));
abline(0, 1, col="red");
mtext("a", side=3, adj=-0.17, line=2.5, cex=2.0, font=2);

plot(df$CADscoreLT_refseq_iface_score, df$CADscoreLT_remapped_iface_score, cex=0.6,
  xlab=name_CADscoreLT_iface, ylab=name_CADscoreLT_remapped_iface, main=paste0(name_CADscoreLT_iface, " vs \n ", name_CADscoreLT_remapped_iface));
abline(0, 1, col="red");
mtext("b", side=3, adj=-0.16, line=1.8, cex=1.5, font=2);

dev.off();

#######################################

png("./plot_comparing_times_of_CADscore_LT_with_and_without_chain_remapping.png", height=6.0, width=12.0, units="in", res=200);
par(mfrow=c(1, 2));

plot(df$CADscoreLT_refseq_global_time, df$CADscoreLT_remapped_global_time, cex=0.6,
  xlab="Seconds (CPU user time) without remapping", ylab="Seconds (CPU user time) with remapping", main="CAD-score-LT whole structure scoring \n times without and with chain remapping");
abline(0, 1, col="red");
mtext("a", side=3, adj=-0.17, line=2.5, cex=2.0, font=2);

plot(df$CADscoreLT_refseq_iface_time, df$CADscoreLT_remapped_iface_time, cex=0.6,
  xlab="Seconds (CPU user time) without remapping", ylab="Seconds (CPU user time) with remapping", main="CAD-score-LT inter-chain interface scoring \n times without and with chain remapping");
abline(0, 1, col="red");
mtext("b", side=3, adj=-0.16, line=1.8, cex=1.5, font=2);

dev.off();

#######################################

quantile(df$total_input_atoms);

cor(df$CADscore2012_global_score, df$CADscoreLT_global_score);
mean(abs(df$CADscore2012_global_score-df$CADscoreLT_global_score));

cor(df$CADscore2012_iface_score, df$CADscoreLT_iface_score);
mean(abs(df$CADscore2012_iface_score-df$CADscoreLT_iface_score));

cor(df$CADscore2021_global_score, df$CADscoreLT_global_score);
mean(abs(df$CADscore2021_global_score-df$CADscoreLT_global_score));

cor(df$CADscore2021_iface_score, df$CADscoreLT_iface_score);
mean(abs(df$CADscore2021_iface_score-df$CADscoreLT_iface_score));

#######################################

sum(df$CADscore2012_global_time)/sum(df$CADscoreLT_global_time);
sum(df$CADscore2012_iface_time)/sum(df$CADscoreLT_iface_time);

sum(df$CADscore2021_global_time)/sum(df$CADscoreLT_global_time);
sum(df$CADscore2021_iface_time)/sum(df$CADscoreLT_iface_time);

sum(df$CADscore2021_LT_global_time)/sum(df$CADscoreLT_global_time);
sum(df$CADscore2021_LT_iface_time)/sum(df$CADscoreLT_iface_time);

sum(df$CADscoreLT_remapped_global_time)/sum(df$CADscoreLT_refseq_global_time);
sum(df$CADscoreLT_remapped_iface_time)/sum(df$CADscoreLT_refseq_iface_time);

df=df[which(df$CADscoreLT_iface_time>0.001),];

quantile(df$CADscore2012_global_time/df$CADscoreLT_global_time);
quantile(df$CADscore2012_iface_time/df$CADscoreLT_iface_time);

quantile(df$CADscore2021_global_time/df$CADscoreLT_global_time);
quantile(df$CADscore2021_iface_time/df$CADscoreLT_iface_time);

quantile(df$CADscore2021_LT_global_time/df$CADscoreLT_global_time);
quantile(df$CADscore2021_LT_iface_time/df$CADscoreLT_iface_time);

quantile(df$CADscoreLT_remapped_global_time/df$CADscoreLT_refseq_global_time);
quantile(df$CADscoreLT_remapped_iface_time/df$CADscoreLT_refseq_iface_time);

#######################################

EOF

################################################################################

exit 0

