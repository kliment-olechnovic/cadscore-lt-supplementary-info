#!/bin/bash

set -u

cd "$WORKROOTDIR"

###########################################################################

mkdir -p ./output/merged_cadscore_tables

find ./output/prepared_input_files/ -type f -name 'rmsd.tsv' | sort \
| xargs cat \
| awk '{if(NR==1 || $1!="model"){print $0}}' \
> ./output/merged_cadscore_tables/rmsd.tsv

find ./output/cadscore/ -type f -name 'global_scores_raw.tsv' | sort \
| xargs cat \
| awk '{if(NR==1 || $1!="target"){print $0}}' \
> ./output/merged_cadscore_tables/global_scores_raw.tsv

find ./output/cadscore/ -type f -name 'global_scores_with_atom_types_info.tsv' | sort \
| xargs cat \
| awk '{if(NR==1 || $1!="target"){print $0}}' \
> ./output/merged_cadscore_tables/global_scores_with_atom_types_info.tsv

cd ./output/merged_cadscore_tables

R --vanilla << 'EOF'
df1=read.table("global_scores_raw.tsv", header=TRUE, stringsAsFactors=FALSE);
nrow(df1);

targets=union(df1$target, df1$target);
target_maxes=c();
for(target in targets)
{
	target_maxes=c(target_maxes, max(df1$residue_residue_cadscore[which(df1$target==target)]));
}
targets=targets[which(target_maxes>0.0)];
df1=df1[which(is.element(df1$target, targets)),];
nrow(df1);

df2=read.table("global_scores_with_atom_types_info.tsv", header=TRUE, stringsAsFactors=FALSE);
nrow(df2);

df3=read.table("rmsd.tsv", header=TRUE, stringsAsFactors=FALSE);
nrow(df3);

df=merge(df1, df2, by=c("target", "model"));
nrow(df);

df=merge(df, df3, by=c("model"));
nrow(df);

df$raw_atom_atom_cadscore=df$atom_atom_cadscore.x;
df$raw_atom_atom_F1_of_areas=df$atom_atom_F1_of_areas.x
df$raw_atom_sites_cadscore=df$atom_sites_cadscore.x;
df$raw_atom_sites_F1_of_areas=df$atom_sites_F1_of_areas.x
df$adjusted_atom_atom_cadscore=df$atom_atom_cadscore.y;
df$adjusted_atom_atom_F1_of_areas=df$atom_atom_F1_of_areas.y
df$adjusted_atom_sites_cadscore=df$atom_sites_cadscore.y;
df$adjusted_atom_sites_F1_of_areas=df$atom_sites_F1_of_areas.y

df=df[,c("target", "model", "raw_atom_atom_cadscore", "raw_atom_atom_F1_of_areas", "raw_atom_sites_cadscore", "raw_atom_sites_F1_of_areas", "adjusted_atom_atom_cadscore", "adjusted_atom_atom_F1_of_areas", "adjusted_atom_sites_cadscore", "adjusted_atom_sites_F1_of_areas", "rmsd")];

write.table(df, file="casf2016_scores.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t");

target_names=sort(union(df$target, df$target));

target_cors_raw_atom_atom_cadscore_vs_rmsd=c();
target_cors_adjusted_atom_atom_cadscore_vs_rmsd=c();

target_cors_raw_atom_sites_cadscore_vs_rmsd=c();
target_cors_adjusted_atom_sites_cadscore_vs_rmsd=c();

target_cors_raw_atom_atom_areas_f1_vs_rmsd=c();
target_cors_adjusted_atom_atom_areas_f1_vs_rmsd=c();

target_rmsd_of_top_atom_atom_cadscore=c();

for(tname in target_names)
{
	sdf=df[which(df$target==tname),];

	target_cors_raw_atom_atom_cadscore_vs_rmsd=c(target_cors_raw_atom_atom_cadscore_vs_rmsd, cor(sdf$rmsd, sdf$raw_atom_atom_cadscore, method="spearman"));
	target_cors_adjusted_atom_atom_cadscore_vs_rmsd=c(target_cors_adjusted_atom_atom_cadscore_vs_rmsd, cor(sdf$rmsd, sdf$adjusted_atom_atom_cadscore, method="spearman"));

	target_cors_raw_atom_sites_cadscore_vs_rmsd=c(target_cors_raw_atom_sites_cadscore_vs_rmsd, cor(sdf$rmsd, sdf$raw_atom_sites_cadscore, method="spearman"));
	target_cors_adjusted_atom_sites_cadscore_vs_rmsd=c(target_cors_adjusted_atom_sites_cadscore_vs_rmsd, cor(sdf$rmsd, sdf$adjusted_atom_sites_cadscore, method="spearman"));

	target_cors_raw_atom_atom_areas_f1_vs_rmsd=c(target_cors_raw_atom_atom_areas_f1_vs_rmsd, cor(sdf$rmsd, sdf$raw_atom_atom_F1_of_areas, method="spearman"));
	target_cors_adjusted_atom_atom_areas_f1_vs_rmsd=c(target_cors_adjusted_atom_atom_areas_f1_vs_rmsd, cor(sdf$rmsd, sdf$adjusted_atom_atom_F1_of_areas, method="spearman"));

	target_rmsd_of_top_atom_atom_cadscore=c(target_rmsd_of_top_atom_atom_cadscore, sdf$rmsd[order(0-sdf$adjusted_atom_atom_cadscore)[1]]);
}

mean(target_cors_raw_atom_atom_cadscore_vs_rmsd);
mean(target_cors_adjusted_atom_atom_cadscore_vs_rmsd);

mean(target_cors_raw_atom_sites_cadscore_vs_rmsd);
mean(target_cors_adjusted_atom_sites_cadscore_vs_rmsd);

mean(target_cors_raw_atom_atom_areas_f1_vs_rmsd);
mean(target_cors_adjusted_atom_atom_areas_f1_vs_rmsd);

target_names[order(0-abs(target_cors_adjusted_atom_atom_cadscore_vs_rmsd))[1:5]];

interesting_target=target_names[order(0-(target_rmsd_of_top_atom_atom_cadscore))[1]];
sdf=df[which(df$target==interesting_target),];
intersting_taget_sel=c(order(sdf$rmsd)[1], order(0-sdf$adjusted_atom_atom_cadscore)[1]);
sdf[intersting_taget_sel, c("target", "model", "adjusted_atom_atom_cadscore", "rmsd", "adjusted_atom_atom_F1_of_areas")];

###############################################################

png("CASF2016_raw_atom_atom_cadscores_vs_atom_atom_cadscores.png", height=5, width=5, units="in", res=200);
plot(df$raw_atom_atom_cadscore, df$adjusted_atom_atom_cadscore,
  xlab="raw atom-atom CAD-score", ylab="symmetry-adjusted atom-atom CAD-score", main="symmetry-adjusted atom-atom CAD-score \n vs raw atom-atom CAD-score",
  col=densCols(df$raw_atom_atom_cadscore, df$adjusted_atom_atom_cadscore));
abline(0, 1);
dev.off();

###############################################################

png("CASF2016_atom_atom_cadscores_vs_atom_sites_cadscores.png", height=5, width=5, units="in", res=200);
plot(df$adjusted_atom_atom_cadscore, df$adjusted_atom_sites_cadscore,
  xlab="symmetry-adjusted atom-atom CAD-score", ylab="atom-sites CAD-score", main="symmetry-adjusted atom-atom CAD-score \n vs atom-sites CAD-score",
  col=densCols(df$adjusted_atom_atom_cadscore, df$adjusted_atom_sites_cadscore));
abline(0, 1);
dev.off();

###############################################################

png("CASF2016_contacts_cadscores_vs_rmsd.png", height=5, width=10, units="in", res=200);

par(mfrow=c(1, 2));

plot(df$adjusted_atom_atom_cadscore, df$rmsd,
  xlab="symmetry-adjusted atom-atom CAD-score", ylab="symmetry-adjusted RMSD", main="symmetry-adjusted atom-atom CAD-score \n vs symmetry-adjusted RMSD",
  col=densCols(df$adjusted_atom_atom_cadscore, df$rmsd));
mtext("a", side=3, adj=-0.16, line=1.8, cex=2.0, font=2);

plot(df$adjusted_atom_sites_cadscore, df$rmsd, xlab="atom sites CAD-score",
  ylab="symmetry-adjusted RMSD", main="atom-sites CAD-score \n vs symmetry-adjusted RMSD",
  col=densCols(df$adjusted_atom_sites_cadscore, df$rmsd));
mtext("b", side=3, adj=-0.16, line=1.8, cex=2.0, font=2);

dev.off();

###############################################################

png("CASF2016_hist_of_spearman_cor_of_atom_atom_cadscores_vs_rmsd.png", height=10, width=7, units="in", res=200);
par(mfrow=c(2, 1));

hist(target_cors_raw_atom_atom_cadscore_vs_rmsd, xlim=c(-1.0, 0.0), breaks=seq(-1, 0, 0.05),
  xlab="raw atom-atom CAD-score", main="Histogram of per-target Spearman cor. coefficients \n between raw atom-atom CAD-score and RMSD");
mtext("a", side=3, adj=-0.11, line=1.8, cex=2.0, font=2);

hist(target_cors_adjusted_atom_atom_cadscore_vs_rmsd, xlim=c(-1.0, 0.0), breaks=seq(-1, 0, 0.05),
  xlab="symmetry-adjusted atom-atom CAD-score", main="Histogram of per-target Spearman cor. coefficients \n between symmetry-adjusted atom-atom CAD-score and RMSD");
mtext("b", side=3, adj=-0.11, line=1.8, cex=2.0, font=2);

dev.off();

###############################################################

png("CASF2016_hist_of_spearman_cor_of_atom_sites_cadscores_vs_rmsd.png", height=10, width=7, units="in", res=200);
par(mfrow=c(2, 1));

hist(target_cors_raw_atom_sites_cadscore_vs_rmsd, xlim=c(-1.0, 0.0), breaks=seq(-1, 0, 0.05),
  xlab="raw atom-atom CAD-score", main="Histogram of per-target Spearman cor. coefficients \n between raw atom sites CAD-score and RMSD");
mtext("a", side=3, adj=-0.11, line=1.8, cex=2.0, font=2);

hist(target_cors_adjusted_atom_sites_cadscore_vs_rmsd, xlim=c(-1.0, 0.0), breaks=seq(-1, 0, 0.05),
  xlab="symmetry-adjusted atom-atom CAD-score", main="Histogram of per-target Spearman cor. coefficients \n between symmetry-adjusted atom sites CAD-score and RMSD");
mtext("b", side=3, adj=-0.11, line=1.8, cex=2.0, font=2);

dev.off();

###############################################################

png("CASF2016_atom_atom_cadscores_vs_atom_atom_f1_of_areas.png", height=5, width=10, units="in", res=200);
par(mfrow=c(1, 2));

plot(df$adjusted_atom_atom_cadscore, df$adjusted_atom_atom_F1_of_areas,
  xlab="symmetry-adjusted atom-atom CAD-score", ylab="symmetry-adjusted atom-atom F1 of areas",
  main="symmetry-adjusted atom-atom CAD-score \n vs symmetry-adjusted atom-atom F1 of areas", col=densCols(df$adjusted_atom_atom_cadscore, df$adjusted_atom_atom_F1_of_areas));
abline(0, 1);
mtext("a", side=3, adj=-0.16, line=1.8, cex=2.0, font=2);

plot(df$adjusted_atom_sites_cadscore, df$adjusted_atom_sites_F1_of_areas,
  xlab="symmetry-adjusted atom sites CAD-score", ylab="symmetry-adjusted atom sites F1 of areas",
  main="symmetry-adjusted atom sites CAD-score \n vs symmetry-adjusted atom sites F1 of areas", col=densCols(df$adjusted_atom_sites_cadscore, df$adjusted_atom_sites_F1_of_areas));
abline(0, 1);
mtext("b", side=3, adj=-0.16, line=1.8, cex=2.0, font=2);

dev.off();

###############################################################

EOF

exit 0

