#### Microbiome analysis for biginner by fleuriane DJONDJI 

#### This protocol should be adapted by users according to what they are looking for

#### last modification: 26/04/2023

#### change the working directory to the analysis folder where all your data is located

setwd("E:/projet/These/16S_sequencing/16s_crid/Microbiome_analyses/data_analysis/Analyse_1")

#### loading packages
library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("qiime2R")
library("grid")
library("phyloseq")
library("devtools")
library("ggplot2")
library("dplyr")
library("tibble")
library("microbiome")
library("vegan")
library("DESeq2")
library("knitr")
library("ape")
library("ggpubr")
library("dendextend")
library("extrafont")

#### import data
otu <- read.csv("OTU_Table_rev.csv", header = T, row.names = 1, sep = ";")
tax <- read.csv("Taxonomy_Table_fun.csv", header = T, row.names = 1, sep = ";")
samples_df <- read.csv("Sample_Table_rev.csv", header = T, row.names = 1, sep = ";")

#### Let's have a look at the different tables:
otu[1:2, 1:6]
tax[1:2, ]
samples_df[1:2, ]

#### Transform into matrixes otu and tax tables (sample table can be left as data frame)
otu_mat <- as.matrix(otu)
tax_mat <- as.matrix(tax)

#### build the phyloseq object
mach <- phyloseq(otu_table(otu_mat, taxa_are_rows = TRUE),
                 tax_table(tax_mat),
                 sample_data(samples_df))

#### Visualize the phyloseq object 
mach
sample_names(mach)
rank_names(mach)
sample_variables(mach)

#### Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(mach))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(mach, standf)
carbom

#### Taxonomic filtering to get only otus present in more than 3 sample
## Create table, number of features for each Genus
table(tax_table(carbom)[, "Genus"], exclude = NULL)


## ensures that features with ambiguous Genus annotation are also removed
ps0 <- subset_taxa(carbom, !is.na(Genus) & !Genus %in% c("", "uncharacterized"))

## Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

## Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))

##Compute the total and average prevalences of the features in each genus.
plyr::ddply(prevdf, "Genus", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

## Rugosibacter appeared in just 2 samples total so we need to Define genus to filter
filterPhyla = c("Rugosibacter")
ps1 = subset_taxa(ps0, !Genus %in% filterPhyla)
ps1

## Subset to the remaining Genus
prevdf1 = subset(prevdf, Genus %in% get_taxa_unique(ps1, "Genus"))

##  Define prevalence threshold of at least 3 sample on the total samples
prevalenceThreshold = 0.142857142857142 * nsamples(ps1)
prevalenceThreshold

## Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps0)
ps2

##compile all the data in one file and exporte it in excel sheet
pseq <- ps2
meta <- meta(pseq)
taxonomy <- tax_table(pseq)
otu.absolute <- abundances(pseq)
otu.relative <- abundances(pseq, "compositional")
reads_sample <- readcount(pseq)
reads_sample[1:5]
sample_data(pseq)$reads_sample <- reads_sample
df <- psmelt(pseq)
kable(head(df))

write.csv2(df, file = "metafleuriane.csv")

#### you can go back into taxonomy and otu file and delete otus which have been removed in meta2
##file then sart to import data again or continuing with ps2 data
#### rarefaction curve
col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
rarecurve(t(otu_table(ps2)), "otu", col= col, step=50, cex=0.5, xlim=c(0, 500000))

## Bacterial abundance
library(phyloseq)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)

# 1. Calcul des abondances relatives
ps_rel_abund <- transform_sample_counts(ps2, function(x) x / sum(x))

# 2. Agglomération au rang "Genus"
ps_genus <- tax_glom(ps_rel_abund, taxrank = "Genus")

# 3. Extraction au format data.frame
df_melt <- psmelt(ps_genus)
df_melt$Genus <- as.character(df_melt$Genus)

# 4. Identification du top 15
top10_genera <- df_melt %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

# 5. Regroupement des autres en "Others"
df_melt$Genus[!(df_melt$Genus %in% top10_genera)] <- "Others"

# 6. Palette de 16 couleurs (10 + Others)
palette_11 <- c(
  "#1b9e77", "#d95f02", "#e7298a", "#7570b3", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#00AFBB", "#FC4E07", "#CCCCCC"  # gris pour "Other"
)

# 7. Tracer le graphique
ggplot(df_melt, aes(x = Phenotypes, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill", width = 0.9) +
  facet_wrap(~Species, scales = "free_x", nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = palette_11) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(
    title = "Microbial composition: Top 10 Bacterial genera",
    x = NULL,
    y = "Relative abundance (%)",
    fill = "Genus"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.y = element_line(color = "gray80"),
    panel.spacing = unit(1.2, "lines"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )


## Anopheles gambiae
carbomgamb <- subset_samples(ps2, Species =="Anopheles_gambiae")
table(phyloseq::tax_table(carbomgamb)[, "Genus"])
phyloseq::otu_table(carbomgamb)[1:5, 1:5]
ps_rel_abundgamb = phyloseq::transform_sample_counts(carbomgamb, function(x){x / sum(x)})
phyloseq::otu_table(ps_rel_abundgamb)[1:5, 1:5]
phyloseq::plot_bar(ps_rel_abundgamb, fill = "Genus") +
  geom_bar(stat="identity", position="fill", width = 0.2) +
  facet_wrap(~Phenotypes, nrow=1, scales="free_x",strip.position ="bottom")+
  theme_bw() +
  scale_y_continuous(labels=scales::percent)+
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90, hjust=1), panel.background = element_rect(fill = 'white', colour = 'white')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance\n") +
  ggtitle("a. Bacterial Composition") +
  labs(fill = "Genus") +
  theme(legend.position = "right") +
  theme(strip.text = element_text(face="bold", size=8,lineheight=5.0))

## Anopheles funestus
carbomfun <- subset_samples(ps2, Species =="Anopheles_funestus")
table(phyloseq::tax_table(carbomfun)[, "Genus"])
phyloseq::otu_table(carbomfun)[1:5, 1:5]
ps_rel_abundfun = phyloseq::transform_sample_counts(carbomfun, function(x){x / sum(x)})
phyloseq::otu_table(ps_rel_abundfun)[1:5, 1:5]
phyloseq::plot_bar(ps_rel_abundfun, fill = "Genus") +
  geom_bar(stat="identity", position="fill", width = 0.2) +
  facet_wrap(~Phenotypes, nrow=1, scales="free_x",strip.position ="bottom")+
  theme_bw() +
  scale_y_continuous(labels=scales::percent)+
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90, hjust=1), panel.background = element_rect(fill = 'white', colour = 'white')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance\n") +
  ggtitle("Bacterial Composition") +
  labs(fill = "Genus") +
  theme(legend.position = "right") +
  theme(strip.text = element_text(face="bold", size=8,lineheight=5.0))

#### Cumulative Relative abundance of Genus
## Per Specie

Top <- names(sort(taxa_sums(ps2), decreasing=TRUE))
dat.aglo = tax_glom(ps2, taxrank = "Genus")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
prune.dat.two = prune_taxa(Top, dat.trans)
dat.dataframe = psmelt(prune.dat.two)
dat.agr = aggregate(Abundance~Species+Genus, data=dat.dataframe, FUN=mean)
ggplot(dat.agr, aes(x=Species, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity", position="fill") +
  theme_bw() +
  scale_y_continuous(labels=scales::percent) +
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.x = element_blank(), panel.background = element_rect(fill = 'white', colour = 'white')) + 
  ylab("Relative Abundance\n") +
  ggtitle("Bacterial Composition") +
  labs(fill = "Genus") +
  theme(legend.position = "right") +
  theme(strip.text = element_text(face="bold", size=8,lineheight=5.0))

## An. gambiae
Topgamb<- names(sort(taxa_sums(carbomgamb), decreasing=TRUE))
dat.aglo.gamb = tax_glom(carbomgamb, taxrank = "Genus")
dat.trans.gamb= transform_sample_counts(dat.aglo.gamb, function(x) x/sum(x))
prune.dat.two.gamb = prune_taxa(Topgamb, dat.trans.gamb)
dat.dataframe.gamb = psmelt(prune.dat.two.gamb)
dat.agr.gamb = aggregate(Abundance~Phenotypes+Genus, data=dat.dataframe.gamb, FUN=mean)
ggplot(dat.agr.gamb, aes(x=Phenotypes, y=Abundance, fill=Genus))+
  geom_bar(stat="identity", position="fill") +
  theme_bw() +
  scale_y_continuous(labels=scales::percent) +
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.x = element_blank(), panel.background = element_rect(fill = 'white', colour = 'white')) + 
  ylab("Relative Abundance\n") +
  ggtitle("Bacterial Composition") +
  labs(fill = "Genus") +
  theme(legend.position = "right") +
  theme(strip.text = element_text(face="bold", size=8,lineheight=5.0))

## An. funestus
Topfun<- names(sort(taxa_sums(ps2), decreasing=TRUE))
dat.aglo.fun = tax_glom(carbomfun, taxrank = "Genus")
dat.trans.fun= transform_sample_counts(dat.aglo.fun, function(x) x/sum(x))
prune.dat.two.fun = prune_taxa(Topfun, dat.trans.fun)
dat.dataframe.fun = psmelt(prune.dat.two.fun)
dat.agr.fun = aggregate(Abundance~Phenotypes+Genus, data=dat.dataframe.fun, FUN=mean)
ggplot(dat.agr.fun, aes(x=Phenotypes, y=Abundance, fill=Genus))+
  geom_bar(stat="identity", position="fill") +
  theme_bw() +
  scale_y_continuous(labels=scales::percent) +
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.x = element_blank(), panel.background = element_rect(fill = 'white', colour = 'white')) + 
  ylab("Relative Abundance\n") +
  ggtitle("Bacterial Composition") +
  labs(fill = "Genus") +
  theme(legend.position = "right") +
  theme(strip.text = element_text(face="bold", size=8,lineheight=5.0))



ALPHADIVERSITY


####Test whether the observed number of OTUs differs significantly between groups (Alpha diversity parameters)
ps3 <- prune_taxa(taxa_sums(ps2) > 0, ps2)
tab <- microbiome::alpha(ps2, index = "all")
ps3.meta <- meta(ps2)
ps3.meta$Shannon <- tab$diversity_shannon 
ps3.meta$inversesimpson <- tab$diversity_inverse_simpson
hist(ps3.meta$Shannon, main="Shannon diversity", xlab="", breaks=10)

plot_richness(ps2, x="Species", color = "Species", measures=c("Observed", "Shannon", "Chao", "ace")) + 
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="", axis.text.x=element_text(angle=0,hjust=1,vjust=1,size=8))

d <- meta(ps2)
d$diversity <- microbiome::diversity(ps2, "Shannon")$shannon
spl <- split(d$diversity, d$Species)
bmi <- levels(ps3.meta$Species)
bmi.pairs <- combn(seq_along(spl), 2, simplify = FALSE, FUN = function(i)bmi[i])
p <- ggboxplot(ps3.meta, x = "Species", y = c("Shannon"),
               color = "Species", palette =c("#E7B800", "#00AFBB"),
               add = "jitter", shape = "Species")
p + stat_compare_means(comparisons = bmi.pairs) +
  stat_compare_means(label.y = 2) 

## Anopheles gambiae
ps4 <- prune_taxa(taxa_sums(carbomgamb) > 0, carbomgamb)
tab <- microbiome::alpha(carbomgamb, index = "all")
ps4.meta <- meta(carbomgamb)
ps4.meta$Shannon <- tab$diversity_shannon 
ps4.meta$InverseSimpson <- tab$diversity_inverse_simpson
hist(ps4.meta$Shannon, main="Shannon diversity", xlab="", breaks=10)
dg<- meta(carbomgamb)
dg$diversity <- microbiome::diversity(carbomgamb, "Shannon")$Shannon
splg <- split(d$diversity, d$phenotypes)
bmig <- levels(ps4.meta$phenotypes)
bmig.pairs <- combn(seq_along(splg), 2, simplify = FALSE, FUN = function(i)bmig[i])
p4 <- ggboxplot(ps4.meta, x = "Phenotypes", y = "Shannon",
                color = "Phenotypes", palette =c("#E7B800", "#00AFBB", "#00FF00", "#0000FF", "#00FFFF", "#0000A0", "#ADD8E6", "#800080", "#FFFF00", "#00ff00", "#808000", "#008000", "#800000", "#A52A2A"),
                add = "jitter", shape = "Phenotypes")
p4 + stat_compare_means(method = "kruskal.test", label.y = 2) + 
  stat_compare_means(aes(label = after_stat(p.signif)), ref.group = "Alive_1X")

## Anopheles funestus
ps5 <- prune_taxa(taxa_sums(carbomfun) > 0, carbomfun)
tab <- microbiome::alpha(carbomfun, index = "all")
ps5.meta <- meta(carbomfun)
ps5.meta$Shannon <- tab$diversity_shannon 
ps5.meta$InverseSimpson <- tab$diversity_inverse_simpson
hist(ps5.meta$Shannon, main="Shannon diversity", xlab="", breaks=10)
dfu<- meta(carbomfun)
dfu$diversity <- microbiome::diversity(carbomfun, "Shannon")$Shannon
splfu <- split(d$diversity, d$phenotypes)
bmifu <- levels(ps5.meta$phenotypes)
bmifu.pairs <- combn(seq_along(splfu), 2, simplify = FALSE, FUN = function(i)bmifu[i])
p5 <- ggboxplot(ps5.meta, x = "Phenotypes", y = "Shannon",
                color = "Phenotypes", palette =c("#E7B800", "#00AFBB"),
                add = "jitter", shape = "Phenotypes")
p5 + stat_compare_means(method = "wilcox.test", label.y = 2)+ 
  stat_compare_means(aes(label = after_stat(p.signif)), ref.group = "Fumoz_Unselect")


####PHYLOGENETIC TREE

####We can create a random phylogenetic tree with the ape package                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           dataset. Make sure its tip labels match your OTU_table.
##Species
random_tree = rtree(ntaxa(ps2), rooted=TRUE, tip.label=taxa_names(ps2))
plot(random_tree)
physeq = merge_phyloseq(ps2, random_tree)
physeq
plot_tree(physeq, color = "Genus", label.tips="Genus", ladderize="right", plot.margin=0.3)

##Anopheles gambiae
random_treegamb = rtree(ntaxa(carbomgamb), rooted=TRUE, tip.label=taxa_names(carbomgamb))
plot(random_treegamb)
physeqgamb = merge_phyloseq(carbomgamb, random_treegamb)
physeqgamb
plot_tree(physeqgamb, color = "Genus", label.tips="Genus", ladderize="right", plot.margin=0.3)

##Anopheles funestus
random_treefun = rtree(ntaxa(carbomfun), rooted=TRUE, tip.label=taxa_names(carbomfun))
plot(random_treefun)
physeqfun = merge_phyloseq(carbomfun, random_treefun)
physeqfun
plot_tree(physeqfun, color = "Genus", label.tips="Genus", ladderize="right", plot.margin=0.3)

##Hierarchical clustering based on Bray-Curtis Index values, showing the relationship between different samples and group
##acccording to Anopheles species
ps_rel_abundgamb = phyloseq::transform_sample_counts(physeqgamb, function(x){x / sum(x)})
phyloseq::otu_table(physeqgamb)[1:5, 1:5]
ps_rel_otugamb <- data.frame(phyloseq::otu_table(ps_rel_abundgamb))
ps_rel_otugamb <- t(ps_rel_otugamb)
bc_distgamb <- vegan::vegdist(ps_rel_otugamb, method = "bray")
as.matrix(bc_distgamb)[1:5, 1:5]
wardgamb <- as.dendrogram(hclust(bc_distgamb, method = "ward.D2"))
metagamb <- data.frame(phyloseq::sample_data(ps_rel_abundgamb))
colorCode <- c(`Gpp` = "#0048BA", `Gtach` = "#B0BF1A", `Gcal` = "#7CB9E8", `Gpal` = "#C0E8D5")
labels_colors(wardgamb) <- colorCode[meta$Phenotypes][order.dendrogram(wardgamb)]
plot(wardgamb)

#We can create a random phylogenetic tree with the ape package in An. funestus                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           dataset. Make sure its tip labels match your OTU_table.
random_treefun = rtree(ntaxa(carbomfun), rooted=TRUE, tip.label=taxa_names(carbomfun))
plot(random_treefun)
physeqfun = merge_phyloseq(carbomfun, random_treefun)
physeqfun
plot_tree(physeqfun, color = "Genus", label.tips="Genus", ladderize="left", plot.margin=0.3)
plot_tree(physeqfun, color="Phenotypes", shape="Phenotypes", label.tips="Genus", ladderize="left", plot.margin=0.3)

##Hierarchical clustering based on Bray-Curtis Index values, showing the relationship between different samples and group
##acccording to Anopheles species
ps_rel_abundfun = phyloseq::transform_sample_counts(physeqfun, function(x){x / sum(x)})
phyloseq::otu_table(physeqfun)[1:5, 1:5]
ps_rel_otufun <- data.frame(phyloseq::otu_table(ps_rel_abundfun))
ps_rel_otufun <- t(ps_rel_otufun)
bc_distfun <- vegan::vegdist(ps_rel_otufun, method = "bray")
as.matrix(bc_distfun)[1:5, 1:5]
wardfun <- as.dendrogram(hclust(bc_distfun, method = "ward.D2"))
metafun <- data.frame(phyloseq::sample_data(ps_rel_abundfun))
colorCode <- c(`Gpp` = "#0048BA", `Gtach` = "#B0BF1A", `Gcal` = "#7CB9E8", `Gpal` = "#C0E8D5")
labels_colors(wardfun) <- colorCode[meta$Phenotypes][order.dendrogram(wardfun)]
plot(wardfun)



BETADIVERSITY

library("phyloseq")
library("ggplot2")
library("dplyr")
library("ggpubr")
library("Matrix")
library("reshape2")
library("vegan")

## Compare community dissimilarities in each group using bray-curtis index
relab_genera = transform_sample_counts(ps2, function(x) x / sum(x) * 100) 
ord = ordinate(relab_genera, method="PCoA", distance = "bray")
plot_ordination(relab_genera, ord, color = "Species", shape="Test") + 
  geom_point(size=3) + 
  stat_ellipse(aes(Species=Species))

##PERMANOVA significance test for group-level differences using bray-curtis index
pseq.rel <- microbiome::transform(ps2, "compositional")
otu <- abundances(ps2)
meta <- meta(pseq.rel)
permanova <- adonis2(t(otu) ~ Species,
                     data = meta, permutations=999, method = "bray")
permanova

## Compare community dissimilarities in each phenotype using bray-curtis index
##Anopheles gambiae
relab_genera1 = transform_sample_counts(carbomgamb, function(x) x / sum(x) * 100) 
ord1 = ordinate(relab_genera1, method="PCoA", distance = "bray")
plot_ordination(relab_genera1, ord1, color = "Phenotypes", shape="Test") + 
  geom_point(size=3) + 
  stat_ellipse(aes(Phenotypes=Phenotypes))
otu1 <- abundances(carbomgamb)
pseq.rel1 <- microbiome::transform(carbomgamb, "compositional")
meta1 <- meta(pseq.rel1)
permanova1 <- adonis2(t(otu1) ~ Phenotypes, 
                      data = meta1, permutations=999, method = "bray")
permanova1

##Anopheles funestus
relab_genera2 = transform_sample_counts(carbomfun, function(x) x / sum(x) * 100) 
ord2 = ordinate(relab_genera2, method="PCoA", distance = "bray")
plot_ordination(relab_genera2, ord2, color = "Phenotypes", shape="Test") + 
  geom_point(size=3) + 
  stat_ellipse(aes(Phenotypes=Phenotypes))
otu2 <- abundances(carbomfun)
pseq.rel2 <- microbiome::transform(carbomfun, "compositional")
meta2 <- meta(pseq.rel2)
permanova2 <- adonis2(t(otu2) ~ Phenotypes,
                      data = meta2, permutations=999, method = "bray")
permanova2






ds2 <- phyloseq_to_deseq2(carbomgamb, ~ Phenotypes)
dds <- DESeq(ds2)
res <- results(dds)
deseq.results <- as.data.frame(res)
df <- deseq.results
df$taxon <- rownames(df)
df1 <- df %>% arrange(log2FoldChange, padj)
df2 <- df1 %>% filter(pvalue < 0.05 & log2FoldChange < -1.5) %>%
  arrange(pvalue, log2FoldChange)
df3 <- df1 %>% filter(pvalue < 0.05 & log2FoldChange > 1.5) %>%
  arrange(pvalue, log2FoldChange)
df3
kable(df2, digits = 5)
kable(df3, digits = 5)















library(lefser)
res <- lefser(ps2, groupCol = "codes", blockCol = "codes")
sample_data(ps2)$code <- as.factor(sample_data(ps2)$code)
ps.taxa <- tax_glom(ps2, taxrank = 'Genus', NArm = FALSE)
ps.taxa.pse <- ps2
otu_table(ps.taxa.pse) <- otu_table(ps.taxa) + 1
ps.taxa.pse.sub <- subset_samples(ps.taxa.pse, code %in% c("fzu", "fzse"))
ds = phyloseq_to_deseq2(ps.taxa.pse.sub, ~ code)
ds = DESeq(ds, test="Wald", fitType="non parametric")
ds
alpha = 0.05 
res = results(ds, alpha=alpha)
res = res[order(res$padj, na.last=NA), ]

res

library(BiocManager)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "4.3")




















#### Plot Heatmap abundance plot of the Genus to capture the best candidate bacteria
carbom_abund <- filter_taxa(ps2, function(x) sum(x > total*0.05) > 0, TRUE)
otu_table(carbom_abund)
plot_heatmap(carbom_abund, method = "MDS", distance = "bray", 
             taxa.label = "Genus", taxa.order = "Genus", 
             trans=NULL, low="beige", high="red", na.value="beige")

## Anopheles gambiae 
carbomgamb_abund <- filter_taxa(carbom_fraction, function(x) sum(x > total*0.05) > 0, TRUE)
otu_table(carbomgamb_abund)
plot_heatmap(carbomgamb_abund, method = "MDS", distance = "bray", 
             taxa.label = "Genus", taxa.order = "Genus", 
             trans=NULL, low="beige", high="red", na.value="beige")

## Anopheles funestus
carbom_abundfun <- filter_taxa(carbom_fractionfun, function(x) sum(x > total*0.05) > 0, TRUE)
otu_table(carbom_abundfun)
plot_heatmap(carbom_abundfun, method = "MDS", distance = "bray", 
             taxa.label = "Genus", taxa.order = "Genus", 
             trans=NULL, low="beige", high="red", na.value="beige")

