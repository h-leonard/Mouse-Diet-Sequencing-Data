## Heatmap comparison of each time point separately


source("https://bioconductor.org/biocLite.R")
library(dada2)
library(ggplot2)
library(DESeq2)
library(phyloseq)
library(randomForest)
library(dplyr)
library(knitr)




####################### Diets at day 5



sample_df_day5 <- data.frame(sample=samples.out, description=diet_description, day=days_on_diet)
rownames(sample_df_day5) <- samples.out



# remove antibiotic time point
sample_df_day5 <- sample_df_day5[sample_df_day5$description != "dZD, 23 days, post_Abx",]
sample_df_day5 <- sample_df_day5[sample_df_day5$description != "dN, 23 days, post_Abx",]
sample_df_day5 <- sample_df_day5[sample_df_day5$description != "dPD, 23 days, post_Abx",]



sample_df_day5 <- sample_df_day5[!grepl("weaned", sample_df_day5$description),]
sample_df_day5 <- sample_df_day5[!grepl("8 days", sample_df_day5$description),]
sample_df_day5 <- sample_df_day5[!grepl("12 days", sample_df_day5$description),]
sample_df_day5 <- sample_df_day5[!grepl("15 days", sample_df_day5$description),]



ps_day5 <- phyloseq(otu_table(seqtab.nochim.taxa, taxa_are_rows=FALSE), 
                    sample_data(sample_df_day5), 
                    tax_table(taxa.table))



# Filter sequences that have less than 10 counts in any sample
ps.count_day5 <- filter_taxa(ps_day5, function(x) max(x) > 10, TRUE)
ps.count_day5



#convert to factor for DESeq
sample_data(ps.count_day5)$description <- factor(sample_data(ps.count_day5)$description)



# Day 0 comparison ############
deseq2_input_d0 <- phyloseq_to_deseq2(ps.count_day5, ~description)
deseq2_output_d0 <- DESeq(deseq2_input_d0, test="Wald", fitType="parametric")
deseq2_results_d0 <- results(deseq2_output_d0, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d0 = deseq2_results_d0[which(deseq2_results_d0$padj < alpha), ]
sigtab_d0 = cbind(as(sigtab_d0, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d0), ], "matrix"))
sigtab_d0


seqs <- rownames(sigtab_d0)




# Transform sequence table from counts to relative abundance
ps.rel <- transform_sample_counts(ps.count_day5, function(x) x/sum(x))
ps.rel

p <- plot_heatmap(prune_taxa(seqs,ps.rel),"NMDS","bray",sample.label = "description","Genus", first.sample = "Plate3-A4")
p + theme (axis.text.x = element_text(size=6))








####################### Diets at day 8


sample_df_day8 <- data.frame(sample=samples.out, description=diet_description, day=days_on_diet)
rownames(sample_df_day8) <- samples.out



# remove antibiotic time point
sample_df_day8 <- sample_df_day8[sample_df_day8$description != "dZD, 23 days, post_Abx",]
sample_df_day8 <- sample_df_day8[sample_df_day8$description != "dN, 23 days, post_Abx",]
sample_df_day8 <- sample_df_day8[sample_df_day8$description != "dPD, 23 days, post_Abx",]



sample_df_day8 <- sample_df_day8[!grepl("weaned", sample_df_day8$description),]
sample_df_day8 <- sample_df_day8[!grepl("5 days", sample_df_day8$description),]
sample_df_day8 <- sample_df_day8[!grepl("12 days", sample_df_day8$description),]
sample_df_day8 <- sample_df_day8[!grepl("15 days", sample_df_day8$description),]



ps_day8 <- phyloseq(otu_table(seqtab.nochim.taxa, taxa_are_rows=FALSE), 
                    sample_data(sample_df_day8), 
                    tax_table(taxa.table))



# Filter sequences that have less than 10 counts in any sample
ps.count_day8 <- filter_taxa(ps_day8, function(x) max(x) > 10, TRUE)
ps.count_day8



#convert to factor for DESeq
sample_data(ps.count_day8)$description <- factor(sample_data(ps.count_day8)$description)



# Day 0 comparison ############
deseq2_input_d0 <- phyloseq_to_deseq2(ps.count_day8, ~description)
deseq2_output_d0 <- DESeq(deseq2_input_d0, test="Wald", fitType="parametric")
deseq2_results_d0 <- results(deseq2_output_d0, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d0 = deseq2_results_d0[which(deseq2_results_d0$padj < alpha), ]
sigtab_d0 = cbind(as(sigtab_d0, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d0), ], "matrix"))
sigtab_d0


seqs <- rownames(sigtab_d0)




# Transform sequence table from counts to relative abundance
ps.rel <- transform_sample_counts(ps.count_day8, function(x) x/sum(x))
ps.rel

p <- plot_heatmap(prune_taxa(seqs,ps.rel),"NMDS","bray",sample.label = "description","Genus", first.sample = "Plate3-A4")
p + theme (axis.text.x = element_text(size=6))












####################### Diets on day 12


sample_df_day12 <- data.frame(sample=samples.out, description=diet_description, day=days_on_diet)
rownames(sample_df_day12) <- samples.out



# remove antibiotic time point
sample_df_day12 <- sample_df_day12[sample_df_day12$description != "dZD, 23 days, post_Abx",]
sample_df_day12 <- sample_df_day12[sample_df_day12$description != "dN, 23 days, post_Abx",]
sample_df_day12 <- sample_df_day12[sample_df_day12$description != "dPD, 23 days, post_Abx",]



sample_df_day12 <- sample_df_day12[!grepl("weaned", sample_df_day12$description),]
sample_df_day12 <- sample_df_day12[!grepl("5 days", sample_df_day12$description),]
sample_df_day12 <- sample_df_day12[!grepl("8 days", sample_df_day12$description),]
sample_df_day12 <- sample_df_day12[!grepl("15 days", sample_df_day12$description),]



ps_day12 <- phyloseq(otu_table(seqtab.nochim.taxa, taxa_are_rows=FALSE), 
                     sample_data(sample_df_day12), 
                     tax_table(taxa.table))



# Filter sequences that have less than 10 counts in any sample
ps.count_day12 <- filter_taxa(ps_day12, function(x) max(x) > 10, TRUE)
ps.count_day12



#convert to factor for DESeq
sample_data(ps.count_day12)$description <- factor(sample_data(ps.count_day12)$description)



# Day 0 comparison ############
deseq2_input_d0 <- phyloseq_to_deseq2(ps.count_day12, ~description)
deseq2_output_d0 <- DESeq(deseq2_input_d0, test="Wald", fitType="parametric")
deseq2_results_d0 <- results(deseq2_output_d0, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d0 = deseq2_results_d0[which(deseq2_results_d0$padj < alpha), ]
sigtab_d0 = cbind(as(sigtab_d0, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d0), ], "matrix"))
sigtab_d0


seqs <- rownames(sigtab_d0)




# Transform sequence table from counts to relative abundance
ps.rel <- transform_sample_counts(ps.count_day12, function(x) x/sum(x))
ps.rel

p <- plot_heatmap(prune_taxa(seqs,ps.rel),"NMDS","bray",sample.label = "description","Genus", first.sample = "Plate3-A10")
p + theme (axis.text.x = element_text(size=6))










####################### Diets on day 15


sample_df_day15 <- data.frame(sample=samples.out, description=diet_description, day=days_on_diet)
rownames(sample_df_day15) <- samples.out



# remove antibiotic time point
sample_df_day15 <- sample_df_day15[sample_df_day15$description != "dZD, 23 days, post_Abx",]
sample_df_day15 <- sample_df_day15[sample_df_day15$description != "dN, 23 days, post_Abx",]
sample_df_day15 <- sample_df_day15[sample_df_day15$description != "dPD, 23 days, post_Abx",]



sample_df_day15 <- sample_df_day15[sample_df_day15$day == "15",]



ps_day15 <- phyloseq(otu_table(seqtab.nochim.taxa, taxa_are_rows=FALSE), 
                     sample_data(sample_df_day15), 
                     tax_table(taxa.table))



# Filter sequences that have less than 10 counts in any sample
ps.count_day15 <- filter_taxa(ps_day15, function(x) max(x) > 10, TRUE)
ps.count_day15



#convert to factor for DESeq
sample_data(ps.count_day15)$description <- factor(sample_data(ps.count_day15)$description)



# Day 0 comparison ############
deseq2_input_d0 <- phyloseq_to_deseq2(ps.count_day15, ~description)
deseq2_output_d0 <- DESeq(deseq2_input_d0, test="Wald", fitType="parametric")
deseq2_results_d0 <- results(deseq2_output_d0, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d0 = deseq2_results_d0[which(deseq2_results_d0$padj < alpha), ]
sigtab_d0 = cbind(as(sigtab_d0, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d0), ], "matrix"))
sigtab_d0


seqs <- rownames(sigtab_d0)




# Transform sequence table from counts to relative abundance
ps.rel <- transform_sample_counts(ps.count_day15, function(x) x/sum(x))
ps.rel

p <- plot_heatmap(prune_taxa(seqs,ps.rel),"NMDS","bray",sample.label = "description","Genus", first.sample = "Plate3-H12")
p + theme (axis.text.x = element_text(size=6))











