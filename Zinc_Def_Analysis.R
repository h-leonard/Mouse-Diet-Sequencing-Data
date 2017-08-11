## Zinc Deficient


source("https://bioconductor.org/biocLite.R")
library(dada2)
library(ggplot2)
library(DESeq2)
library(phyloseq)
library(randomForest)
library(dplyr)
library(knitr)



sample_df_zinc <- data.frame(sample=samples.out, description=diet_description, day=days_on_diet)
rownames(sample_df_zinc) <- samples.out



# remove anitbiotic time point
sample_df_zinc <- sample_df_zinc[sample_df_zinc$description != "dZD, 23 days, post_Abx",]


sample_df_zinc <- sample_df_zinc[!grepl("dPD", sample_df_zinc$description),]

sample_df_zinc <- sample_df_zinc[!grepl("dN", sample_df_zinc$description),]





ps_zinc <- phyloseq(otu_table(seqtab.nochim.taxa, taxa_are_rows=FALSE), 
                sample_data(sample_df_zinc), 
                tax_table(taxa.table))



# Filter sequences that have less than 10 counts in any sample
ps.count_zinc <- filter_taxa(ps_zinc, function(x) max(x) > 10, TRUE)
ps.count_zinc


#convert to factor for DESeq
sample_data(ps.count_zinc)$description <- factor(sample_data(ps.count_zinc)$description)






# Day 0 comparison ############
deseq2_input_d0 <- phyloseq_to_deseq2(ps.count_zinc, ~description)
deseq2_output_d0 <- DESeq(deseq2_input_d0, test="Wald", fitType="parametric")
deseq2_results_d0 <- results(deseq2_output_d0, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d0 = deseq2_results_d0[which(deseq2_results_d0$padj < alpha), ]
sigtab_d0 = cbind(as(sigtab_d0, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d0), ], "matrix"))
sigtab_d0





seqs <- rownames(sigtab_d0)


# Transform sequence table from counts to relative abundance
ps.rel <- transform_sample_counts(ps.count_zinc, function(x) x/sum(x))
ps.rel

p <- plot_heatmap(prune_taxa(seqs,ps.rel),"NMDS","bray",sample.label = "description","Genus", first.sample = "Plate3-A3")
p + theme (axis.text.x = element_text(size=6))





 