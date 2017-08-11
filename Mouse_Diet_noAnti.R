# Analysis with the Antibiotic Time Points Removed

source("https://bioconductor.org/biocLite.R")
library(dada2)
library(ggplot2)
library(DESeq2)
library(phyloseq)
library(randomForest)
library(dplyr)
library(knitr)



samdf <- data.frame(sample=samples.out, description=diet_description, day=days_on_diet)
rownames(samdf) <- samples.out



# remove antibiotic time point

samdf <- samdf[samdf$description != "dZD, 23 days, post_Abx",]
samdf <- samdf[samdf$description != "dN, 23 days, post_Abx",]
samdf <- samdf[samdf$description != "dPD, 23 days, post_Abx",]






# Construct phyloseq object 
ps <- phyloseq(otu_table(seqtab.nochim.taxa, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa.table))
ps





# Filter sequences that have less than 10 counts in any sample
ps.count <- filter_taxa(ps, function(x) max(x) > 10, TRUE)
ps.count


#convert day values to factor for DESeq
sample_data(ps.count)$day <- factor(sample_data(ps.count)$day)

# Create separate object for comparisons at each time point
d0_subset <- prune_samples(sample_data(ps.count)$day == 0, ps.count)
d5_subset <- prune_samples(sample_data(ps.count)$day == 5, ps.count)
d8_subset <- prune_samples(sample_data(ps.count)$day == 8, ps.count)
d12_subset <- prune_samples(sample_data(ps.count)$day == 12, ps.count)
d15_subset <- prune_samples(sample_data(ps.count)$day == 15, ps.count)




# Day 0 comparison ############
deseq2_input_d0 <- phyloseq_to_deseq2(d0_subset,~description)
deseq2_output_d0 <- DESeq(deseq2_input_d0, test="Wald", fitType="parametric")
deseq2_results_d0 <- results(deseq2_output_d0, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d0 = deseq2_results_d0[which(deseq2_results_d0$padj < alpha), ]
sigtab_d0

#No results less than specified alpha, indicates no initial differences in groups



# Day 5 comparison
deseq2_input_d5 <- phyloseq_to_deseq2(d5_subset,~description)
deseq2_output_d5 <- DESeq(deseq2_input_d5, test="Wald", fitType="parametric")
deseq2_results_d5 <- results(deseq2_output_d5, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d5 = deseq2_results_d5[which(deseq2_results_d5$padj < alpha), ]
sigtab_d5 = cbind(as(sigtab_d5, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d5), ], "matrix"))
head(sigtab_d5)
sigtab_d5

# Day 8 comparison
deseq2_input_d8 <- phyloseq_to_deseq2(d8_subset,~description)
deseq2_output_d8 <- DESeq(deseq2_input_d8, test="Wald", fitType="parametric")
deseq2_results_d8 <- results(deseq2_output_d8, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d8 = deseq2_results_d8[which(deseq2_results_d8$padj < alpha), ]
sigtab_d8 = cbind(as(sigtab_d8, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d8), ], "matrix"))
head(sigtab_d8)


# Day 12 comparison
deseq2_input_d12 <- phyloseq_to_deseq2(d12_subset,~description)
deseq2_output_d12 <- DESeq(deseq2_input_d12, test="Wald", fitType="parametric")
deseq2_results_d12 <- results(deseq2_output_d12, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d12 = deseq2_results_d12[which(deseq2_results_d12$padj < alpha), ]
sigtab_d12 = cbind(as(sigtab_d12, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d12), ], "matrix"))
head(sigtab_d12)


# Day 15 comparison
deseq2_input_d15 <- phyloseq_to_deseq2(d15_subset,~description)
deseq2_output_d15 <- DESeq(deseq2_input_d15, test="Wald", fitType="parametric")
deseq2_results_d15 <- results(deseq2_output_d15, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d15 = deseq2_results_d15[which(deseq2_results_d15$padj < alpha), ]
sigtab_d15 = cbind(as(sigtab_d15, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d15), ], "matrix"))
head(sigtab_d15)




# Get sequences that were significantly different in each comparison
seqs <- union(rownames(sigtab_d5),rownames(sigtab_d8))
seqs <- union(seqs,rownames(sigtab_d12))
seqs <- union(seqs,rownames(sigtab_d15))



# Transform sequence table from counts to relative abundance
ps.rel <- transform_sample_counts(ps.count, function(x) x/sum(x))
ps.rel

p <- plot_heatmap(prune_taxa(seqs,ps.rel),"NMDS","bray",sample.label = "description","Genus", first.sample = "Plate3-A3")
p + theme (axis.text.x = element_text(size=6))



# Random Forest Modeling

# make dataframe of training data from OTU table with samples as rows and OTUs as columns

ps.rel_pruned <- prune_taxa(seqs,ps.rel)


predictors <- (otu_table(ps.rel_pruned))

dim(predictors)

# make response variable column

response <- as.factor(sample_data(ps.rel_pruned)$description)
otu_table(ps.rel_pruned)
tax_table(ps.rel_pruned)
# combine into one frame

rf_data <- data.frame(response, predictors)


# model






set.seed(1000)
diet.classify <- randomForest(response~., data = rf_data, ntree = 214)
print(diet.classify)
plot(diet.classify)


names(diet.classify)

# predictor names and importance


imp <- importance(diet.classify)
imp <- data.frame(predictors = rownames(imp), imp)

# order by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# plot top 20 most important 
imp.20 <- imp.sort[1:20, ]

ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs")




# match names in tax table and otu table
otunames <- imp.20$predictors
otunames <- gsub("X", "",  otunames)


r <- rownames(tax_table(ps.rel)) %in% otunames



k <- as.data.frame((tax_table(ps.rel)[r, ]))

t <- k[match(otunames, rownames(k)),]

kable(t)






