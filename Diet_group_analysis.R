## Comparison of only each diet group


source("https://bioconductor.org/biocLite.R")
library(dada2)
library(ggplot2)
library(DESeq2)
library(phyloseq)
library(randomForest)
library(dplyr)
library(knitr)


sample_df <- data.frame(sample=samples.out, description=diet_description, day=days_on_diet)
rownames(sample_df) <- samples.out



# remove antibiotic time point

sample_df <- sample_df[sample_df$description != "dZD, 23 days, post_Abx",]
sample_df <- sample_df[sample_df$description != "dN, 23 days, post_Abx",]
sample_df <- sample_df[sample_df$description != "dPD, 23 days, post_Abx",]





sample_df$description <- sub(",.*", "", sample_df$description)

ps2 <- phyloseq(otu_table(seqtab.nochim.taxa, taxa_are_rows=FALSE), 
               sample_data(sample_df), 
               tax_table(taxa.table))

ps2

# Filter sequences that have less than 10 counts in any sample
ps.count2 <- filter_taxa(ps2, function(x) max(x) > 10, TRUE)
ps.count2


#convert diets to factor for DESeq
sample_data(ps.count2)$description <- factor(sample_data(ps.count2)$description)
ps.count2



# Create separate object for comparisons of each diet
dn_subset <- prune_samples(sample_data(ps.count2)$description == "dN", ps.count2)
dPD_subset <- prune_samples(sample_data(ps.count2)$description == "dPD", ps.count2)
dZD_subset <- prune_samples(sample_data(ps.count2)$description == "dZD", ps.count2)




# DN comparison ############
deseq2_input_dn <- phyloseq_to_deseq2(dn_subset, ~1)
deseq2_output_dn <- DESeq(deseq2_input_dn, test="Wald", fitType="parametric")
deseq2_results_dn <- results(deseq2_output_dn, cooksCutoff = FALSE)

alpha = 0.05
sigtab_dn = deseq2_results_dn[which(deseq2_results_dn$padj < alpha), ]
sigtab_dn = cbind(as(sigtab_dn, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_dn), ], "matrix"))
sigtab_dn



# DPD comparison
deseq2_input_dpd <- phyloseq_to_deseq2(dPD_subset,~1)
deseq2_output_dpd <- DESeq(deseq2_input_dpd, test="Wald", fitType="parametric")
deseq2_results_dpd <- results(deseq2_output_dpd, cooksCutoff = FALSE)

alpha = 0.05
sigtab_dpd = deseq2_results_dpd[which(deseq2_results_dpd$padj < alpha), ]
sigtab_dpd = cbind(as(sigtab_dpd, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_dpd), ], "matrix"))
head(sigtab_dpd)
sigtab_dpd




# DZD comparison
deseq2_input_dzd <- phyloseq_to_deseq2(dZD_subset,~1)
deseq2_output_dzd <- DESeq(deseq2_input_dzd, test="Wald", fitType="parametric")
deseq2_results_dzd <- results(deseq2_output_dzd, cooksCutoff = FALSE)

alpha = 0.05
sigtab_dzd = deseq2_results_dzd[which(deseq2_results_dzd$padj < alpha), ]
sigtab_dzd = cbind(as(sigtab_dzd, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_dzd), ], "matrix"))
head(sigtab_dzd)






# Get sequences that were significantly different in each comparison
seqs <- union(rownames(sigtab_dn),rownames(sigtab_dpd))
seqs <- union(seqs,rownames(sigtab_dzd))



# Transform sequence table from counts to relative abundance
ps.rel <- transform_sample_counts(ps.count2, function(x) x/sum(x))
ps.rel

p <- plot_heatmap(prune_taxa(seqs,ps.rel),"NMDS","bray",sample.label = "description","Genus", first.sample = "Plate3-A3")
p + theme (axis.text.x = element_text(size=6))





# Random Forest Modeling

# Transform sequence table from counts to relative abundance
ps.rel2 <- transform_sample_counts(ps.count2, function(x) x/sum(x))
ps.rel2


# make dataframe of training data from OTU table with samples as rows and OTUs as columns

predictors2 <- (otu_table(ps.rel2))

dim(predictors2)

p<- as.data.frame(predictors2)

# make response variable column

response2 <- as.factor(sample_data(ps.rel2)$description)


# combine into one frame

rf_data2 <- data.frame(response2, predictors2)

# model

set.seed(1000)
diet.classify2 <- randomForest(response2~., data = rf_data2, ntree = 59)
print(diet.classify2)
plot(diet.classify2)

# Gradient Boosted Model

c <- xgboost(data = data.matrix(p), label = response2, max.depth = 400, nrounds = 55, verboose = 1)
print(c)
importancec <- xgb.importance(feature_names = colnames(p), model = c)


# predictor names and importance


imp2 <- importance(diet.classify2)
imp2 <- data.frame(predictors = rownames(imp2), imp2)

# order by importance
imp.sort2 <- arrange(imp2, desc(MeanDecreaseGini))
imp.sort2$predictors <- factor(imp.sort2$predictors, levels = imp.sort2$predictors)

# plot top 20 most important 
imp.202 <- imp.sort2[1:20, ]

ggplot(imp.202, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs")




# match names in tax table and otu table
otunames <- imp.202$predictors
otunames <- gsub("X", "",  otunames)


r <- rownames(tax_table(ps.rel2)) %in% otunames



k <- as.data.frame((tax_table(ps.rel2)[r, ]))

t <- k[match(otunames, rownames(k)),]

kable(t)

