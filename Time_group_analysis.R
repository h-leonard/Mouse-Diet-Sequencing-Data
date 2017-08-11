## Random Forest model of only the time points 


source("https://bioconductor.org/biocLite.R")
library(dada2)
library(ggplot2)
library(DESeq2)
library(phyloseq)
library(randomForest)
library(dplyr)
library(knitr)




sample_df2 <- data.frame(sample=samples.out, description=diet_description, day=days_on_diet)
rownames(sample_df2) <- samples.out



# remove antibiotic time point
sample_df2 <- sample_df2[sample_df2$description != "dZD, 23 days, post_Abx",]
sample_df2 <- sample_df2[sample_df2$description != "dN, 23 days, post_Abx",]
sample_df2 <- sample_df2[sample_df2$description != "dPD, 23 days, post_Abx",]


sample_df2$description <- sub("^([^ ]*)\\s", "", sample_df2$description)

ps3 <- phyloseq(otu_table(seqtab.nochim.taxa, taxa_are_rows=FALSE), 
                sample_data(sample_df2), 
                tax_table(taxa.table))

# Filter sequences that have less than 10 counts in any sample
ps.count3 <- filter_taxa(ps3, function(x) max(x) > 10, TRUE)
ps.count3



#convert days to factor for DESeq
sample_data(ps.count3)$description <- factor(sample_data(ps.count3)$description)






# Random Forest Modeling

# Transform sequence table from counts to relative abundance
ps.rel3 <- transform_sample_counts(ps.count3, function(x) x/sum(x))
ps.rel3



# make dataframe of training data from OTU table with samples as rows and OTUs as columns

predictors3 <- (otu_table(ps.rel3))

dim(predictors3)


# make response variable column

response3 <- as.factor(sample_data(ps.rel3)$description)



# combine into one frame

rf_data3 <- data.frame(response3, predictors3)

# model

set.seed(1000)
diet.classify3 <- randomForest(response3~., data = rf_data3, ntree = 800)
print(diet.classify3)
plot(diet.classify3)


# predictor names and importance


imp3 <- importance(diet.classify3)
imp3 <- data.frame(predictors = rownames(imp3), imp3)

# order by importance
imp.sort3 <- arrange(imp3, desc(MeanDecreaseGini))
imp.sort3$predictors <- factor(imp.sort3$predictors, levels = imp.sort3$predictors)

# plot top 20 most important 
imp.203 <- imp.sort3[1:20, ]

ggplot(imp.203, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs")




# match names in tax table and otu table
otunames <- imp.203$predictors
otunames <- gsub("X", "",  otunames)


r <- rownames(tax_table(ps.rel3)) %in% otunames

k <- as.data.frame((tax_table(ps.rel3)[r, ]))

t <- k[match(otunames, rownames(k)),]

kable(t)



