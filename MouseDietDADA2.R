# Hampton Leonard

source("https://bioconductor.org/biocLite.R")
library(dada2)
biocLite("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2)
biocLite("phyloseq")
biocLite("DESeq2")
library("DESeq2")
library(phyloseq)
biocLite("S4Vectors")
biocLite("data.table")
packageVersion("phyloseq")



path <- "C:/Users/hll4c/R/Mouse Data/Unzipped_final"

fns <- list.files(path)
fns


fastqs <- fns[grepl(".fastq$", fns)]

fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order

fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files


# Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
sample.names


# Fully specify the path for the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)


par(mfrow=c(2,3))
#forward
plotQualityProfile(fnFs[[1]])
plotQualityProfile(fnFs[[5]])
plotQualityProfile(fnFs[[13]])
plotQualityProfile(fnFs[[17]])

#reverse
plotQualityProfile(fnRs[[1]])
plotQualityProfile(fnRs[[5]])
plotQualityProfile(fnRs[[13]])
plotQualityProfile(fnRs[[17]])


# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))


# Filter
for(i in seq_along(fnFs)) {
  fastqFilter(fnFs[i], filtFs[i],
                    trimLeft=20, truncLen=200, 
                    maxN=0, maxEE=2, truncQ=2, 
                    compress=TRUE, verbose=TRUE)
}


derepFs <- derepFastq(filtFs, verbose=TRUE)


# Name the derep-class objects by the sample names
names(derepFs) <- sample.names


derepFs[[1]]


# Sample Inference
dadaFs <- dada(derepFs, err=NULL, selfConsist = TRUE)


#Visualize error rates
plotErrors(dadaFs[[2]], nominalQ=TRUE)


# construct sequence table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
table(nchar(colnames(seqtab)))


# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
#0.8267492


# assign taxonomy
path2 <- "C:/Users/hll4c/R/Mouse Data"

taxa <- assignTaxonomy(seqtab.nochim, paste(path2,"/rdp_train_set_14.fa.gz",sep = ""))
taxa.plus <- addSpecies(taxa, paste(path2,"/rdp_species_assignment_14.fa.gz",sep = ""))
colnames(taxa.plus) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")
unname(taxa.plus)



# Phyloseq analysis

# load metadata for sorting
sample <- read.csv("C:/Users/hll4c/R/SampleMouse.csv")


# Make a data.frame holding the sample data
rownames(sample) <- sample$SAMPLE_ID
samples.out <- rownames(seqtab.nochim)
sample <- sample[samples.out,]
diet_description = sample$Diet_and_des
days_on_diet = sample$Days.on.diet


# save seqtab.nochim and taxa.plus for analysis w/ metagenomeSeq
write.table(t(seqtab.nochim), file="C:/Users/hll4c/R/Mouse Data/seqtab.nochim.tsv", quote=FALSE, sep='\t', col.names = NA)
write.table(taxa.plus, file="C:/Users/hll4c/R/Mouse Data/taxa.plus.tsv", quote=FALSE, sep='\t', col.names = NA)


samdf <- data.frame(sample=samples.out, description=diet_description, day=days_on_diet)
rownames(samdf) <- samples.out


# copy to create a new object that has taxa information
seqtab.nochim.taxa <- seqtab.nochim
colnames(seqtab.nochim.taxa) <- 1:nrow(taxa.plus)
taxa.table <- taxa.plus
#######rownames(taxa.table) <- 1:nrow(taxa.plus)#?

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa.table))
ps

import::from(phyloseq, reconcile_categories)

# Filter sequences that don't have more than 10 counts in any sample
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
d23_subset <- prune_samples(sample_data(ps.count)$day == 23, ps.count)



# D0 comparison ############
deseq2_input_d0 <- phyloseq_to_deseq2(d0_subset,~description)
deseq2_output_d0 <- DESeq(deseq2_input_d0, test="Wald", fitType="parametric")
deseq2_results_d0 <- results(deseq2_output_d0, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d0 = deseq2_results_d0[which(deseq2_results_d0$padj < alpha), ]
sigtab_d0 = cbind(as(sigtab_d0, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d0), ], "matrix"))
head(sigtab_d0)

ps.count
sigtab_d0
deseq2_results_d0$padj


# D5 comparison
deseq2_input_d5 <- phyloseq_to_deseq2(d5_subset,~description)
deseq2_output_d5 <- DESeq(deseq2_input_d5, test="Wald", fitType="parametric")
deseq2_results_d5 <- results(deseq2_output_d5, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d5 = deseq2_results_d5[which(deseq2_results_d5$padj < alpha), ]
sigtab_d5 = cbind(as(sigtab_d5, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d5), ], "matrix"))
head(sigtab_d5)


# D8 comparison
deseq2_input_d8 <- phyloseq_to_deseq2(d8_subset,~description)
deseq2_output_d8 <- DESeq(deseq2_input_d8, test="Wald", fitType="parametric")
deseq2_results_d8 <- results(deseq2_output_d8, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d8 = deseq2_results_d8[which(deseq2_results_d8$padj < alpha), ]
sigtab_d8 = cbind(as(sigtab_d8, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d8), ], "matrix"))
head(sigtab_d8)


# D12 comparison
deseq2_input_d12 <- phyloseq_to_deseq2(d12_subset,~description)
deseq2_output_d12 <- DESeq(deseq2_input_d12, test="Wald", fitType="parametric")
deseq2_results_d12 <- results(deseq2_output_d12, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d12 = deseq2_results_d12[which(deseq2_results_d12$padj < alpha), ]
sigtab_d12 = cbind(as(sigtab_d12, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d12), ], "matrix"))
head(sigtab_d12)


# D15 comparison
deseq2_input_d15 <- phyloseq_to_deseq2(d15_subset,~description)
deseq2_output_d15 <- DESeq(deseq2_input_d15, test="Wald", fitType="parametric")
deseq2_results_d15 <- results(deseq2_output_d15, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d15 = deseq2_results_d15[which(deseq2_results_d15$padj < alpha), ]
sigtab_d15 = cbind(as(sigtab_d15, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d15), ], "matrix"))
head(sigtab_d15)


# D23 comparison
deseq2_input_d23 <- phyloseq_to_deseq2(d23_subset,~description)
deseq2_output_d23 <- DESeq(deseq2_input_d23, test="Wald", fitType="parametric")
deseq2_results_d23 <- results(deseq2_output_d23, cooksCutoff = FALSE)

alpha = 0.05
sigtab_d23 = deseq2_results_d23[which(deseq2_results_d23$padj < alpha), ]
sigtab_d23 = cbind(as(sigtab_d23, "data.frame"), as(tax_table(ps.count)[rownames(sigtab_d23), ], "matrix"))
head(sigtab_d23)


# Get sequences that were significantly different in each comparison
seqs <- union(rownames(sigtab_d5),rownames(sigtab_d8))
seqs <- union(seqs,rownames(sigtab_d12))
seqs <- union(seqs,rownames(sigtab_d15))
seqs <- union(seqs,rownames(sigtab_d23))


# Transform sequence table from counts to relative abundance
ps.rel <- transform_sample_counts(ps.count, function(x) x/sum(x))


plot_heatmap(prune_taxa(seqs,ps.rel),"NMDS","bray",sample.label = "Description","Genus", first.sample = "Plate3-A1")


# save image as svg





