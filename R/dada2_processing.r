############################################################
############################################################
##  BACTERIAL READS DADA2 PROCESSING  ######################
############################################################
############################################################

library(dada2); packageVersion("dada2")
path <- "~/path_to_sequence_directory" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
# Only need forward reads for this analysis
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
#fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
#filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
#names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, trimLeft = 23, truncLen=c(240),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
errF <- learnErrors(filtFs, multithread=TRUE)
#errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
#dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
#mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
#head(mergers[[1]])
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim.16s <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim.16s)
sum(seqtab.nochim.16s)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim.16s))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
taxa.16s <- assignTaxonomy(seqtab.nochim.16s, "~/silva_nr_v132_train_set.fa.gz", tryRC = TRUE, multithread=TRUE)
taxa.16s <- addSpecies(taxa.16s, "~/silva_species_assignment_v132.fa.gz", tryRC = TRUE)

taxa.16s.print <- taxa.16s # Removing sequence rownames for display only
rownames(taxa.16s.print) <- NULL
head(taxa.16s.print)
save(seqtab.nochim.16s,taxa.16s, file = "Bacterial_Dada2_seq_tax_output.RData")

############################################################
############################################################
##  FUNGAL READS DADA2 PROCESSING  #########################
############################################################
############################################################

library(dada2); packageVersion("dada2")
path <- "~/path_to_sequence_directory" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
#fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
#filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
#names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, trimLeft = 23, truncLen=c(240),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
errF <- learnErrors(filtFs, multithread=TRUE)
#errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
#dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
#mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
#head(mergers[[1]])
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim.its <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim.its)
sum(seqtab.nochim.its)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim.its))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
taxa.its <- assignTaxonomy(seqtab.nochim.its, "~/sh_general_release_dynamic_s_04.02.2020.fasta", tryRC = TRUE, multithread=TRUE)
colnames(taxa.its) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa.its= gsub("^.*?__","",taxa.its)
taxa.its.print <- taxa.its # Removing sequence rownames for display only
rownames(taxa.its.print) <- NULL
head(taxa.its.print)
save(seqtab.nochim.its,taxa.its, file = "Fungal_Dada2_seq_tax_output.RData")
