# Install dada2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.8")

## DADA2 pipeline tutorial: https://benjjneb.github.io/dada2/tutorial.html
## Microbiome workflow http://web.stanford.edu/class/bios221/MicrobiomeWorkflowII.html

# loading dada2
library(dada2); packageVersion("dada2")

# path arrange
path <- "D:/R/korea_microbiome/FASTQ"
list.files(path)

# Forward fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect read quality profiles
plotQualityProfile(fnFs[1:52])
## Result: quality median score decreased after 250 -> applied to Filter and trim

# Filter and trim
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, truncLen=c(250),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

## FAQ (http://benjjneb.github.io/dada2/faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data)
## We recommend additional (relative to the tutorial) filtering of 454 sequences by maximum length:
## filterAndTrim(..., maxLen=XXX) # XXX depends on the chemistry - was not applied

# Learn the error rates
errF <- learnErrors(filtFs, multithread=FALSE)
plotErrors(errF, nominalQ=TRUE)

# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
## Name the derep-class objects by the sample names
names(derepFs) <- sample.names

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=FALSE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
dadaFs[[1]]

# Merge paired reads - skipped in R454 data

# Construct sequence table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

## Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=F, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "D:/R/korea_microbiome/tax/rdp_train_set_16.fa.gz", multithread=F, tryRC=T)
## used taxonomy database: RDP v16, see other database: https://benjjneb.github.io/dada2/training.html
## may use another methods in the tuturial, see 'assign taxonomy' section

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

### Reference: https://f1000research.com/articles/5-1492/v2
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# Construct phylogenetic tree - DECIPHER package needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER", version = "3.8")

library(DECIPHER); packageVersion("DECIPHER")

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

# phangorn tree asignment was NOT working 
# RAxML instead: https://github.com/amkozlov/raxml-ng/

# Combine data into a phyloseq object
mimarks_path <- "D:/R/korea_microbiome/clinical_dataset.csv"
samdf <- read.csv(mimarks_path, header=TRUE)
all(rownames(seqtab.nochim) %in% samdf$Name) # TRUE
rownames(samdf) <- samdf$Name

## MOVE to phyloseq.R
