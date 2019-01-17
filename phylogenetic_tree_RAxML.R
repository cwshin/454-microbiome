# Phylogenetic tree reconstruction

# Preparation
install.packages("seqinr")
library("DECIPHER")
library("phangorn")

# made SV tables
seqtab.nochim1 <- seqtab.nochim
colnames(seqtab.nochim1) <- paste0("SV", 1:ncol(seqtab.nochim1))

# and change the names again (the order was the same as in seqtab.nochim)
taxtab1 <- taxa
rownames(taxtab1) <- paste0("SV", 1:nrow(taxtab1))

# The SV table and SV sequences were exported also:
  export_taxa_table_and_seqs = function(seqtab.nochim, file_seqtab, file_seqs) {
    seqtab.t = as.data.frame(t(seqtab.nochim))
    seqs = row.names(seqtab.t)
    row.names(seqtab.t) = paste0("SV", 1:nrow(seqtab.t))
    outlist = list(data_loaded = seqtab.t)
    mctoolsr::export_taxa_table(outlist, file_seqtab)
    seqs = as.list(seqs)
    seqinr::write.fasta(seqs, row.names(seqtab.t), file_seqs)
  }
export_taxa_table_and_seqs(seqtab.nochim,"/path_to your_file/SV_table.txt",
                           "/path_to your_file/SV_seqs.fa")

## submitted to https://www.phylo.org/portal2/home.action
## RAxML Blackbox used

# Then you incorporate everything to your phyloseq object:
all(rownames(seqtab.nochim) %in% samdf$Name)
rownames(samdf) <- samdf$Name

tree <- read_tree("D:/R/korea_microbiome/raxml/korea_microbiota1/RAxML_bipartitionsBranchLabels.result")
ps <- phyloseq(tax_table(taxtab1), sample_data(samdf),otu_table(seqtab.nochim1, taxa_are_rows = FALSE), phy_tree(tree))
saveRDS(ps, file="ps.rds")

### Move to phyloseq.R
