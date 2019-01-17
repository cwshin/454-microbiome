### Phylosequence after DADA2 ###
# Front page: http://joey711.github.io/phyloseq/index.html
# Workflow reference: https://f1000research.com/articles/5-1492/v2

# Install phyloseq, ggplot

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5", version = "3.8") # when phyloseq needs 'rhdf5' package
BiocManager::install("phyloseq")

# loading phyloseq, ggplot

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

ps # ps from the phylogenetic tree reconstruction (RAxML)

# Filtering
# Taxonomic Filtering
rank_names(ps)
table(tax_table(ps)[, "Phylum"], exclude = NULL) # Create table, number of features for each phyla
## <NA> 80 - probable artifacts
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)}) # Compute prevalence of each feature, store as data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0)) # Add taxonomy and total read counts to this data.frame
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

filterPhyla = c("Acidobacteria", "Candidatus_Saccharibacteria","Cyanobacteria/Chloroplast",
                "Fusobacteria","Synergistetes","Tenericutes","Verrucomicrobia") # Define phyla to filter
ps1 = subset_taxa(ps0, !Phylum %in% filterPhyla) # Filter entries with unidentified Phylum.
ps1

### Prevalence filtering - not performed here
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps0)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps0)

# Agglomerate taxa - skipped prevalence filtering
# How many genera would be present after filtering?
length(get_taxa_unique(ps1, taxonomic.rank = "Genus"))
ps3 = tax_glom(ps1, "Genus", NArm = TRUE)
h1 = 0.4
ps4 = tip_glom(ps1, h = h1)

# Plot tree
multiPlotTitleTextSize = 8
p1tree = plot_tree(ps1, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(ps4, method = "treeonly",
                   ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))

# group plots together - skip, needs gridExtra function
grid.arrange(nrow = 1, p1tree, p3tree, p4tree)
