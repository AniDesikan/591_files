# Load necessary libraries
library(tibble)
library(dplyr)
library(fgsea)
library(GSEABase)

# Load GMT file
gmt_pathways <- getGmt("~/Documents/GitHub/591_files/R_Shiny_app/rnaseq_rewritten/gmt_files/h.all.v2023.2.Hs.symbols.gmt.txt")
print(class(gmt_pathways))  # Verify it's GeneSetCollection
print(gmt_pathways)  # Verify contents

# Convert GeneSetCollection to a list of pathways
gmt_list <- lapply(gmt_pathways, geneIds)
names(gmt_list) <- names(gmt_pathways)
print(head(gmt_list))  # Print the first few pathways to check the structure

# Simulate EA data
ea_data <- tibble(
  Gene = c("ACAA2", "ACADL", "ACADM", "ACADS", "ACADS", "ACO2"),
  baseMean = c(1924.73336, 1401.85324, 1437.38440, 1639.32133, 1890.94223, 92.78162),
  log2FoldChange = c(0.7456007, -0.5697995, -0.3232443, -0.4887181, 0.7568135, 1.4765208),
  lfcSE = c(0.2478956, 0.1809261, 0.1588603, 0.1957328, 0.2457001, 0.6108521),
  stat = c(3.007721, -3.149348, -2.034771, -2.496864, 3.080232, 2.417149),
  pvalue = c(0.002632145, 0.001636350, 0.041873944, 0.012529707, 0.002068391, 0.015642598),
  padj = c(0.01834634, 0.01376234, 0.11666562, 0.05184247, 0.01582544, 0.05968297)
)

make_ranked_log2fc <- function(dataf) {
  # Remove duplicates by averaging log2FoldChange for duplicate genes
  dataf <- dataf %>%
    group_by(Gene) %>%
    summarize(log2FoldChange = mean(log2FoldChange, na.rm = TRUE)) %>%
    ungroup()
  
  new_table <- dataf %>%
    mutate(rank = rank(log2FoldChange)) %>%
    arrange(rank) %>%
    dplyr::select(Gene, log2FoldChange) %>%
    filter(is.finite(log2FoldChange))
  
  new_vector <- pull(new_table, log2FoldChange, Gene)
  return(new_vector)
}

# Create ranked log2fc list
rnk_list <- make_ranked_log2fc(ea_data)
print(rnk_list)  # Verify it's a named numeric vector

# Check common identifiers between the ranked list and the gene sets
common_genes <- intersect(names(rnk_list), unique(unlist(gmt_list)))
print(common_genes)

# Check if there are common genes
if (length(common_genes) == 0) {
  stop("No common genes found between the ranked list and the gene sets.")
}

# Subset the ranked list to include only common genes
subset_rnk_list <- rnk_list[names(rnk_list) %in% common_genes]
print(subset_rnk_list)

# Run fgsea with the subset ranked list
fgsea_results <- fgsea(pathways = gmt_list, stats = subset_rnk_list, minSize = 15, maxSize = 500)
print(fgsea_results)

# Convert to tibble
fgsea_tibble <- as_tibble(fgsea_results)
print(fgsea_tibble)
