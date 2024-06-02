library(fgsea)
library(data.table)
library(ggplot2)

data(examplePathways)
data(exampleRanks)
set.seed(42)

gmt_pathways <- getGmt("~/Documents/GitHub/591_files/R_Shiny_app/rnaseq_rewritten/gmt_files/h.all.v2023.2.Hs.symbols.gmt.txt")
print(class(gmt_pathways))  # Verify it's GeneSetCollection
print(head(gmt_pathways)) # Verify contents

print(class(examplePathways))
print(head(examplePathways))

print(head(exampleRanks))

fgseaRes <- fgsea(pathways = gmtPathways, 
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 500)

head(fgseaRes[order(pval), ])