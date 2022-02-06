library(magrittr)
library(data.table)
library(BarCluster)
library(ggplot2)
library(WebGestaltR)

# NJ Neuro-directed iPSC alevel = 0.8

# vector of short names
sn_v <- c(paste("YG", 1:3, sep = ""), "Kleind2ws", "CJ",
paste("MDA", 1:2, sep = ""), paste("WM983", 1:2, sep = ""))

# vector of long names
tv <- c("High dose BRAF inhibitor\nReplicate 1", "High dose BRAF inhibitor\nReplicate 2",
"Low dose BRAF inhibitor", "Weinreb et al. (2020)\nIn vitro hematopoiesis",
"Cardio-directed iPSC", "Paclitaxel-treated breast cancer cells\nReplicate 1", "Paclitaxel-treated breast cancer cells\nReplicate 2",
"WM983 BRAF inhibitor\nReplicate 1", "WM983 BRAF inhibitor\nReplicate 2")

names(tv) <- sn_v

# cutoffs determined by max alpha before cluster number increases
mt <- data.table::data.table(sample_id = sn_v,
  alevel = c(0.55, 0.6, 0.6, 0.75, 0.55, 0.4, 0.7, 0.7, 0.65))

mt[, alevel2 := (alevel / 2) %>% round(digits = 1)]

set.seed(42)

# figs dir path
fdir <- "Paper/plots/"

pdir <- "Paper/extractedData/"

dbdir <- "Data/Data_barcodes/"

dgdir <- "Data/Data_genes/"
