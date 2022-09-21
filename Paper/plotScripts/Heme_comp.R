library(magrittr)
library(data.table)
library(ClonoCluster)
library(ggplot2)
library(entropy)

# find optimal alpha for predicing various cell types
# 2 column sankey with cells in purple

set.seed(42)

source("Paper/extractionScripts/Constants.R")

meta <- "Paper/extractedData/GSM4185642_stateFate_inVitro_metadata.txt" %>% fread

cl <- "Paper/extractedData/Kleind2ws_cluster_assignments.txt" %>% fread

m6 <- meta[`Time point` == 2, .SD %>% unique,
.SDcols = c("Cell barcode", "Cell type annotation")]

names(m6) <- c("rn", "label")

# blacklist barcodes with indeterminant day6 phenotype from multiple libraries
blacklist <- m6[, .N, by = "rn"] %>% .[N > 1, rn]

m6 <- m6[!rn %chin% blacklist]

cl <- merge(cl, m6, by = "rn")

# entropy
et <- lapply(cl[, alpha %>% unique], function(a){

  d <- cl[alpha == a]

  et <- lapply(1:1000, function(seed){

    set.seed(seed)

    d <- d[sample(1:nrow(d), size = floor(nrow(d) / 3))]

    et <- table(d$Group, d$label) %>% data.table::as.data.table(keep.rownames = TRUE)

    et[, ent := entropy(N, units = "log2", method = "Laplace"), by = "V2"]

    et[, alpha := a]

    et[, seed := seed]

    return(et[, .SD %>% unique, .SDcols = c("V2", "alpha", "ent", "seed")])

  }) %>% data.table::rbindlist()

  return(et)

}) %>% data.table::rbindlist()

et[, m := mean(ent), by = c("V2", "alpha")]

et[, lower := ent %>% quantile(0.025), by = c("V2", "alpha")]

et[, upper := ent %>% quantile(0.975), by = c("V2", "alpha")]

pt <- et[alpha %in% c(0, 0.4, 0.75), .(m, upper, lower, V2, alpha)] %>% unique

pt[V2 == "Baso", V2 := "Basophil"]

pt[V2 == "Mast", V2 := "Mast cell"]

pt[V2 == "Meg", V2 := "Megakaryocyte"]

pt[V2 == "Ccr7_DC", V2 := "Migratory DC"]

# make this pretty
p <- ggplot(pt, aes(x = as.factor(alpha), y = m)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), color = "black", width = 0.2) +
  geom_point(color = "dodgerblue", size = 2) +
  facet_wrap(~V2, scales = "free_y") +
  theme_bw() +
  ClonoCluster::ttheme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) +
  scale_x_discrete(labels =
    # use @ sign for easy alpha swap in inkscape because pdf hates symbols
    c("Transcriptome\n@ = 0",
    "Low @\n@ = 0.4",
    "High @\n@ = 0.75")) +
  ylim(0.8454706, 2.5) +
  xlab(NULL) +
  ylab("Entropy (95% confidence level)")

ggsave(plot = p, filename = file.path(fdir, "entropy", "entropy.pdf"), width = 7.5, height = 7)
