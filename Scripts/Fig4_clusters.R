library(magrittr)
library(data.table)
library(BarCluster)
library(ggplot2)

set.seed(42)

sn_v <- c(paste("YG", 1:3, sep = ""), paste("Kleind", c(2,4,6), sep = ""), "CJ", "NJ")

sn <- "YG1"

res <- 1

if (sn %like% "Klein") res <- 0.6

gt <- data.table::fread(paste("../Data_genes/", sn, "_genes.txt", sep = ""))

rn <- gt[, rn]

gt <- gt[, .SD, .SDcols = 2:ncol(gt)] %>% as.matrix %>% t

colnames(gt) <- rn

if(sn == "CJ") rownames(gt) %<>% stringr::str_replace("-[0-9]$", "")

if(sn == "NJ") rownames(gt) <- rownames(gt) %>%
    stringr::str_replace("^S[0-9]_", "") %>%
    stringr::str_replace("-[0-9]$", "")

bt <- paste("../Data_barcodes/", sn, "_barcodes.tsv", sep = "") %>% data.table::fread()

bt <- bt[cellID %chin% rownames(gt)]

gt <- gt[bt[, cellID], ]

bt %>% data.table::setnames("cellID", "rn")

irl <- irlba_wrap(gt, npc = 100)

dl <- BarCluster::barcluster(irl, bt, alpha = c(0,0.25,0.5, 0.75, 1), beta = 0.1, res = res)

gl <- lapply(c(0,0.25,0.5, 0.75, 1), function(w){

  mo <- barcode_warp(irl, bt, w)

  return(mo)

})

names(gl) <- c(0,0.25,0.5, 0.75, 1)

gl <- lapply(gl %>% seq_along, function(a){

    g <- umap_matrix(gl[[a]])

    g[, UMAP_1 := UMAP_1 %>% scale]

    g[, UMAP_2 := UMAP_2 %>% scale]

    g[, s := names(gl)[a] %>% as.numeric]

}) %>% data.table::rbindlist()

gl <- merge(gl, bt, by = "rn")

gl[, bcs := rn %>% unique %>% length, by = "Barcode"]

bct <- gl[order(-bcs), Barcode %>% unique %>% .[3:5]]

cv <- c("grey", c25[4])

names(cv) <- cv

pl <- lapply(bct %>% seq_along, function(n){

  bc <- bct[n]

  dt <- gl %>% data.table::copy()

  dt[Barcode != bc, paint := cv[1]]

  dt[Barcode == bc, paint := cv[2]]

  dt %>% setkey(paint)

  dt %<>% .[order(paint)]

  p <- ggplot(dt, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(col = paint), size = 0.01) +
    geom_point(data = dt[Barcode == bc], aes(col = paint), size = 0.01) +
    facet_wrap(~s, nrow = 1) +
    scale_color_manual(values = cv) +
    theme_void() +
    theme(legend.position = "none", strip.text = element_blank(),
    plot.title = element_text(size = 8, face = "bold", family = "sans", hjust = 0.5)) +
    ggtitle(paste("Barcode #", n))

  return(p)

})

plot <- cowplot::plot_grid(plotlist = pl, ncol = 1, scale = 0.9)

ggsave(plot = plot, "../Figs/Fig4_clusters.png", width = 6, height = 3 * 6/5)
