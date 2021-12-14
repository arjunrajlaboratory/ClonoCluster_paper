library(magrittr)
library(data.table)
library(BarCluster)
library(ggplot2)
library(WebGestaltR)

sn_v <- c(paste("YG", 1:3, sep = ""), "Kleind2ws", "CJ")

tv <- c("High dose BRAFi - 1", "High dose BRAFi - 2", "Low dose BRAFi", "In Vitro Hematopoiesis Day 2", "Cardio-directed iPSC")

names(tv) <- sn_v

# cutoffs determined by max alpha before cluster number increases
mt <- data.table::data.table(sample_id = sn_v,
  alevel = c(0.55, 0.6, 0.6, 0.75, 0.55))

mt[, alevel2 := (alevel / 2) %>% round(digits = 1)]

fl <- list.files("../Processed_data/", full.names = TRUE, pattern = "rearrangement_ORA")

hml <- lapply(1:nrow(mt), function(r){

  sn <- mt[r, sample_id]

  a <- mt[r, alevel]

  a2 <- mt[r, alevel2]

  dt <- fl[fl %like% sn] %>% data.table::fread()

  dt <- dt[alpha %in% c(a,a2,1)]

  dt[FDR > 0.05, enrichmentRatio := 1]

  dt[, lER := log(enrichmentRatio, base = 2)]

  dt[, max := lER == max(lER), by = c("lin", "alpha")]

  wl <- dt[max == TRUE & FDR < 0.05 & alpha %in% c(a,a2,1), description %>% unique]

  return(data.table::data.table(sn, wl))

}) %>% data.table::rbindlist()

wl <- hml[, sn %>% unique %>% length, by = wl] %>% .[V1 > 2, wl] %>% sort

hml <- lapply(1:nrow(mt), function(r){

  sn <- mt[r, sample_id]

  a <- mt[r, alevel]

  a2 <- mt[r, alevel2]

  dt <- fl[fl %like% sn] %>% data.table::fread()

  dt <- dt[alpha %in% c(a, a2, 1)]

  dt[FDR > 0.05, enrichmentRatio := 1]

  dt[, lER := log(enrichmentRatio, base = 2)]

  dt[, max := lER == max(lER), by = c("lin", "alpha")]

  dt %<>% .[order(-lER)]

  m <- dt[description %chin% wl, .SD, .SDcols = c("description", "alpha", "lER")] %>%
    unique(by = c("description", "alpha"))

  m <- m %>% dcast(description ~ alpha, value.var = "lER") %>% as.data.table

  rn <- m[, description]

  m <- m %>% as.data.table %>% .[, 2:ncol(m)] %>% as.matrix

  rownames(m) <- rn

  while(ncol(m) < 4){

    m <- cbind(m, replicate(n = nrow(m), as.numeric(NA)))

  }

  while(nrow(m) < 7){

    m <- rbind(m, replicate(n = ncol(m), as.numeric(NA)))

  }

  colnames(m)[colnames(m) == 0] <- "Transcriptome"
  colnames(m)[colnames(m) == 1] <- "Transcriptome to Barcodes"
  colnames(m)[colnames(m) == a2] <- "Transcriptome to Low alpha"
  colnames(m)[colnames(m) == a] <- "Transcriptome to High alpha"

  cv <- c("Transcriptome", "Transcriptome to Barcodes",
  "Transcriptome to Low alpha", "Transcriptome to High alpha")

  cn <- cv[!cv %chin% colnames(m)]

  colnames(m)[colnames(m) == ""] <- cn

  m <- m[, c("Transcriptome to Low alpha", "Transcriptome to High alpha", "Transcriptome to Barcodes")]

  rv <- wl[!wl %chin% rownames(m)]

  rownames(m)[rownames(m) == ""] <- rv

  m <- m[wl, ]

  l <- TRUE

  #if (r == 3) l <- TRUE

  hm <- pheatmap::pheatmap(m, cluster_cols = FALSE, cluster_rows = FALSE,
    breaks = seq(0,10, by = 0.1),
    fontsize_row = 8,
    fontsize_col = 8,
    cellheight = 10,
    cellwidth = 15,
    color = colorRampPalette(RColorBrewer::brewer.pal(n = 7,
    name = "Reds"))(100),
    na_col = "grey", main = tv[names(tv) == sn],
    legend = l,
    fontsize = 8)

  # rotation of x axis labels
  hm$gtable$grobs[[3]]$rot <- 90
  hm$gtable$grobs[[3]]$hjust <- 1
  hm$gtable$grobs[[3]]$vjust <- 0.5

  return(hm$gtable)

})

hml <- hml[c(3,1:2)]

for (i in 1:(length(hml) - 1)){

  hml[[i]]$grobs[[3]]$label <- c("", "", "")
}

plot <- cowplot::plot_grid(plotlist = hml, ncol = 1, scale = 1)

ggsave(plot = plot, filename = "../Figs/Rearrangement_genesets_hm.png", width = 3.5, height = 9, limitsize = FALSE)
