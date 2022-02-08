library(magrittr)
library(data.table)
library(ClonoCluster)
library(ggplot2)

set.seed(42)

source("Paper/extractionScripts/Constants.R")

lapply(sn_v, function(sn){

  set.seed(42)

  res <- 1

  if (sn %like% "Klein") res <- 0.6

  # prep gene counts
  gt <- data.table::fread(paste(dgdir, sn, "_genes.txt", sep = ""))

  rn <- gt[, rn]

  gt <- gt[, .SD, .SDcols = 2:ncol(gt)] %>% as.matrix %>% t

  colnames(gt) <- rn

  if(sn == "CJ") rownames(gt) %<>% stringr::str_replace("-[0-9]$", "")

  if(sn == "NJ") rownames(gt) <- rownames(gt) %>%
      stringr::str_replace("^S[0-9]_", "") %>%
      stringr::str_replace("-[0-9]$", "")

  # prep barcodes
  bt <- paste(dbdir, sn, "_barcodes.tsv", sep = "") %>% data.table::fread()

  bt <- bt[cellID %chin% rownames(gt)]

  gt <- gt[bt[, cellID], ]

  bt %>% data.table::setnames("cellID", "rn")

  # pca
  irl <- irlba_wrap(gt, npc = 100)

  # warp pca at these levels
  gl <- lapply(c(0,2.5,5, 7.5, 10), function(w){

    mo <- barcode_warp(irl, bt, w)

    return(mo)

  })

  names(gl) <- c(0,0.25,0.5, 0.75, 1)

  # generate warped umaps
  gl <- lapply(gl %>% seq_along, function(a){

      g <- umap_matrix(gl[[a]])

      g[, UMAP_1 := UMAP_1 %>% scale]

      g[, UMAP_2 := UMAP_2 %>% scale]

      g[, s := names(gl)[a] %>% as.numeric]

  }) %>% data.table::rbindlist()

  gl <- merge(gl, bt, by = "rn")

  gl[, bcs := rn %>% unique %>% length, by = "Barcode"]

  # mark singlets
  gl[bcs == 1, Barcode := "Singlet"]

  bct <- gl[order(-bcs), Barcode %>% unique %>% .[1:3]]

  bct %<>% c("Singlet")

  cv <- c("grey", c25[4])

  names(cv) <- cv

  title_vector <- c(paste("Barcode #", 1:3), "Singlets")

  # Plot 3 barcodes and singlets for each sample
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
      ggtitle(title_vector[n])

    return(p)

  })

  plot <- cowplot::plot_grid(plotlist = pl, ncol = 1, scale = 0.9)

  fn <- paste(fdir, "warped/", sn, "_warped.pdf", sep = "")

  ggsave(plot = plot, fn, width = 6, height = 4 * 6/5)

})
