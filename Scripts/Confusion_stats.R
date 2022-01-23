library(magrittr)
library(data.table)
library(BarCluster)
library(ggplot2)

set.seed(42)

source("Constants.R")

# cluster assignments
fl <- list.files("../Processed_data", pattern = "assignments", full.names = TRUE)

# barcodes
bcs <- list.files("../Data_barcodes", full.names = TRUE)

# loop to get confusion metrics from all data then subset on top 15
ct <- lapply(sn_v, function(sn){

  dt <- fl[fl %like% sn] %>% data.table::fread()

  bt <- bcs[bcs %like% sn] %>% data.table::fread()

  bt %>% setnames("cellID", "rn")

  wl <- bt[, .N, by = "Barcode"] %>% .[order(-N)] %>%
    .[, Barcode %>% unique %>% .[1:15]]

  ctt <- lapply(dt[, alpha %>% unique], function(a){

    d <- dt[alpha == a, .SD, .SDcols = c("rn", "Group")]

    ct <- cast_confusion(d, bt)

    ct[, alpha := a]

    return(ct)

  }) %>% data.table::rbindlist()

  ctt[, sample_id := sn]

  return(ctt[barcode %chin% wl])

}) %>% data.table::rbindlist()

# plot
pll <- lapply(sn_v, function(sn){

  # whitelisted cohens_k
  bt <- bcs[bcs %like% sn] %>% data.table::fread()

  wl <- bt[, .N, by = "Barcode"] %>% .[order(-N)] %>%
    .[, Barcode %>% unique %>% .[1:15]]

  c <- ct[sample_id == sn & barcode %chin% wl]

  p2 <- ggplot(c, aes(x = as.factor(alpha), y = cohens_k)) +
          geom_boxplot(fill = "dodgerblue", size = 0.5) +
          theme_bw() +
          theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(color = "black", size = 8),
          axis.title = element_text(face = "bold", color = "black", size = 8),
          panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_line(size = 0.5),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
          xlab("\u03B1 value") +
          ylab("Cohen's \u03BA\nbarcodes to clusters\n") +
          ggtitle(tv[names(tv) == sn])

  fn <- paste("../Figs/Confusion/", sn, "_cohens_k_supplemental.pdf", sep = "")

  ggsave(plot = p2, fn, width = 3.5, height = 2.5, limitsize = FALSE)

  return(p2)

})
