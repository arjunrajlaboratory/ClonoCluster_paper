library(magrittr)
library(data.table)
library(ClonoCluster)
library(ggplot2)

set.seed(42)

source("Paper/extractionScripts/Constants.R")

# get hybrid cluster marker files
fl <- list.files(pdir, pattern = "auc_linclust", full.names = TRUE)

lapply(fl, function(f){

  dt <- data.table::fread(f)

  sn <- basename(f) %>% stringr::str_extract("^[^_]+")

  dt[, sample_id := sn]

  # not doing barcodes on these plots
  dt <- dt[alpha < 1]

  # take the top marker per cluster
  topmarkers <- dt %>% split(by = c("alpha", "lin")) %>%

    lapply(., function(t){

      return(t[auc == max(auc)])

  }) %>% data.table::rbindlist()

  # subset on top markers
  dt <- dt[rn %chin% topmarkers[, rn]]

  # get max value for any cluster for that marker
  dt <- dt[, max(auc), by = c("rn", "alpha")]

  # convert to decimal
  dt[, auc := V1 / 100]

  # high and low alpha levels
  al <- mt[sample_id == sn, .(alevel2, alevel)] %>% unlist

  # label
  ll <- c("Transcriptome", paste("Low @\n@ = ", al[1], sep = ""),
  paste("High @\n@ = ", al[2], sep = ""))

  dt[alpha == 0, lab := ll[1]]

  dt[alpha == al[1], lab := ll[2]]

  dt[alpha == al[2], lab := ll[3]]

  # make wide table
  dtw <- dt %>% dcast(rn ~ lab, value.var = "auc")

  # order columns
  dtw <- dtw[, .SD, .SDcols = c("rn", ll)]

  # convert to matrix
  m <- ClonoCluster::dt2m(dtw)

  # heatmap
  hm <- pheatmap::pheatmap(m,
    cluster_cols = FALSE,
    breaks = seq(0.5,1, length.out = 101),
    fontsize_row = 6,
    fontsize_col = 6,
    cellheight = 10,
    cellwidth = 20,
    # color = colorRampPalette(RColorBrewer::brewer.pal(n = 7,
    # name = "RdYlBu"))(100),
    na_col = "grey",
    main = tv[names(tv) == sn],
    legend = TRUE,
    fontsize = 8)

  # rotation of x axis labels
  hm$gtable$grobs[[4]]$rot <- 45
  hm$gtable$grobs[[4]]$hjust <- 1
  hm$gtable$grobs[[4]]$vjust <- 1


  ggsave(plot = hm$gtable,
    filename = paste(fdir, "Turnover_hm/", sn, "_top_hm.pdf", sep = ""),
    width = 3.5, height = 8, limitsize = FALSE)

})
