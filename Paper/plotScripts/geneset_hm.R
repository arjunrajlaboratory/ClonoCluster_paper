library(magrittr)
library(data.table)
library(ClonoCluster)
library(ggplot2)
library(WebGestaltR)

set.seed(42)

source("Paper/extractionScripts/Constants.R")

# read overrepresentation analysis results
fl <- list.files(pdir, full.names = TRUE, pattern = "rearrangement_ORA")

hml <- lapply(1:nrow(mt), function(r){

  sn <- mt[r, sample_id]

  a <- mt[r, alevel]

  a2 <- mt[r, alevel2]

  dt <- fl[fl %like% sn] %>% data.table::fread()

  # subset on alphas of interest
  dt <- dt[alpha %in% c(a,a2,1)]

  # if FDR is nonsignif, call enrichment 1 (not enriched)
  dt[FDR > 0.05, enrichmentRatio := 1]

  # log transform for better visualization
  dt[, lER := log(enrichmentRatio, base = 2)]

  # take max value of log enrichment ratio
  # (computed by hybrid to contributing transcriptome pair)
  # for any hybrid cluster at each alpha level
  dt[, max := lER == max(lER), by = c("lin", "alpha")]

  # white list of significantly enriched genesets at any level in this dataset
  wl <- dt[max == TRUE & FDR < 0.05 & alpha %in% c(a,a2,1), description %>% unique]

  return(data.table::data.table(sn, wl))

}) %>% data.table::rbindlist()

# gene sets with at least 3 samples with any enrichment
wl <- hml[, sn %>% unique %>% length, by = wl] %>% .[V1 > 2, wl] %>% sort

wl1 <- wl

# negative control genesets
wl2 <- c("transmembrane receptor protein kinase activity",
        "bHLH transcription factor binding",
        "activating transcription factor binding",
        "transcription factor complex",
        "cyclin binding",
        "cyclin-dependent protein kinase activity",
        "protein tyrosine kinase binding")

wl <- c(wl1, wl2)

# low alpha heatmap
# loop over samples and cbind matrices
hm <- lapply(1:nrow(mt), function(r){

  print(r)

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

  # fill in insignificant comparisons with NA
  while(ncol(m) < 4){

    m <- cbind(m, replicate(n = nrow(m), as.numeric(NA)))

  }

  while(nrow(m) < length(wl)){

    m <- rbind(m, replicate(n = ncol(m), as.numeric(NA)))

  }

  cv <- c("Transcriptome", "Barcode\nreorganization\nmarkers",
  "Low \u03B1\nreorganization\nmarkers", "High \u03B1\nreorganization\nmarkers")

  colnames(m)[colnames(m) == 0] <- cv[1]
  colnames(m)[colnames(m) == 1] <- cv[2]
  colnames(m)[colnames(m) == a2] <- cv[3]
  colnames(m)[colnames(m) == a] <- cv[4]

  cn <- cv[!cv %chin% colnames(m)]

  colnames(m)[colnames(m) == ""] <- cn

  m <- m[, c(cv[3])]

  m[m == 0] <- as.numeric(NA)

  m %<>% as.matrix

  colnames(m) <- sn

  rv <- wl[!wl %chin% rownames(m)]

  rownames(m)[rownames(m) == ""] <- rv

  m <- m[wl, ]

  return(m)

}) %>% purrr::reduce(cbind)

colnames(hm) <- mt[, sample_id]

hm <- hm[, !colnames(hm) %like% "NJ"]

colnames(hm) <- tv[1:ncol(hm)]

# heatmap of gene sets
hm1 <- pheatmap::pheatmap(hm[wl1, ], cluster_cols = FALSE, cluster_rows = FALSE,
    breaks = seq(0,10, by = 0.1),
    fontsize_row = 6,
    fontsize_col = 6,
    cellheight = 10,
    cellwidth = 35,
    color = colorRampPalette(RColorBrewer::brewer.pal(n = 7,
    name = "Reds"))(100),
    na_col = "grey",
    main = "Low \u03B1 Reorganization Gene Sets",
    legend = TRUE,
    fontsize = 8)

# heatmap of control genesets
hm2 <- pheatmap::pheatmap(hm[wl2, ], cluster_cols = FALSE, cluster_rows = FALSE,
    breaks = seq(0,10, by = 0.1),
    fontsize_row = 6,
    fontsize_col = 6,
    cellheight = 10,
    cellwidth = 35,
    color = colorRampPalette(RColorBrewer::brewer.pal(n = 7,
    name = "Reds"))(100),
    na_col = "grey",
    legend = TRUE,
    fontsize = 8)

# rotation of x axis labels
hm1$gtable$grobs[[3]]$rot <- 45
hm1$gtable$grobs[[3]]$hjust <- 1
hm1$gtable$grobs[[3]]$vjust <- 1

hm2$gtable$grobs[[2]]$rot <- 45
hm2$gtable$grobs[[2]]$hjust <- 1
hm2$gtable$grobs[[2]]$vjust <- 1

hmp <- cowplot::plot_grid(hm1$gtable, hm2$gtable, ncol = 1)

ggsave(plot = hmp,
    filename = paste(fdir, "reorg_hm/low_rearrangement_genesets_hm.pdf", sep = ""),
    width = 7.5, height = 8, limitsize = FALSE)

# same as above loop but for high alpha heatmap
hm <- lapply(1:nrow(mt), function(r){

  print(r)

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

  while(nrow(m) < length(wl)){

    m <- rbind(m, replicate(n = ncol(m), as.numeric(NA)))

  }

  cv <- c("Transcriptome", "Barcode\nreorganization\nmarkers",
  "Low \u03B1\nreorganization\nmarkers", "High \u03B1\nreorganization\nmarkers")

  colnames(m)[colnames(m) == 0] <- cv[1]
  colnames(m)[colnames(m) == 1] <- cv[2]
  colnames(m)[colnames(m) == a2] <- cv[3]
  colnames(m)[colnames(m) == a] <- cv[4]

  cn <- cv[!cv %chin% colnames(m)]

  colnames(m)[colnames(m) == ""] <- cn

  m <- m[, c(cv[4])]

  m[m == 0] <- as.numeric(NA)

  m %<>% as.matrix

  colnames(m) <- sn

  rv <- wl[!wl %chin% rownames(m)]

  rownames(m)[rownames(m) == ""] <- rv

  m <- m[wl, ]

  return(m)

}) %>% purrr::reduce(cbind)

colnames(hm) <- mt[, sample_id]

hm <- hm[, !colnames(hm) %like% "NJ"]

colnames(hm) <- tv[1:ncol(hm)]

hm <- pheatmap::pheatmap(hm, cluster_cols = FALSE, cluster_rows = FALSE,
    breaks = seq(0,10, by = 0.1),
    fontsize_row = 6,
    fontsize_col = 6,
    cellheight = 10,
    cellwidth = 35,
    color = colorRampPalette(RColorBrewer::brewer.pal(n = 7,
    name = "Reds"))(100),
    na_col = "grey",
    main = "High \u03B1 Reorganization Gene Sets",
    legend = TRUE,
    fontsize = 8)

# rotation of x axis labels
hm$gtable$grobs[[3]]$rot <- 45
hm$gtable$grobs[[3]]$hjust <- 1
hm$gtable$grobs[[3]]$vjust <- 1

ggsave(plot = hm$gtable,
    filename = paste(fdir, "reorg_hm/high_rearrangement_genesets_hm.pdf", sep = ""),
    width = 7.5, height = 4, limitsize = FALSE)
