library(magrittr)
library(data.table)
library(BarCluster)
library(ggplot2)
library(randomcoloR)

set.seed(42)

source("Constants.R")

dcols <- c("#629B60", "#98B3DF", "#E7DFCC", "#56E884", "#EDE644", "#C2ECA5",
"#D44DA1", "#E68D9E", "#D4F0E6", "#D6D5E6", "#546790", "#5C99E6",
"#9BDFC0", "#E4F2C2", "#7E37D9", "#64EC41", "#5EEABC", "#6465CF",
"#DBAC48", "#A7B73A", "#E8B8E1", "#E4E97F", "#DFCF89", "#5EE7E1",
"#B691A5", "#DFB59C", "#52BEE4", "#E74B64", "#78A59C", "#A797E0",
"#A6E770", "#96DAE6", "#A6AB79", "#E48449", "#E050E6", "#DF90DB",
"#B47EE6", "#98E699", "#9C7058", "#B5EA3B")

lapply(sn_v, function(sn){

  set.seed(42)

  res <- 1

  if (sn %like% "Klein") res <- 0.6

  # prep gene counts
  gt <- data.table::fread(paste("../Data_genes/", sn, "_genes.txt", sep = ""))

  rn <- gt[, rn]

  gt <- gt[, .SD, .SDcols = 2:ncol(gt)] %>% as.matrix %>% t

  colnames(gt) <- rn

  if(sn == "CJ") rownames(gt) %<>% stringr::str_replace("-[0-9]$", "")

  if(sn == "NJ") rownames(gt) <- rownames(gt) %>%
      stringr::str_replace("^S[0-9]_", "") %>%
      stringr::str_replace("-[0-9]$", "")

  # prep barcodes
  bt <- paste("../Data_barcodes/", sn, "_barcodes.tsv", sep = "") %>% data.table::fread()

  bt <- bt[cellID %chin% rownames(gt)]

  gt <- gt[bt[, cellID], ]

  bt %>% data.table::setnames("cellID", "rn")

  # pca
  irl <- irlba_wrap(gt, npc = 100)

  als <- c(0, mt[sample_id == sn, alevel] %>% unlist)

  if (!sn %chin% c("CJ", "NJ", paste("MDA", 1:2, sep = ""), paste("WM983", 1:2, sep = ""))){

    method <- "index"

  } else {

    method <- "fast"

  }

  dl <- barcluster(irl, bt, alpha = als, res = res, method = method)

  # warp pca at these levels
  gl <- lapply(c(0,5), function(w){

    mo <- engage_warp(irl, bt, w)

    return(mo)

  }) %>% data.table::rbindlist()

  gl <- merge(gl, bt, by = "rn")

  gl[, bcs := rn %>% unique %>% length, by = "Barcode"]

  # mark singlets
  gl[bcs == 1, Barcode := "Singlet"]

  bct <- gl[order(-bcs), Barcode %>% unique %>% .[1:3]]

  marks <- Find_Markers_ROC(dl, gt, n_threads = 6)

  topm <- marks[, .(rn[auc == max(auc)], auc[auc == max(auc)]), by = c("alpha", "Group")]

  topm[, mauc := max(V2), by = c("V1", "alpha")]

  topm <- topm[mauc == V2]

  grt <- topm[alpha == als[2], .SD %>% unique, .SDcols = c("V1", "Group")]

  topm %<>% dcast(V1 ~ alpha, value.var = "V2")

  topm %>% setnames(as.character(als), c("Transcriptome", "High @"))

  topm <- topm[is.na(Transcriptome)]

  topm <- merge(topm, grt, by = "V1")

  topm %<>% .[order(-`High @`)]

  topm <- topm[`High @` > 80]

  merget <- merge(gl, dl, by = "rn", allow.cartesian = TRUE)

  merget %<>% .[order(alpha, warp)]

  merget[, alphaf := paste("@ =", alpha)]

  merget[, warpf := paste("Warp ", warp)]

  merget[, alphaf := factor(alphaf, ordered = TRUE, levels = paste("@ =", als))]

  merget[, warpf := factor(warpf, ordered = TRUE, levels = paste("Warp ", c(0,5)))]

  merget %>% setkey(alphaf, warpf)

  topmt <- marks[, .(rn[auc == max(auc)], auc[auc == max(auc)]), by = c("alpha", "Group")]

  merget <- merge(merget, topmt, by = c("alpha", "Group"))

  subdir <- paste("../Figs/Gene_umaps/", sn, "/", sep = "")

  if(!dir.exists(subdir)) dir.create(subdir)

  pl <- lapply(1:nrow(topm), function(r){

    rnl <- merget[alpha == als[2] & Group == topm[r, Group], rn]

    mgt <- merget %>% data.table::copy()

    mgt[rn %chin% rnl, vcol := V1]

    mgt[is.na(vcol), vcol := "Other"]

    blacklist <- mgt[, .N, by = "vcol"] %>% .[N < 10, vcol]

    mgt[vcol %chin% blacklist, vcol := "Small contributor"]

    dcols <- mgt[, vcol %>% unique]

    names(dcols)[dcols == "Other"] <- "grey"

    names(dcols)[dcols == "Small contributor"] <- "black"

    names(dcols)[dcols == topm[r, V1]] <- "#aa00ff"

    cw_colors2 <- cw_colors[c(1,3:length(cw_colors))]

    names(dcols)[is.na(names(dcols))] <- cw_colors2[1:length(dcols[is.na(names(dcols))])]

    vco <- names(dcols)

    names(vco) <- dcols

    fl <- c(dcols[c("grey", "black", "#aa00ff")], dcols[!names(dcols) %chin% c("grey", "black", "#aa00ff")])

    mgt[, vcolf := factor(vcol, ordered = TRUE, levels = fl)]

    p <- ggplot(mgt, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(col = vcolf), size = 0.1, alpha = 0.4) +
      geom_point(data = mgt[vcol != "grey"], aes(col = vcolf), size = 0.1, alpha = 0.4) +
      facet_wrap(warpf~alphaf, nrow = 2) +
      scale_color_manual(values = vco) +
      ttheme +
      theme_void() +
      theme(legend.position = "right",
      legend.title = element_text(hjust = 0.5),
      legend.text = element_text(size = 6),
      strip.text = element_text(face = "bold", size = 6),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
      guides(color = guide_legend(title = "Top cluster marker",
      ncol = 2, override.aes = list(size = 2, alpha = 1))) +
      ggtitle(paste(tv[names(tv) == sn], "\n", topm[r, V1], " cluster", sep = ""))

    fn <- paste(subdir, sn, "_", topm[r, V1], "_gene_umap.pdf", sep = "")

    ggsave(plot = p, fn, height = 3, width = 4)

    return(p)

  })

  fn <- paste("../Figs/Gene_umaps/", sn, "_gene_umap.pdf", sep = "")

  plot <- cowplot::plot_grid(plotlist = pl, ncol = 1)

  ggsave(plot = plot, fn, width = 6, height = 5 * length(pl), limitsize = FALSE)

})
