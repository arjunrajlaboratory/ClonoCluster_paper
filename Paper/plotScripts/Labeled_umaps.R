library(magrittr)
library(data.table)
library(ClonoCluster)
library(ggplot2)

set.seed(42)

source("Paper/extractionScripts/Constants.R")

# discrete color vector
dcols <- c("#629B60", "#98B3DF", "#E7DFCC", "#56E884", "#EDE644", "#C2ECA5",
"#D44DA1", "#E68D9E", "#D4F0E6", "#D6D5E6", "#546790", "#5C99E6",
"#9BDFC0", "#E4F2C2", "#7E37D9", "#64EC41", "#5EEABC", "#6465CF",
"#DBAC48", "#A7B73A", "#E8B8E1", "#E4E97F", "#DFCF89", "#5EE7E1",
"#B691A5", "#DFB59C", "#52BEE4", "#E74B64", "#78A59C", "#A797E0",
"#A6E770", "#96DAE6", "#A6AB79", "#E48449", "#E050E6", "#DF90DB",
"#B47EE6", "#98E699", "#9C7058", "#B5EA3B")

lapply(sn_v, function(sn){

  set.seed(42)

  # community detection resolution
  res <- 1

  # set lower for large Weinreb et al. data
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

  # alpha levels of interest, 0 and high alpha
  als <- c(0, mt[sample_id == sn, alevel] %>% unlist)

  # choose barcode matrix method to be fastest empirically
  if (!sn %chin% c("CJ", "NJ", paste("MDA", 1:2, sep = ""), paste("WM983", 1:2, sep = ""))){

    method <- "index"

  } else {

    method <- "fast"

  }

  # get hybrid clusters
  dl <- clonocluster(irl, bt, alpha = als, res = res, method = method)

  # warp pca at these levels
  gl <- lapply(c(0,5), function(w){

    mo <- engage_warp(irl, bt, w)

    return(mo)

  }) %>% data.table::rbindlist()

  # merge warped umaps to barcodes
  gl <- merge(gl, bt, by = "rn")

  # get barcode size
  gl[, bcs := rn %>% unique %>% length, by = "Barcode"]

  # mark singlets
  gl[bcs == 1, Barcode := "Singlet"]

  # get markers for clusters
  marks <- Find_Markers_ROC(dl, gt, n_threads = 6)

  # get top cluster per gene
  topm <- marks[, .(rn[auc == max(auc)], auc[auc == max(auc)]), by = c("alpha", "Group")]

  # get top marker for each hyrbid cluster
  topm[, mauc := max(V2), by = c("V1", "alpha")]

  topm <- topm[mauc == V2]

  # get the high alpha groups of interest and its marker
  grt <- topm[alpha == als[2], .SD %>% unique, .SDcols = c("V1", "Group")]

  topm %<>% dcast(V1 ~ alpha, value.var = "V2")

  topm %>% setnames(as.character(als), c("Transcriptome", "High @"))

  # only markers that weren't present at transcriptome level
  topm <- topm[is.na(Transcriptome)]

  # merge table of top marker names to group numbers
  topm <- merge(topm, grt, by = "V1")

  topm %<>% .[order(-`High @`)]

  # only take markers > 80
  topm <- topm[`High @` > 80]

  # merge top markers to long format cluster assignments
  merget <- merge(gl, dl, by = "rn", allow.cartesian = TRUE)

  # make factors so orders look right and facets are ordered
  merget %<>% .[order(alpha, warp)]

  merget[, alphaf := paste("@ =", alpha)]

  merget[, warpf := paste("Warp ", warp)]

  merget[, alphaf := factor(alphaf, ordered = TRUE, levels = paste("@ =", als))]

  merget[, warpf := factor(warpf, ordered = TRUE, levels = paste("Warp ", c(0,5)))]

  merget %>% setkey(alphaf, warpf)

  topmt <- marks[, .(rn[auc == max(auc)], auc[auc == max(auc)]), by = c("alpha", "Group")]

  merget <- merge(merget, topmt, by = c("alpha", "Group"))

  subdir <- paste(fdir, "Gene_umaps/", sn, "/", sep = "")

  if(!dir.exists(subdir)) dir.create(subdir)

  # loop over selected markers
  pl <- lapply(1:nrow(topm), function(r){

    # get cluster members
    rnl <- merget[alpha == als[2] & Group == topm[r, Group], rn]

    mgt <- merget %>% data.table::copy()

    mgt[rn %chin% rnl, vcol := V1]

    mgt[is.na(vcol), vcol := "Other"]

    # find small contributing clusters and lump together
    blacklist <- mgt[, .N, by = "vcol"] %>% .[N < 10, vcol]

    mgt[vcol %chin% blacklist, vcol := "Small contributor"]

    # get colors right
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

    # plot umap faceted by alpha and warp, colored by cluster top markers if part of high α cluster of interest
    p <- ggplot(mgt, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(col = vcolf), size = 0.1, alpha = 0.4) +
      geom_point(data = mgt[vcol != "Other"], aes(col = vcolf), size = 0.1, alpha = 0.4) +
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

  return(NULL)
  # fn <- paste("../Figs/Gene_umaps/", sn, "_gene_umap.pdf", sep = "")
  #
  # plot <- cowplot::plot_grid(plotlist = pl, ncol = 1)
  #
  # ggsave(plot = plot, fn, width = 6, height = 5 * length(pl), limitsize = FALSE)

})


###

# sankey of Mpo and Prtn3 specifically and increase umap warp for Kleind2ws

sn <- "Kleind2ws"

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

als <- c(0, mt[sample_id == sn, alevel] %>% unlist)

if (!sn %chin% c("CJ", "NJ", paste("MDA", 1:2, sep = ""), paste("WM983", 1:2, sep = ""))){

  method <- "index"

} else {

  method <- "fast"

}

dl <- clonocluster(irl, bt, alpha = als, res = res, method = method)

# warp pca at these levels
gl <- lapply(c(0,9.9), function(w){

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

merget[, warpf := factor(warpf, ordered = TRUE, levels = paste("Warp ", c(0,9.9)))]

merget %>% setkey(alphaf, warpf)

topmt <- marks[, .(rn[auc == max(auc)], auc[auc == max(auc)]), by = c("alpha", "Group")]

merget <- merge(merget, topmt, by = c("alpha", "Group"))

subdir <- paste(fdir, "Gene_umaps/", sn, "/", sep = "")

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
    geom_point(data = mgt[vcol != "Other"], aes(col = vcolf), size = 0.1, alpha = 0.4) +
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

# make marker sankeys
# get our marker table again
topm <- marks[, .(rn[auc == max(auc)], auc[auc == max(auc)]), by = c("alpha", "Group")]

topm[, mauc := max(V2), by = c("V1", "alpha")]

topm <- topm[mauc == V2]

# top cluster markers table
grt <- topm[, .SD %>% unique, .SDcols = c("V1", "Group", "alpha")]

grt %>% setnames("V1", "gname")

dl <- merge(dl, grt, by = c("Group", "alpha"))

# replacing group number with group top cluster marker
dl %>% setnames(c("Group", "gname"), c("gname", "Group"))

dl %<>% .[, .SD, .SDcols = c("rn", "Group", "alpha")]

# genes of interest
#         V1 Transcriptome   High @ Group
# 1: Ccdc62            NA 83.79476     8
# 2:   Cpa3            NA 95.20596     4
# 3:    F2r            NA 82.10204     5
# 4:   Lyz2            NA 82.45051     3
# 5:    Mpo            NA 84.45219     2
# 6:  Prtn3            NA 83.28624     6
gv <- c("Mpo", "Prtn3")

# loop over genes
pli <- lapply(gv %>% seq_along, function(n){

  g <- gv[n]

  # get high alpha cluster members
  rnlc <- dl[Group == g & alpha == als[2], rn]

  # run ROC and get threshold
  c <- ROCR_wrap(gt[rownames(gt) %chin% rnlc, g],
  gt[!rownames(gt) %chin% rnlc, g], return_curve = TRUE)

  # find minimum euclidean distance from AUC = 1
  c[, dist := sqrt(((1 - tpr)^2) + ((fpr)^2))]

  c <- c[dist == min(dist)]

  # get cells marked by marker
  v <- gt[, g] > c[, thresh]

  if (c[, is_flipped]) v <- v == FALSE

  # all marker pos cells at this threshold
  rnl <- rownames(gt)[v]

  # cells marked by marker and in cluster of interest (true pos)
  rnl2 <- rnl[rnl %chin% rnlc]

  # cells marked but out of cluster (false pos)
  rnlo <- rnl[!rnl %chin% rnl2]

  ttheme <- theme(axis.title = element_text(size = 8, face = "bold", color = "black"),
                  axis.text = element_text(size = 8, color = "black"),
                  axis.text.x = element_text(size = 8, color = "black", angle = 35, hjust = 1),
                  plot.title = element_text(size = 8, face = "bold", color = "black", hjust = 0.5),
                  plot.subtitle = element_text(size = 8, color = "black"))

  # plot
  p <- Plot_alluvia_track(dl,
                          ids = list(rnlo, rnl2),
                            title = paste("Contributing transcriptome clusters to\n", g, " high @ cluster", sep = ""),
                            ylab = "Number of cells",
                            xlab = NULL,
                            border_size = 0.5,
                            label_nodes = TRUE,
                            flow_alpha = 1,
                            cols = c("gray", ClonoCluster::cw_colors[2]),
                            alluvia_cols = c("gray100", ClonoCluster::cw_colors[2]),
                            col2 = "gray100",
                            label_size = 1.75,
                            )

  # x axis labels
  lvec <- c(paste("Transcriptome\n", "@", " = 0", sep = ""),
  paste("High @\n", "@", " = ", als[2], sep = ""))

  # theming and label x vals
  p <- p + ttheme +
    scale_x_discrete(labels = lvec)

  fn <- paste(fdir, "Gene_umaps/Kleind2ws/", sn, "_", g, "_auc_sankey.pdf", sep = "")

  ggsave(plot = p, fn, width = 2.75, height = 2.75)

  return(p)

})
