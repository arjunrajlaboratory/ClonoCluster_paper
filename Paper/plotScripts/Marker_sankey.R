library(magrittr)
library(data.table)
library(ClonoCluster)
library(ggplot2)

set.seed(42)

source("Paper/extractionScripts/Constants.R")

# loop to get wide table of gene aucs
dtc <- lapply(sn_v, function(sn){

  # get alpha levels
  al <- c(mt[sample_id == sn, alevel2], mt[sample_id == sn, alevel])

  # get hybrid cluster markers
  dt <- data.table::fread(paste(pdir, sn, "_auc_linclust.txt", sep = "")) %>%
  setnames("lin", "Group")

  # get cluster assignments
  dl <- data.table::fread(paste(pdir, sn, "_cluster_assignments.txt", sep = ""))

  # get cluster sizes
  sizet <- dl[, rn %>% length, by = c("Group", "alpha")]

  # looks crappy if you only mark a small cluster on sankey, though valid!
  sc <- 300

  # smaller size for smaller samples
  if (sn %in% c("NJ", "CJ")) sc <- 100

  sizet <- sizet[V1 > sc]

  dt <- merge(dt, sizet, by = c("Group", "alpha"))

  # get max auc per gene for any cluster at that alpha level
  dt[, max_auc := max(auc), by = c("rn", "alpha")]

  dt[alpha == 0, lab := "Transcriptome"]

  dt[alpha == 1, lab := "Barcodes"]

  dt[alpha == al[1], lab := "Low alpha"]

  dt[alpha == al[2], lab := "High alpha"]

  # wide table of max auc at each alpha
  dtc <- dt[!is.na(lab), .SD %>% unique, .SDcols = c("rn", "lab", "max_auc")] %>% dcast(rn ~ lab, value.var = "max_auc")

  dtc[, sample_id := sn]

  return(dtc)

}) %>% data.table::rbindlist(fill = TRUE, use.names = TRUE)

# loop to get genes of interest
dtcl <- lapply(sn_v, function(sn){

  al <- c(mt[sample_id == sn, alevel2], mt[sample_id == sn, alevel])

  d <- dtc[sample_id == sn] %>% data.table::copy()

  d[, gmax := `High alpha` - Transcriptome]

  # Highest delta between high alpha and transcriptome
  g1 <- d[gmax == max(gmax), rn[1]]

  d[, gmax := `Low alpha` - Transcriptome]

  # highest delta between low alpha and transcriptome
  g2 <- d[rn != g1] %>% .[gmax == max(gmax), rn[1]]

  # marker that falls off the most from transcriptome to high alpha
  d[, gmax := Transcriptome - `High alpha`]

  g3 <- d[gmax == max(gmax), rn[1]]

  # marker that is consistently good at all levels
  d[, gmax := `High alpha` + `Low alpha` + Transcriptome]

  g4 <- d[gmax == max(gmax), rn[1]]

  d <- data.table::data.table(sample_id = sn,
    type = c("High @", "Low @", "Transcriptome", "Transcriptome"),
    gl = c(g1,g2,g3, g4),
    al = c(al[2], al[1], 0, 0))

  return(d)

}) %>% data.table::rbindlist()

# add extra genes that I want for later
# dtcl <- list(dtcl,
#   data.table::data.table(sample_id = "Kleind2ws",
#   type = "High @", gl = c("Mpo", "Prtn3"), al = 0.75)
#             ) %>% data.table::rbindlist()


# loop to plot
lapply(sn_v, function(sn){

  print(sn)

  al <- c(mt[sample_id == sn, alevel2], mt[sample_id == sn, alevel])

  lvec <- c(paste("Transcriptome\n", "@", " = 0", sep = ""),
  paste("Low @\n", "@", " = ", al[1], sep = ""),
  paste("High @\n", "@", " = ", al[2], sep = ""),
  paste("Barcodes\n", "@", " = 1", sep = ""))

  dtci <- dtcl[sample_id == sn]

  dt <- data.table::fread(paste(pdir, sn, "_auc_linclust.txt", sep = "")) %>%
  setnames("lin", "Group")

  dl <- data.table::fread(paste(pdir, sn, "_cluster_assignments.txt", sep = ""))

  sizet <- dl[, rn %>% length, by = c("Group", "alpha")]

  # looks crappy if you only mark a small cluster on sankey, though valid!
  sc <- 300

  if (sn %in% c("NJ", "CJ")) sc <- 100

  sizet <- sizet[V1 > sc]

  dt <- merge(dt, sizet, by = c("Group", "alpha"))

  dt[, max_auc := max(auc), by = c("rn", "alpha")]

  dt[alpha == 0, lab := "Transcriptome"]

  dt[alpha == 1, lab := "Barcodes"]

  dt[alpha == al[1], lab := "Low @"]

  dt[alpha == al[2], lab := "High @"]

  # get gene counts for the genes of interest in the sample
  gt <- list.files(dgdir, full.names = TRUE) %>% .[. %like% sn] %>%
    data.table::fread() %>% .[rn %chin% dtci[sample_id == sn, gl]]

  mgt <- gt[, .SD, .SDcols = 2:ncol(gt)] %>% as.matrix

  rownames(mgt) <- gt[, rn]

  mgt %<>% t

  if(sn == "CJ") rownames(mgt) %<>% stringr::str_replace("-[0-9]$", "")

  if(sn == "NJ") rownames(mgt) <- rownames(mgt) %>%
    stringr::str_replace("^S[0-9]_", "") %>%
    stringr::str_replace("-[0-9]$", "")

  # loop over the 4 genes of interest
  pli <- lapply(1:nrow(dtci), function(r){

    # alpha level of interest
    t <- dtci[r, type]

    # gene name
    g <- dtci[r, gl]

    # get cluster at this level
    lin <- dt[rn == g & lab == t & auc == max_auc, Group]

    # get cluster members
    rnl <- dl[alpha == dtci[r, al] & Group == lin, rn]

    # redo ROC so we get the threshold
    c <- ROCR_wrap(mgt[rownames(mgt) %chin% rnl, g],
    mgt[!rownames(mgt) %chin% rnl, g], return_curve = TRUE)

    # find minimum euclidean distance from AUC = 1
    c[, dist := sqrt(((1 - tpr)^2) + ((fpr)^2))]

    c <- c[dist == min(dist)]

    # get cells marked by marker
    v <- mgt[, g] > c[, thresh]

    if (c[, is_flipped]) v <- v == FALSE

    # all marker pos cells at this threshold
    rnl <- rownames(mgt)[v]

    # cells marked by marker and in cluster of interest (true pos)
    rnl2 <- rnl[rnl %chin% dl[alpha == dtci[r, al] & Group == lin, rn]]

    # cells marked but out of cluster (false pos)
    rnlo <- rnl[!rnl %chin% rnl2]

    ttheme <- theme(axis.title = element_text(size = 8, face = "bold", color = "black"),
                    axis.text = element_text(size = 8, color = "black"),
                    axis.text.x = element_text(size = 8, color = "black", angle = 35, hjust = 1),
                    plot.title = element_text(size = 8, face = "bold", color = "black", hjust = 0.5),
                    plot.subtitle = element_text(size = 8, color = "black"))

    # plot
    p <- Plot_alluvia_track(dl[alpha %in% c(0, al)],
                              ids = list(rnlo, rnl2),
                              title = paste(tv[names(tv) == sn], "\n",
                              g, "-positive", sep = ""),
                              ylab = "Number of cells",
                              xlab = NULL,
                              border_size = 0.25,
                              label_nodes = FALSE,
                              flow_alpha = 1,
                              cols = c("gray", ClonoCluster::cw_colors[2]),
                              alluvia_cols = c("gray100", ClonoCluster::cw_colors[2]),
                              col2 = "gray100"
                              )

    # annotation with AUC
    l1 <- dtc[sample_id == sn & rn == g, .(Transcriptome, `Low alpha`, `High alpha`)] %>% unlist

    l1 <- (l1 / 100) %>% sprintf(fmt = "%.2f")

    p <- p + annotate("text", size = 2, x = 1:3, y = nrow(dl[alpha == 0]) * 1.1,
                        label = paste0("AUC: ", l1), fontface = "bold")

    # theming and label x vals
    p <- p + ttheme +
      scale_x_discrete(labels = lvec)

    fn <- paste(fdir, "AUC_sankey/", sn, "_", g, "_auc_sankey.pdf", sep = "")

    ggsave(plot = p, fn, width = 2.75, height = 2.75)

    return(p)

  })

  return(NULL)

  # plot2 <- cowplot::plot_grid(plotlist = pli, nrow = 1, scale = 0.9)
  #
  # ggsave(plot = plot2, paste("../Figs/", sn,"_marker_alluvia.png", sep = ""), height = 2.75, width = 2.75 * length(pli))

})
