library(magrittr)
library(data.table)
library(ClonoCluster)
library(ggplot2)

set.seed(42)

source("Paper/extractionScripts/Constants.R")

# reorg sankeys, get table of our clusters of interest
dtc <- lapply(sn_v, function(sn){

  al <- c(mt[sample_id == sn, alevel2], mt[sample_id == sn, alevel])

  # get reorg markers
  dt <- data.table::fread(paste(pdir, sn, "_auc_barcodes.txt", sep = "")) %>%
  setnames("lin", "Group")

  # get cluster assignments
  dl <- data.table::fread(paste(pdir, sn, "_cluster_assignments.txt", sep = ""))

  # sort by size to get easily visualized representative plots with small figures
  sizet <- dl[, rn %>% length, by = c("Group", "alpha")]

  # looks crappy if you only mark a small cluster on sankey, though valid!
  sc <- 300

  if (sn %in% c("NJ", "CJ")) sc <- 100

  sizet <- sizet[V1 > sc]

  dt <- merge(dt, sizet, by = c("Group", "alpha"))

  # get contributing transcriptome size
  trt <- dl[alpha == 0, rn %>% length, by = c("Group")] %>% setnames(c("Group", "V1"), c("tr", "tr_size"))

  dt <- merge(dt, trt, by = c("tr"))

  dt[, diff := V1 - tr_size]

  # dont want to take clusters that are totally eaten, want to show classification for this figure
  dt <- dt[diff > 10]

  # low alpha level for comparison
  dt <- dt[alpha == al[1]]

  # get max auc for gene at that level
  dt[, max_auc := max(auc), by = c("rn", "alpha")]

  dt <- dt[auc == max_auc]

  dt[, sample_id := sn]

  return(dt)

}) %>% data.table::rbindlist(fill = TRUE, use.names = TRUE)

# # Basic function to convert mouse to human gene names
# convertMGL <- function(x){
#
#   require("biomaRt")
#
#   human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#
#   mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#
#   genes <- biomaRt::getLDS(attributes = c("mgi_symbol"),
#                   filters = "mgi_symbol",
#                   values = x ,
#                   mart = mouse,
#                   attributesL = c("hgnc_symbol"),
#                   martL = human,
#                   uniqueRows = TRUE)
#
#   conv <- genes %>% unique %>% as.data.table
#
#   return(conv)
#
# }
#
# convertHGL <- function(x){
#
#   require("biomaRt")
#
#   human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#
#   mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#
#   genes <- biomaRt::getLDS(attributes = c("hgnc_symbol"),
#                   filters = "hgnc_symbol",
#                   values = x,
#                   mart = human,
#                   attributesL = c("mgi_symbol"),
#                   martL = mouse,
#                   uniqueRows = TRUE)
#
#   conv <- genes %>% unique %>% as.data.table
#
#   return(conv)
#
# }

# get reorger with best auc for any group
dtc[, max_auc := max(auc), by = c("sample_id", "Group", "tr")]

dtc <- dtc[auc == max_auc]

# take top 10 representative from each sample
dtc %<>% .[order(-auc)] %>% split(by = "sample_id") %>% lapply(., function(d){
  return(d[1:10])
}) %>% data.table::rbindlist()

dtc %<>% na.omit()

dtc <- dtc[rn %chin% c("CCNB1", "PCOLCE")]

# loop over samples to make plots
plotl <- lapply(sn_v, function(sn){

  print(sn)

  al <- c(mt[sample_id == sn, alevel2], mt[sample_id == sn, alevel])

  alv <- c(0, al, 1)

  # set up labels
  names(alv) <- c("Transcriptome", "Low @", "High @", "Barcodes")

  lvec <- c(paste("Transcriptome\n", "@", " = 0", sep = ""),
  paste("Low @\n", "@", " = ", al[1], sep = ""),
  paste("High @\n", "@", " = ", al[2], sep = ""),
  paste("Barcodes\n", "@", " = 1", sep = ""))

  dtci <- dtc[sample_id == sn]

  if (nrow(dtci) == 0) return(NULL)

  # get cluster assignments
  dl <- data.table::fread(paste(pdir, sn, "_cluster_assignments.txt", sep = ""))

  # get relevant gene counts
  gt <- list.files(dgdir, full.names = TRUE) %>% .[. %like% sn] %>%
    data.table::fread() %>% .[rn %chin% dtci[sample_id == sn, rn]]

  mgt <- gt[, .SD, .SDcols = 2:ncol(gt)] %>% as.matrix

  rownames(mgt) <- gt[, rn]

  mgt %<>% t

  # conditionals for sample cell IDs
  if(sn == "CJ") rownames(mgt) %<>% stringr::str_replace("-[0-9]$", "")

  if(sn == "NJ") rownames(mgt) <- rownames(mgt) %>%
    stringr::str_replace("^S[0-9]_", "") %>%
    stringr::str_replace("-[0-9]$", "")

  # loop over representative clusters and genes
  pli <- lapply(1:nrow(dtci), function(r){

    print(r)

    # gene
    g <- dtci[r, rn]

    # hybrid clusters
    lin <- dtci[r, Group]

    # contributing transcriptome cluster
    tr <- dtci[r, tr]

    # alpha value of low alpha
    ali <- alv[names(alv) == "Low @"]

    # all members of low alpha cluster
    lac <- dl[alpha == ali & Group == lin, rn]

    # all members of contributing transcriptome cluster
    trc <- dl[alpha == 0 & Group == tr, rn]

    # ingroup is cells from contributing transcriptome clusters in low alpha cluster
    ing <- intersect(lac,trc)

    # outgroup is transcriptome cluster that doesnt go to low alpha cluster
    outg <- setdiff(trc,lac)

    # this shouldnt trigger, will error if it does
    if (length(outg) == 0) {

      rnlp <- ing

      return(g)

    } else {

      print("roc")

      # rerun ROC to get threshold
      c <- ROCR_wrap(mgt[rownames(mgt) %chin% ing, g],
      mgt[rownames(mgt) %chin% outg, g], return_curve = TRUE)

      # get AUC
      auc_c <- c[, auc] %>% unique

      # determine threshold by euclidean distance
      c[, dist := sqrt(((1 - tpr)^2) + ((fpr)^2))]

      c <- c[dist == min(dist)]

      v <- mgt[, g] > c[, thresh]

      if (c[, is_flipped]) v <- v == FALSE

      # all marker positive
      rnlp <- rownames(mgt)[v]

      # all marker negatives
      rnln <- rownames(mgt)[!v]

      # positives in contributing transcriptome cluster
      rnlp <- rnlp[rnlp %chin% trc]

      # negatives in contributing transcriptome cluster
      rnln <- rnln[rnln %chin% trc]

      # cells from transcriptome cluster that are negative but in low alpha cluster, false negative
      unmarked <- trc[!trc %chin% rnlp & trc %chin% lac]

    }

    # marked and from contributing transcriptome cluster to low alpha cluster
    # true positive
    tp <- intersect(rnlp, ing)

    # Marked and from transcriptome cluster but not in low alpha cluster
    # false positive
    fp <- rnlp[!rnlp %chin% ing]

    # unmarked and does go to low alpha cluster
    # true negative
    tn <- intersect(outg, rnln)

    # unmarked but goes to low alpha cluster
    # false negative
    fn <- rnln[!rnln %chin% outg]

    return(list(tp,fp, unmarked, auc_c, tn, fn))

  })

  # add gene names to list
  names(pli) <- dtci[, rn]

  lengthv <- lapply(pli, length) %>% unlist

  # get rid of null results (clusters that are subsets and won't make good figures)
  pli <- pli[lengthv != 1]

  # true pos
  ids <- pli %>% lapply(., function(d){

    return(d[[1]] %>% unique)

  })

  # false pos
  fps <- pli %>% lapply(., function(d){

    return(d[[2]] %>% unique)

  })

  # true negative
  tns <- pli %>% lapply(., function(d){

    return(d[[5]] %>% unique)

  })

  # false neg
  fns <- pli %>% lapply(., function(d){

    return(d[[6]] %>% unique)

  })

  # unmarked, not used
  ums <- pli %>% lapply(., function(d){

    return(d[[3]] %>% unique)

  }) %>% unlist %>% unique

  # auc
  auc_cs <- pli %>% lapply(., function(d){

    return(d[[4]])

  }) %>% unlist

  names(ids) <- names(pli)

  ttheme <- theme(axis.title = element_text(size = 8, face = "bold", color = "black"),
                  axis.text = element_text(size = 8, color = "black"),
                  axis.text.x = element_text(size = 8, color = "black", angle = 35, hjust = 1),
                  plot.title = element_text(size = 8, face = "bold", color = "black", hjust = 0.5),
                  plot.subtitle = element_text(size = 8, color = "black"))

    colv <- c("grey", ClonoCluster::cw_colors)

    cn <- dl[alpha == 0, rn %>% unique %>% length]

    # plot each gene
    cpl <- lapply(ids %>% seq_along, function(n){

      print(n)

      g <- names(ids)[n]

      lac <- dl[alpha == al[1] & Group == dtci[rn == g, Group %>% unique], rn]

      trc <- dl[alpha == 0 & Group == dtci[rn == g, tr %>% unique], rn]

      # vector of colors (colorblind friendly)
      colv <- c("#CC79A7", "#0072B2", "#D55E00", "#009E73")

      test <- lapply(list(fns[[n]], tns[[n]], fps[[n]], ids[[n]]), length) %>% unlist

      colv <- colv[which(test != 0)]

      p <- Plot_alluvia_track(dl[alpha %in% c(0, al[1])],
                                ids = list(fns[[n]], tns[[n]], fps[[n]], ids[[n]]),
                                title = paste(tv[names(tv) == sn], "\n",
                                "@ = ", al[1], " Cluster ",
                                dtci[, Group %>% unique], ": ", names(ids[n]), sep = ""),
                                ylab = paste(cn, "total cells"),
                                xlab = NULL,
                                border_size = 0.1,
                                label_nodes = FALSE,
                                flow_alpha = 1,
                                cols = colv,
                                alluvia_cols = replicate(n = length(colv), "grey"),
                                col2 = "gray100",
                                orientation = "bottom"
                                )

      p <- p + ttheme +
        scale_x_discrete(labels = lvec)

      # get AUCs for annotation
      lac_auc <- ROCR_wrap(mgt[rownames(mgt) %chin% lac, g],
      mgt[!rownames(mgt) %chin% lac, g])

      trc_auc <- ROCR_wrap(mgt[rownames(mgt) %chin% trc, g],
      mgt[!rownames(mgt) %chin% trc, g])

      # set up labels
      l1 <- c(trc_auc, auc_cs[[n]], lac_auc)

      l1 %<>% as.numeric %>% `/`(., 100) %>% sprintf(fmt = "%.2f")

      l1 <- paste(c("Transcriptome\nAUC: ", "Reorganization\nAUC: ", "Low @\nAUC: "), l1)

      p <- p + annotate("text", size = 2, x = c(1, 1.5, 2), y = cn * 1.1,
                          label = l1, fontface = "bold")

      labv <- c("Unmarked", "True Positive", "False Positive")

      p <- p + theme(legend.position = "none", legend.key.size = unit(0.1, "in"),
                    legend.text = element_text(size = 6),
                    axis.ticks = element_blank(), axis.text = element_blank()) +
        # scale_fill_manual(labels = labv, values = colv) +
        # scale_color_manual(labels = labv, values = colv) +
        guides(fill = "none", color = guide_legend(title = NULL))

      #
      fn <- paste(fdir, "Reorg_sankey/", sn, "_", g, "_sample_cluster_sankey.pdf", sep = "")

      ggsave(plot = p, fn, width = 3.5, height = 3.5)

      return(p)

    })

    return(NULL)

})
