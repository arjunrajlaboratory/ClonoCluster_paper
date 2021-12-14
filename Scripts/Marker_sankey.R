library(magrittr)
library(data.table)
library(BarCluster)
library(ggplot2)

sn_v <- c(paste("YG", 1:3, sep = ""), "Kleind2ws", "CJ")

tv <- c("High dose BRAFi - 1", "High dose BRAFi - 2", "Low dose BRAFi", "In Vitro Hematopoiesis Day 2", "Cardio-directed iPSC")

names(tv) <- sn_v

# cutoffs determined by max alpha before cluster number increases
mt <- data.table::data.table(sample_id = sn_v,
  alevel = c(0.55, 0.6, 0.6, 0.75, 0.55))

mt[, alevel2 := (alevel / 2) %>% round(digits = 1)]

# now marker sankeys
dtc <- lapply(sn_v, function(sn){

  al <- c(mt[sample_id == sn, alevel2], mt[sample_id == sn, alevel])

  dt <- data.table::fread(paste("../Processed_data/", sn, "_auc_linclust.txt", sep = "")) %>%
  setnames("lin", "Group")

  dl <- data.table::fread(paste("../Processed_data/", sn, "_cluster_assignments.txt", sep = ""))

  sizet <- dl[, rn %>% length, by = c("Group", "alpha")]

  # looks crappy if you only mark a small cluster on sankey, though valid!
  sc <- 300

  if (sn %in% c("NJ", "CJ")) sc <- 100

  sizet <- sizet[V1 > sc]

  dt <- merge(dt, sizet, by = c("Group", "alpha"))

  dt[, max_auc := max(auc), by = c("rn", "alpha")]

  dt[alpha == 0, lab := "Transcriptome"]

  dt[alpha == 1, lab := "Barcodes"]

  dt[alpha == al[1], lab := "Low alpha"]

  dt[alpha == al[2], lab := "High alpha"]

  dtc <- dt[!is.na(lab), .SD %>% unique, .SDcols = c("rn", "lab", "max_auc")] %>% dcast(rn ~ lab, value.var = "max_auc")

  dtc[, sample_id := sn]

  return(dtc)

}) %>% data.table::rbindlist(fill = TRUE, use.names = TRUE)

dtcl <- lapply(sn_v, function(sn){

  al <- c(mt[sample_id == sn, alevel2], mt[sample_id == sn, alevel])

  d <- dtc[sample_id == sn] %>% data.table::copy()

  d[, gmax := `High alpha` - Transcriptome]

  g1 <- d[gmax == max(gmax), rn]

  d[, gmax := `Low alpha` - Transcriptome]

  g2 <- d[rn != g1] %>% .[gmax == max(gmax), rn]

  d[, gmax := Transcriptome - `High alpha`]

  g3 <- d[gmax == max(gmax), rn]

  d <- data.table::data.table(sample_id = sn,
    type = c("High alpha", "Low alpha", "Transcriptome"),
    gl = c(g1,g2,g3),
    al = c(al[2], al[1], 0))

  return(d)

}) %>% data.table::rbindlist()

lapply(sn_v, function(sn){

  al <- c(mt[sample_id == sn, alevel2], mt[sample_id == sn, alevel])

  dtci <- dtcl[sample_id == sn]

  dt <- data.table::fread(paste("../Processed_data/", sn, "_auc_linclust.txt", sep = "")) %>%
  setnames("lin", "Group")

  dl <- data.table::fread(paste("../Processed_data/", sn, "_cluster_assignments.txt", sep = ""))

  sizet <- dl[, rn %>% length, by = c("Group", "alpha")]

  # looks crappy if you only mark a small cluster on sankey, though valid!
  sc <- 300

  if (sn %in% c("NJ", "CJ")) sc <- 100

  sizet <- sizet[V1 > sc]

  dt <- merge(dt, sizet, by = c("Group", "alpha"))

  dt[, max_auc := max(auc), by = c("rn", "alpha")]

  dt[alpha == 0, lab := "Transcriptome"]

  dt[alpha == 1, lab := "Barcodes"]

  dt[alpha == al[1], lab := "Low alpha"]

  dt[alpha == al[2], lab := "High alpha"]

  gt <- list.files("../Data_genes/", full.names = TRUE) %>% .[. %like% sn] %>%
    data.table::fread() %>% .[rn %chin% dtci[sample_id == sn, gl]]

  mgt <- gt[, .SD, .SDcols = 2:ncol(gt)] %>% as.matrix

  rownames(mgt) <- gt[, rn]

  mgt %<>% t

  if(sn == "CJ") rownames(mgt) %<>% stringr::str_replace("-[0-9]$", "")

  if(sn == "NJ") rownames(mgt) <- rownames(mgt) %>%
    stringr::str_replace("^S[0-9]_", "") %>%
    stringr::str_replace("-[0-9]$", "")


  pli <- lapply(dtci[, type %>% unique], function(t){

    print(t)

    g <- dtci[type == t, gl]

    lin <- dt[rn == g & lab == t & auc == max_auc, Group]

    rnl <- dl[alpha == dtci[type == t, al] & Group == lin, rn]

    c <- ROCR_wrap(mgt[rownames(mgt) %chin% rnl, g],
    mgt[!rownames(mgt) %chin% rnl, g], return_curve = TRUE)

    c[, dist := sqrt(((1 - tpr)^2) + ((fpr)^2))]

    c <- c[dist == min(dist)]

    v <- mgt[, g] > c[, thresh]

    if (c[, is_flipped]) v <- v == FALSE

    rnl <- rownames(mgt)[v]

    rnl <- rnl[rnl %chin% dl[alpha == dtci[type == t, al] & Group == lin, rn]]

    ttheme <- theme(axis.title = element_text(size = 8, face = "bold", color = "black"),
                    axis.text = element_text(size = 8, color = "black"),
                    axis.text.x = element_text(size = 8, color = "black", angle = 35, hjust = 1),
                    plot.title = element_text(size = 8, face = "bold", color = "black", hjust = 0.5),
                    plot.subtitle = element_text(size = 8, color = "black"))

    p <- Plot_alluvia_track(dl[alpha %in% c(0, al)],
                              ids = rnl,
                              title = paste(tv[names(tv) == sn], "\n",
                              g, " marker AUC", sep = ""),
                              ylab = "Number of cells",
                              xlab = NULL,
                              border_size = 0.25,
                              label_nodes = FALSE,
                              cols = BarCluster::cw_colors[2]
                              )

    #annotation
    l1 <- dtc[sample_id == sn & rn == g, .(Transcriptome, `Low alpha`, `High alpha`)] %>% unlist

    l1 <- (l1 / 100) %>% sprintf(fmt = "%.2f")

    p <- p + annotate("text", size = 2.5, x = 1:3, y = nrow(dl[alpha == 0]) * 1.1,
                        label = paste0(l1), fontface = "bold")

    p <- p + ttheme +
      scale_x_discrete(labels = c("Transcriptome", "Low alpha", "High alpha"))

    fn <- paste("../Figs/", sn, "_", g, "_auc_sankey.png", sep = "")

    ggsave(plot = p, fn, width = 2.75, height = 2.75)

    return(p)

  })

  plot2 <- cowplot::plot_grid(plotlist = pli, nrow = 1, scale = 0.9)

  ggsave(plot = plot2, paste("../Figs/", sn,"_marker_alluvia.png", sep = ""), height = 2.75, width = 2.75 * length(pli))

})
