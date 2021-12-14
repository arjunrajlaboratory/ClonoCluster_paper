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

fl <- list.files("../Processed_data", pattern = "assignments", full.names = TRUE)

bcs <- list.files("../Data_barcodes", full.names = TRUE)

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

fl <- list.files("../Processed_data", pattern = "auc_linclust", full.names = TRUE)

st <- lapply(fl, function(f){

  dt <- data.table::fread(f)

  dt[, sample_id := basename(f) %>% stringr::str_extract("^[^_]+")]

  dt[, maxauc := max(auc), by = c("lin", "sample_id", "alpha")]

  dt <- dt[auc == maxauc]

  dt %<>% unique(by = c("maxauc", "lin", "alpha"))

  return(dt)

}) %>% data.table::rbindlist()

pll <- lapply(sn_v[c(3,1:2,4:5)], function(sn){

  a <- mt[sample_id == sn, alevel]

  a2 <- mt[sample_id == sn, alevel2]

  # best marker auc
  p4 <- ggplot(st[sample_id == sn & alpha %in% c(0, a2, a, 1)],
      aes(x = as.factor(alpha), y = maxauc)) +
          geom_boxplot(fill = "dodgerblue", size = 0.5) +
          theme_bw() +
          theme(axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(face = "bold", color = "black", size = 8),
          panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_line(size = 1),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
          scale_x_discrete(labels = c("Transcriptome", "Low alpha", "High alpha", "Barcodes")) +
          scale_y_continuous(limits = c(50,100), labels = function(x) {(x / 100) %>% scales::percent()}) +
          xlab(NULL) +
          ylab("Marker AUC\n") +
          ggtitle(paste(tv[names(tv) == sn], "\ntop marker per cluster", sep = ""))

  return(p4)

})

plg <- cowplot::plot_grid(plotlist = pll, ncol = 1, scale = 0.9)

ggsave(plot = plg, "../Figs/Top_markers_boxplot.png", width = 2.5, height = 2.5 * length(pll), limitsize = FALSE)

pll <- lapply(sn_v[c(3,1:2,4:5)], function(sn){

  a <- mt[sample_id == sn, alevel]

  a2 <- mt[sample_id == sn, alevel2]

  # whitelisted cohens_k
  bt <- bcs[bcs %like% sn] %>% data.table::fread()

  wl <- bt[, .N, by = "Barcode"] %>% .[order(-N)] %>%
    .[, Barcode %>% unique %>% .[1:15]]

  c <- ct[sample_id == sn & barcode %chin% wl & alpha %in% c(0,a,a2,1)]

  p2 <- ggplot(c, aes(x = as.factor(alpha), y = cohens_k)) +
          geom_boxplot(fill = "dodgerblue", size = 0.5) +
          theme_bw() +
          theme(axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(face = "bold", color = "black", size = 8),
          panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_line(size = 1),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
          scale_x_discrete(labels = c("Transcriptome", "Low alpha", "High alpha", "Barcodes")) +
          xlab(NULL) +
          ylab("Cohen's kappa\nbarcodes to clusters\n") +
          ggtitle(tv[names(tv) == sn])

  return(p2)

})

ggsave(plot = cowplot::plot_grid(plotlist = pll, ncol = 1, scale = 0.9),
"../Figs/Cohens_k_boxplot.png", width = 2.5, height = 2.5 * length(pll), limitsize = FALSE)

# venn diagrams for marker turnover
library(VennDiagram)

# ligthen color so  plot can be read
hexs <- colorRampPalette(c("#FFFFFF", "dodgerblue"))(10)[c(2,4,6,8)]

# set up our classification colors
names(hexs) <- c("Transcriptome", "Low alpha", "High alpha", "Barcodes")

fl <- list.files("../Processed_data", pattern = "auc_linclust", full.names = TRUE)

vl <- lapply(c(1:5), function(r){

  sn <- mt[r, sample_id]

  a <- mt[r, alevel]

  a2 <- mt[r, alevel2]

  dt <- fl[fl %like% sn] %>% data.table::fread()

  dt <- dt[alpha %in% c(0, a, a2, 1)]

  dt[, max := auc == max(auc), by = c("lin", "alpha")]

  dt <- dt[max == TRUE]

  low <- dt[alpha == 0 & auc > 70, rn %>% unique]

  mid <- dt[alpha == a2 & auc > 70, rn %>% unique]

  high <- dt[alpha == a & auc > 70, rn %>% unique]

  barc <-  dt[alpha == 1 & auc > 70, rn %>% unique]

  # generate initial figure to edit
  venn.plot <- VennDiagram::venn.diagram(list(low, barc, mid, high), filename = NULL,
    category = c("Transcriptome", "Barcodes", "Low\nalpha", "High\nalpha"),
    lwd  = rep(1, 4), lty =  rep("dashed",  4),
    fill = hexs[c(1,4,2,3)], fontface = rep("bold", 15), cex = rep(0.67, 15),
    alpha = rep (0.8, 4), cat.fontface = rep("bold", 4),
    cat.fontfamily = rep("sans", 4), fontfamily = rep("sans", 15),
    cat.cex = rep(0.5, 4),
    cat.dist = c(0.33, 0.33, 0.15, 0.15),
    main.pos = c(0.5, 1.2),
    margin = 0.03, sigdigs = 2, print.mode = c("raw"),
    main = paste(tv[names(tv) == sn], "\ntop cluster markers", sep = ""),
    main.fontfamily = "sans", main.cex = 0.67, main.fontface = "bold")

    grid::grid.newpage()

    grid::grid.draw(venn.plot)

    return(venn.plot)

})

file.remove(list.files(pattern = "VennDiagram.*log$"))

plot <- cowplot::plot_grid(plotlist = vl[c(3, 1:2, 4:5)], scale = 0.8, ncol = 1)

ggsave(plot = plot, "../Figs/Venn_diagrams.png", width = 2.5, height = 2.5 * length(vl))
