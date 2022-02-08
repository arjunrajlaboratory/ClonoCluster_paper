library(magrittr)
library(data.table)
library(ClonoCluster)
library(ggplot2)

set.seed(42)

source("Paper/extractionScripts/Constants.R")

# hybrid cluster assignments
fl <- list.files(pdir, pattern = "auc_linclust", full.names = TRUE)

# barcode assignments
bcs <- list.files(dbdir, full.names = TRUE)

# get top marker per cluster
st <- lapply(fl, function(f){

  dt <- data.table::fread(f)

  dt[, sample_id := basename(f) %>% stringr::str_extract("^[^_]+")]

  dt[, maxauc := max(auc), by = c("lin", "sample_id", "alpha")]

  dt <- dt[auc == maxauc]

  dt %<>% unique(by = c("maxauc", "lin", "alpha"))

  return(dt)

}) %>% data.table::rbindlist()

# stats
# global tests
pt <- lapply(sn_v, function(sn){

  a <- mt[sample_id == sn, alevel]

  a2 <- mt[sample_id == sn, alevel2]

  a <- st[sample_id == sn & alpha %in% c(0, a2, a, 1)]

  pt <- kruskal.test(a$maxauc, a$alpha)

  return(data.table::data.table(sample_id = sn, p = pt$p.value))

}) %>% data.table::rbindlist()

pt[, pcorr := p %>% p.adjust(method = "bonferroni")]

# samples that passed global KW test
wl <- pt[pcorr < 0.05, sample_id]

# print pairwise wilcox tests against transcriptome level
lapply(wl, function(sn){

  a <- mt[sample_id == sn, alevel]

  a2 <- mt[sample_id == sn, alevel2]

  a <- st[sample_id == sn & alpha %in% c(0, a2, a, 1)]

  pm <- pairwise.wilcox.test(a$maxauc, a$alpha, p.adjust.method = "bonferroni")$p.value

  return(data.table::data.table(sn, rownames(pm), pm[, 1]))

}) %>% data.table::rbindlist() %>% .[V3 < 0.05]

# loop over samples for boxplots
pll <- lapply(sn_v, function(sn){

  a <- mt[sample_id == sn, alevel]

  a2 <- mt[sample_id == sn, alevel2]

  lvec <- c(paste("Transcriptome\n", "@", " = 0", sep = ""),
  paste("Low @\n", "@", " = ", a2, sep = ""),
  paste("High @\n", "@", " = ", a, sep = ""),
  paste("Barcodes\n", "@", " = 1", sep = ""))

  # best top markers boxplot
  p4 <- ggplot(st[sample_id == sn & alpha %in% c(0, a2, a, 1)],
      aes(x = as.factor(alpha), y = maxauc)) +
          geom_boxplot(fill = "dodgerblue", size = 0.5) +
          theme_bw() +
          theme(axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(face = "bold", color = "black", size = 8),
          panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_line(size = 0.5),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
          scale_x_discrete(labels = lvec) +
          scale_y_continuous(limits = c(50,100), labels = function(x) {(x / 100) %>% scales::percent()}) +
          xlab(NULL) +
          ylab("Marker AUC\n") +
          ggtitle(paste(tv[names(tv) == sn], "\ntop cluster markers", sep = ""))

  fn <- paste(fdir, "top_markers_box/", sn, "_top_markers_boxplot.pdf", sep = "")

  ggsave(plot = p4, fn, width = 2.5, height = 2.5, limitsize = FALSE)

  return(p4)

})

# venn diagrams for marker turnover
library(VennDiagram)

# ligthen color so  plot can be read
hexs <- colorRampPalette(c("#FFFFFF", "dodgerblue"))(10)[c(2,4,6,8)]

# set up our classification colors
names(hexs) <- c("Transcriptome", "Low @", "High @", "Barcodes")

fl <- list.files(pdir, pattern = "auc_linclust", full.names = TRUE)

# venn of all cluster markers AUC > 70
vl <- lapply(c(1:10), function(r){

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

  lvec <- c(paste("Transcriptome\n", "@", " = 0", sep = ""),
  paste("Low @\n", "@", " = ", a2, sep = ""),
  paste("High @\n", "@", " = ", a, sep = ""),
  paste("Barcodes\n", "@", " = 1", sep = ""))

  # generate initial figure to edit
  venn.plot <- VennDiagram::venn.diagram(list(low, barc, mid, high), filename = NULL,
    category = lvec[c(1,4,2:3)],
    lwd  = rep(1, 4), lty =  rep("dashed",  4),
    fill = hexs[c(1,4,2,3)], fontface = rep("bold", 15), cex = rep(0.67, 15),
    alpha = rep (0.8, 4), cat.fontface = rep("bold", 4),
    cat.fontfamily = rep("sans", 4), fontfamily = rep("sans", 15),
    cat.cex = rep(0.5, 4),
    cat.dist = c(0.33, 0.33, 0.15, 0.15),
    main.pos = c(0.5, 1.15),
    margin = 0.03, sigdigs = 2, print.mode = c("raw"),
    main = paste(tv[names(tv) == sn], "\ntop cluster markers", sep = ""),
    main.fontfamily = "sans", main.cex = 0.67, main.fontface = "bold")

    grid::grid.newpage()

    grid::grid.draw(venn.plot)

    fn <- paste(fdir, "Venn/", sn, "venn.pdf",sep = "")

    ggsave(plot = venn.plot, fn, width = 3, height = 2.5)

    return(venn.plot)

})

# get rid of log files
file.remove(list.files(pattern = "VennDiagram.*log$"))

plot <- cowplot::plot_grid(plotlist = vl, scale = 0.8, ncol = 1)

ggsave(plot = plot, paste(fdir, "Venn/Venn_diagrams.pdf", sep = ""), width = 3, height = 2.5 * length(vl))
