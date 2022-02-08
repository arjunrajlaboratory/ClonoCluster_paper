library(data.table)
library(magrittr)
library(ClonoCluster)
library(ggplot2)

set.seed(42)

# make 3 normal distributions
b <- rnorm(n = 1000, mean = -3, sd = 3)

a <- rnorm(n = 1000, mean = 6, sd = 5)

c <- rnorm(n = 1000, mean = 0, sd = 4)

# 3000 cells per barcode
bt <- lapply(c("A", "B", "C"), function(x) replicate(1000, x)) %>% unlist

bt <- data.table::data.table(rn = 1:length(bt),  Barcode = bt)

# make 500 cells barcode D but will get same PCs as C
bt[rn %in% 2501:3000, Barcode := "D"]

# Make some of barcode C B actually
bt[rn %in% 2300:2500, Barcode := "B"]

# make 6 PCs
pc1 <- c(a,b,c)

pc2 <- c(a,b,b)

pc3 <- c(c,a,b[1:500], c[1:250], a[1:250])

pc4 <- c(c,a,a)

pc5 <- c(b,a,a)

pc6 <- c(b,c,c)

# fake PC matrix
irl <- cbind(pc1,pc2,pc3,pc4,pc5,pc6)

rownames(irl) <- bt[, rn]

# warp the PCs
wfs <- c(0, 2.5, 5, 7.5, 10)

irl_array <- lapply(wfs, function(a){

  mo <- barcode_warp(irl, bt, a)

  return(mo)

})

# loop over warped PCAs and umap
umaps <- lapply(irl_array %>% seq_along, function(i){

  um <- umap_matrix(irl_array[[i]])

  um[, PC1 := irl_array[[i]][, "pc1"]]

  um[, PC2 := irl_array[[i]][, "pc3"]]

  um[, wf := wfs[i]]

  return(um)

}) %>% data.table::rbindlist()

umaps[, rn := as.numeric(rn)]

umaps <- merge(umaps, bt, by = "rn")

umaps %<>% melt(id.vars = names(umaps)[!names(umaps) %like% "^PC"])

# get mean PC per barcode
umaps[, linex := mean(value), by = c("Barcode", "variable")]

# plot UMAP facets
up <- ggplot(umaps, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(col = Barcode), alpha = 0.1, size = 0.001) +
  facet_wrap(~wf, nrow = 1) +
  theme_void() +
  theme(legend.title = element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 6),
      panel.spacing = unit(2, "lines")) +
  scale_color_manual(values = cw_colors[c(2:4,1)], labels = paste("Barcode", 1:4)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

fn <- paste(fdir, "simulations/", "Fig4_umap.pdf", sep = "")

ggsave(plot = up, filename = fn, height = 6/5, width = 6)

# Plot PC histograms
pcp <- ggplot(umaps, aes(x = value)) +
  geom_histogram(aes(fill = Barcode), bins = 100) +
  geom_vline(aes(xintercept = linex, color = Barcode), linetype = "dashed", size = 0.5) +
  facet_wrap(variable ~ wf, nrow = 2, scales = "free_y") +
  theme_bw() +
  theme(legend.title = element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 6),
        legend.position = "none") +
  scale_color_manual(values = cw_colors[c(2:4,1)], labels = paste("Barcode", 1:4), guide = "none") +
  scale_fill_manual(values = cw_colors[c(2:4,1)], labels = paste("Barcode", 1:4)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 2)))

fn <- paste(fdir, "simulations/", "Fig4_pcs.pdf", sep = "")

ggsave(plot = pcp, filename = fn, height = 12/5, width = 6)

# get variation within vs between all barcodes
vart <- lapply(irl_array %>% seq_along, function(wf){

  irl <- irl_array[[wf]]

  d <- cbind(bt, irl)

  d %<>% melt(id.vars = c("rn", "Barcode"))

  varw <- d[, sd(value) ^ 2, by = c("Barcode", "variable")]

  varb <- d[, sd(value) ^ 2, by = "variable"]

  varb[, variation := "between"]

  varb[, wf := wfs[[wf]]]

  varb[, Barcode := NA]

  varw[, wf := wfs[[wf]]]

  varw[, variation := "within"]

  vart <- rbind(varw, varb)

  return(vart)

}) %>% data.table::rbindlist()

# normalize to value at warp 0
normt <- vart[wf == 0, .(variable, Barcode, variation, V1)] %>% setnames("V1", "init")

vart <- merge(vart, normt, by = c("variable", "Barcode", "variation"))

vart[, rel := V1 / init]

vart[, mean_rel := mean(rel), by = c("wf", "variation")]

ttheme <- theme(axis.title = element_text(size = 8, face = "bold", color = "black"),
                axis.text = element_text(size = 8, color = "black"),
                plot.title = element_text(size = 8, face = "bold", color = "black", hjust = 0.5),
                plot.subtitle = element_text(size = 8, color = "black"))

# plot variance within vs between
p <- ggplot(vart, aes(x = wf, y = mean_rel)) +
  geom_point(aes(col = variation)) +
  geom_line(aes(col = variation)) +
  theme_bw() +
  ttheme +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#D55E00", "#009E73") %>% rev) +
  xlab("\nWarp factor") +
  ylab("% of initial variance\n")

fn <- paste(fdir, "simulations/", "Variance_plot.pdf", sep = "")

ggsave(plot = p, filename = fn, height = 2.5, width = 2.5)
