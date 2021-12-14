library(data.table)
library(magrittr)
library(BarCluster)
library(ggplot2)

set.seed(42)

b <- rnorm(n = 1000, mean = -3, sd = 3)

a <- rnorm(n = 1000, mean = 6, sd = 5)

c <- rnorm(n = 1000, mean = 0, sd = 4)

bt <- lapply(c("A", "B", "C"), function(x) replicate(1000, x)) %>% unlist

bt <- data.table::data.table(rn = 1:length(bt),  Barcode = bt)

bt[rn %in% 2501:3000, Barcode := "D"]

bt[rn %in% 2300:2500, Barcode := "B"]

pc1 <- c(a,b,c)

pc2 <- c(a,b,b)

pc3 <- c(c,a,b)

pc4 <- c(c,a,a)

pc5 <- c(b,a,a)

pc6 <- c(b,c,c)

irl <- cbind(pc1,pc2,pc3,pc4,pc5,pc6)

rownames(irl) <- bt[, rn]

alphas <- c(0, 0.25, 0.5, 0.75, 1)

irl_array <- lapply(alphas, function(a){

  mo <- barcode_warp(irl, bt, a)

  return(mo)

})

umaps <- lapply(irl_array %>% seq_along, function(i){

  um <- umap_matrix(irl_array[[i]])

  um[, PC1 := irl_array[[i]][, "pc1"]]

  um[, PC2 := irl_array[[i]][, "pc3"]]

  um[, alpha := alphas[i]]

  return(um)

}) %>% data.table::rbindlist()

umaps[, rn := as.numeric(rn)]

umaps <- merge(umaps, bt, by = "rn")

umaps %<>% melt(id.vars = names(umaps)[!names(umaps) %like% "^PC"])

umaps[, linex := mean(value), by = c("Barcode", "variable")]

up <- ggplot(umaps, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(col = Barcode), alpha = 0.1, size = 0.001) +
  facet_wrap(~alpha, nrow = 1) +
  theme_void() +
  theme(legend.title = element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 6),
      panel.spacing = unit(2, "lines")) +
  scale_color_manual(values = cw_colors[c(2:4,1)], labels = paste("Barcode", 1:4)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

ggsave(plot = up, "../Figs/Fig4_umap.png", height = 6/5, width = 6)

pcp <- ggplot(umaps, aes(x = value)) +
  geom_histogram(aes(fill = Barcode), bins = 100) +
  geom_vline(aes(xintercept = linex, color = Barcode), linetype = "dashed", size = 0.5) +
  facet_wrap(variable ~ alpha, nrow = 2, scales = "free_y") +
  theme_bw() +
  theme(legend.title = element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 6),
        legend.position = "none") +
  scale_color_manual(values = cw_colors[c(2:4,1)], labels = paste("Barcode", 1:4), guide = "none") +
  scale_fill_manual(values = cw_colors[c(2:4,1)], labels = paste("Barcode", 1:4)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 2)))

ggsave(plot = pcp, "../Figs/Fig4_pcs.png", height = 12/5, width = 6)
