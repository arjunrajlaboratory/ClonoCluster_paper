library(magrittr)
library(data.table)
library(ClonoCluster)
library(ggplot2)
library(scales)

# make data table of iters
dt <- data.table::CJ(lin = c(0,1),
                alpha = seq(0, 1, by = 0.1),
                beta = 0.1,
                nn = c(0, 5, 10, 15, 20))

dt[, linf := ifelse(lin == 1, "within barcode", "between barcodes")]

  # compute jaccard index from nearest neighbors
dt[, ji := nn / (40 - nn)]

dt[, jif := paste(scales::percent(nn / 20), "transcriptome\nneighbors")]

dt[, jif := factor(jif, ordered = TRUE, levels = dt[, jif %>% unique])]

dt[, betaf := paste("Beta:", beta)]

dt[, betaf := factor(betaf, ordered = TRUE, levels = dt[, betaf %>% unique])]

# the model
dt[, edge_value := ((alpha ^ beta) * ((lin) - as.numeric(ji)) + as.numeric(ji))]

# plot
p <- ggplot(dt, aes(x = alpha, y = edge_value)) +
  geom_point(aes(col = linf)) +
  geom_line(aes(col = linf)) +
  facet_wrap(~jif, scales = "free", ncol = 1) +
  ylim(0,1) +
  theme_bw() +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(legend.position = "bottom", legend.title = element_blank(),
  legend.text = element_text(size = 6),
  text = element_text(color = "black", family = "Helvetica"),
  axis.text = element_text(color = "black", family = "Helvetica", size = 6),
  axis.title = element_text(size = 6),
  strip.text = element_text(face = "bold", size = 6, family = "Helvetica"),
  strip.background = element_blank()) +
  xlab("\u03B1 value") +
  ylab("Edge Weight (W)") +
  guides(color = guide_legend(nrow = 2))

ggsave(plot = p, paste(fdir, "simulations/Model_iters.pdf", sep = ""), height = 1.5*5, width = 1.5)
