library(magrittr)
library(data.table)
library(BarCluster)
library(ggplot2)

sn_v <- c(paste("YG", 1:3, sep = ""), "Kleind2ws", "CJ")

tv <- c("High dose BRAFi - 1", "High dose BRAFi - 2", "Low dose BRAFi", "In Vitro Hematopoiesis Day 2", "Cardio-directed iPSC", "Neuro-directed iPSC")

names(tv) <- sn_v

# cutoffs determined by max alpha before cluster number increases
mt <- data.table::data.table(sample_id = sn_v,
  alevel = c(0.55, 0.6, 0.6, 0.75, 0.55))

mt[, alevel2 := (alevel / 2) %>% round(digits = 1)]

lapply(sn_v, function(sn){

  print(sn)

  fl <- list.files("../Processed_data/", full.names = TRUE, pattern = "cluster_assignment")

  dt <- fl[fl %like% sn] %>% data.table::fread()

  # cluster number
  n <- dt[, Group %>% unique %>% length, by = c("alpha")]

  ttheme <- theme(axis.title = element_text(size = 8, face = "bold", color = "black"),
                  axis.text = element_text(size = 8, color = "black"),
                  plot.title = element_text(size = 8, face = "bold", color = "black", hjust = 0.5),
                  plot.subtitle = element_text(size = 8, color = "black"))

  p1 <- ggplot(n, aes(x = alpha, y = V1)) +
      geom_line(linetype = "solid", color = "grey", size = 0.75) +
      geom_point(color = "dodgerblue", size = 1, alpha = 0.6) +
      theme_bw() +
      ttheme +
      xlab("Alpha value") +
      ylab("# of clusters") +
      scale_y_log10() +
      ggtitle(tv[names(tv) == sn]) +
      theme(axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1)) +
      scale_x_continuous(breaks = seq(0,1, by = 0.2), limits = c(-0.15,1.1))

  p1 <- p1 +
  annotate("text", x = mt[sample_id == sn, alevel2],
  y = n[alpha == mt[sample_id == sn, alevel2], V1 + 20] ,label = "Low\nalpha", size = 2)

  p1 <- p1 +
  annotate("text", x = mt[sample_id == sn, alevel],
  y = n[alpha == mt[sample_id == sn, alevel], V1 + 20] ,label = "High\nalpha", size = 2)

  p1 <- p1 +
  annotate("text", x = 0,
  y = n[alpha == 0, V1 + 20] ,label = "Transcriptome\nclusters", size = 2)

  p1 <- p1 +
  annotate("text", x = 1,
  y = n[alpha == 1, V1 * 1.2] ,label = "Barcodes", size = 2)

  ggsave(plot = p1, paste0("../Figs/", sn, "_cluster_size_plot.png"), width = 2.5, height = 2.5, limitsize = FALSE)

})
