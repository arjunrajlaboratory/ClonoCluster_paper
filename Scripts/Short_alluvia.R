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

lapply(sn_v, function(sn){

  al <- c(mt[sample_id == sn, alevel2], mt[sample_id == sn, alevel])

  dl <- list.files("../Processed_data", pattern = "assignments", full.names = TRUE) %>%
  .[. %like% sn] %>% data.table::fread()

  bcs <- list.files("../Data_barcodes", full.names = TRUE)

  bt <- bcs[bcs %like% sn] %>% data.table::fread()

  bt %>% setnames("cellID", "rn")

  wl <- bt[, .N, by = "Barcode"] %>% .[order(-N)] %>%
      .[, Barcode %>% unique %>% .[1:15]]

  ttheme <- theme(axis.title = element_text(size = 8, face = "bold", color = "black"),
                  axis.text = element_text(size = 8, color = "black"),
                  axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1),
                  plot.title = element_text(size = 8, face = "bold", color = "black", hjust = 0.5),
                  plot.subtitle = element_text(size = 8, color = "black"))

  rnl <- bt[Barcode %chin% wl, rn]

  p <- BarCluster::Plot_alluvia(dl[alpha %in% c(0, al) & rn %chin% rnl],
                    bt[Barcode %chin% wl],
                    ylab = "Number of cells",
                    xlab = NULL,
                    title = paste(tv[names(tv) == sn], "clusters"),
                    label_nodes = FALSE,
                    ltype = "text",
                    border_size = 0.75)

  p[[1]] <- p[[1]] + ttheme +
    scale_x_discrete(labels = c("Transcriptome", "Low alpha", "High alpha", "Barcodes"))

  p[[2]] <- p[[2]] + ttheme +
    scale_x_discrete(labels = c("Transcriptome", "Low alpha", "High alpha", "Barcodes"))

  plot <- cowplot::plot_grid(plotlist = p, nrow = 1, scale = 0.9)

  ggsave(plot = plot, paste("../Figs/",sn,"_short_alluvia.png", sep = ""), height = 2.75, width = 5.5)

})
