library(magrittr)
library(data.table)
library(BarCluster)
library(ggplot2)

set.seed(42)

source("Constants.R")

lapply(sn_v, function(sn){

  # get alpha levels
  al <- c(mt[sample_id == sn, alevel2], mt[sample_id == sn, alevel])

  # get assignments
  dl <- list.files("../Processed_data", pattern = "assignments", full.names = TRUE) %>%
  .[. %like% sn] %>% data.table::fread()

  # get barcodes
  bcs <- list.files("../Data_barcodes", full.names = TRUE)

  bt <- bcs[bcs %like% sn] %>% data.table::fread()

  bt %>% setnames("cellID", "rn")
  #
  # wl <- bt[, .N, by = "Barcode"] %>% .[order(-N)] %>%
  #     .[, Barcode %>% unique %>% .[1:15]]

  ttheme <- theme(axis.title = element_text(size = 8, face = "bold", color = "black"),
                  axis.text = element_text(size = 8, color = "black"),
                  axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1),
                  plot.title = element_text(size = 8, face = "bold", color = "black", hjust = 0.5),
                  plot.subtitle = element_blank())

  #rnl <- bt[Barcode %chin% wl, rn]

  # plot
  p <- BarCluster::Plot_alluvia(dl[alpha %in% c(0, al)],
                    bt,
                    ylab = "Number of cells",
                    xlab = NULL,
                    title = paste(tv[names(tv) == sn], "clusters"),
                    label_nodes = FALSE,
                    ltype = "text",
                    border_size = 0.75)

  # change x axis labels
  lvec <- c(paste("Transcriptome\n", "@", " = 0", sep = ""),
  paste("Low @\n", "@", " = ", al[1], sep = ""),
  paste("High @\n", "@", " = ", al[2], sep = ""),
  paste("Barcodes\n", "@", " = 1", sep = ""))

  p[[1]] <- p[[1]] + ttheme +
    scale_x_discrete(labels = lvec)

  p[[2]] <- p[[2]] + ttheme +
    scale_x_discrete(labels = lvec)

  plot <- cowplot::plot_grid(plotlist = p, nrow = 1, scale = 0.9)

  ggsave(plot = plot, paste("../Figs/Short_alluvia/",sn,"_short_alluvia.pdf", sep = ""), height = 2.75, width = 5.5)

})
