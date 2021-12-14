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
lapply(sn_v, function(sn){

  al <- c(mt[sample_id == sn, alevel2], mt[sample_id == sn, alevel])

  glo <- c("COL6A1", "Col6a1", "ITGB1", "Itgb1")

  gt <- data.table::fread(paste("../Data_genes/", sn, "_genes.txt", sep = ""))

  lapply(glo, function(gl){

    g <- gl[gl %chin% gt[, rn]]

    if (length(g) == 0) return(NULL)

    gt <- gt[rn == g, .SD, .SDcols = 2:ncol(gt)] %>% t

    colnames(gt) <- g

    if(sn == "CJ") rownames(gt) %<>% stringr::str_replace("-[0-9]$", "")

    if(sn == "NJ") rownames(gt) <- rownames(gt) %>%
      stringr::str_replace("^S[0-9]_", "") %>%
      stringr::str_replace("-[0-9]$", "")

    dl <- data.table::fread(paste("../Processed_data/", sn, "_cluster_assignments.txt", sep = ""))

    ttheme <- theme(axis.title = element_text(size = 8, face = "bold", color = "black"),
                    axis.text = element_text(size = 8, color = "black"),
                    axis.text.x = element_text(size = 8, color = "black", angle = 35, hjust = 1),
                    plot.title = element_text(size = 8, face = "bold", color = "black", hjust = 0.5),
                    plot.subtitle = element_text(size = 8, color = "black"))

    dl <- dl[alpha %in% c(0, al)]

    p2 <- Plot_alluvia_counts(dl,
                              counts = gt,
                              title = paste(tv[names(tv) == sn], ": ", g, sep = ""),
                              ylab = "Number of cells",
                              xlab = NULL,
                              border_size = 0.25,
                              label_nodes = FALSE
                              )

    p2 <- p2 + ttheme +
      scale_x_discrete(labels = c("Transcriptome", "Low alpha", "High alpha", "Barcodes")) +
      scale_color_gradient(low = "gray100", high = "darkblue", breaks = c(0, 5, 10, 15), limits = c(-3,18)) +
      theme(legend.key.size = unit(2.5, "mm"))

    ggsave(plot = p2, paste("../Figs/", sn, "_", g, "_UMI_alluvia.png", sep = ""), height = 3, width = 3)

  })

})
