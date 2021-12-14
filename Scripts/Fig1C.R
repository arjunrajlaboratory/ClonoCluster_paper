library(magrittr)
library(data.table)
library(BarCluster)
library(ggplot2)

set.seed(42)

sn_v <- c(paste("YG", 1:3, sep = ""), "Kleind2ws", "CJ")

tv <- c("High dose BRAFi - 1", "High dose BRAFi - 2", "Low dose BRAFi", "In Vitro Hematopoiesis Day 2", "Cardio-directed iPSC", "Neuro-directed iPSC")

names(tv) <- sn_v

lapply(sn_v, function(sn){

  set.seed(42)

  res <- 1

  if (sn %like% "Klein") res <- 0.6

  gt <- data.table::fread(paste("../Data_genes/", sn, "_genes.txt", sep = ""))

  rn <- gt[, rn]

  gt <- gt[, .SD, .SDcols = 2:ncol(gt)] %>% as.matrix %>% t

  colnames(gt) <- rn

  if(sn == "CJ") rownames(gt) %<>% stringr::str_replace("-[0-9]$", "")

  if(sn == "NJ") rownames(gt) <- rownames(gt) %>%
    stringr::str_replace("^S[0-9]_", "") %>%
    stringr::str_replace("-[0-9]$", "")

  bt <- paste("../Data_barcodes/", sn, "_barcodes.tsv", sep = "") %>% data.table::fread()

  bt <- bt[cellID %chin% rownames(gt)]

  gt <- gt[bt[, cellID], ]

  bt %>% data.table::setnames("cellID", "rn")

  wl <- bt[, .N, by = c("Barcode")] %>% .[order(-N), Barcode %>% unique %>% .[1:15]]

  rnl <- bt[Barcode %chin% wl, rn]

  irl <- irlba_wrap(gt, npc = 100, seed_use = 42)

  nm <- build_barcode_matrix(bt)

  if (!all(sort(rownames(nm)) == sort(rownames(irl)))) stop("barcodes names don't match PCA")

  neighbor.graphs <- Seurat::FindNeighbors(object = irl, k.param = 20,
            compute.SNN = TRUE, prune.SNN = 1/15,
            nn.method = "rann", annoy.metric = "euclidean",
            nn.eps = 0, verbose = TRUE, force.recalc = FALSE)

  m <- neighbor.graphs$snn

  rm("neighbor.graphs")

  beta_val <- 0.1

  dl <- lapply(seq(0, 1, by = 0.05), function(alpha){

    mm2 <- barcluster_model(alpha = alpha, beta = beta_val, m = m, nm = nm)

    ids_b <- RunModularityClustering(SNN = mm2, modularity = 1,
            resolution = res, algorithm = 1, n.start = 1,
            n.iter = 10, random.seed = 42, print.output = TRUE,
            temp.file.location = NULL, edge.file.name = NULL)

    names(ids_b) <- colnames(mm2)

    ids_b <- matrix(ids_b, ncol = 1, dimnames = list(names(ids_b)))

    colnames(ids_b) <- "Group"

    ids_b %<>% as.data.table(keep.rownames = TRUE)

    ids_b[, alpha := alpha]

    ids_b[, resolution := res]

    return(ids_b)

  }) %>% data.table::rbindlist()

  df <- dl %>% data.table::copy()

  df[, sample_id := sn] %>% .[, beta := beta_val]

  df %>% data.table::fwrite(
    paste("../Processed_data/", sn, "_cluster_assignments.txt", sep = ""),
    sep = "\t")

  # fuzzy matching because R is so busted
  alphlist <- dl[, alpha %>% unique %>% sort] %>% .[seq(1, 21, by = 2)]

  p <- Plot_alluvia(dl[alpha %in% alphlist & rn %chin% rnl],
                    bt[Barcode %chin% wl],
                    ylab = "Number of cells",
                    xlab = "Alpha value",
                    title = paste(tv[names(tv) == sn], "clusters"),
                    label_nodes = FALSE,
                    ltype = "text",
                    border_size = 0.75)

  ttheme <- theme(axis.title = element_text(size = 8, face = "bold", color = "black"),
                  axis.text = element_text(size = 8, color = "black"),
                  plot.title = element_text(size = 8, face = "bold", color = "black"),
                  plot.subtitle = element_text(size = 8, color = "black"))

  p[[1]] <- p[[1]] + ttheme

  p[[2]] <- p[[2]] + ttheme

  pl <- cowplot::plot_grid(plotlist = p, ncol = 1, scale = 0.9)

  ggsave(plot = pl, filename = paste("../Figs/", sn, "_alluvia.png", sep = ""), width = 7.5, height = 6)

})
