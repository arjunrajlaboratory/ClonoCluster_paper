library(data.table)
library(magrittr)
library(BarCluster)
library(ggplot2)

set.seed(42)

# scale to size we want to make this (10 to a side)
base_size <- 10

dotsn <- base_size ^ 2

# simulate our data
b <- rnorm(n = 1000, mean = -100, sd = 3)

a <- rnorm(n = 1000, mean = 100, sd = 5)

c <- rnorm(n = 1000, mean = 0, sd = 4)

bt <- lapply(c("A", "B", "C"), function(x) replicate(floor(dotsn/3), x)) %>% unlist

bt <- data.table::data.table(rn = 1:length(bt),  Barcode = bt)

ii <- ceiling(nrow(bt) * 5/6)

ii2 <- floor(nrow(bt) * 9/12)

# do some mixing of data
bt[rn %in% ii:nrow(bt), Barcode := "D"]

bt[rn %in% ii2:(ii - 1), Barcode := "B"]

dt <- data.table::CJ(1:base_size, 1:base_size)

dt[, rn := 1:nrow(dt)]

dt <- merge(dt, bt, by = "rn")

dt %>% setnames(c("V1", "V2"), c("x", "y"))

c1 <- dt[y > 0 & y < quantile(y, 0.33), rn]

c2 <- dt[y >= quantile(y, 0.33) & y <= quantile(y, 0.67), rn]

c3 <- dt[y > quantile(y, 0.67), rn]

pcs1 <- cbind(a,b,c,a,b,c) %>% .[1:length(c1), ]

pcs2 <- cbind(b,c,a,b,c,a) %>% .[1:length(c2), ]

pcs3 <- cbind(c,b,a,c,b,a) %>% .[1:length(c3), ]

# 3 pcs
irl <- rbind(pcs1,pcs2,pcs3)

rownames(irl) <- c(c1,c2,c3)

# reorder it correctly
irl <- irl[dt[, rn], ]

# get "transcriptome" graph
neighbor.graphs <- Seurat::FindNeighbors(object = irl, k.param = 20,
          compute.SNN = TRUE, prune.SNN = 1/15,
          nn.method = "rann", annoy.metric = "euclidean",
          nn.eps = 0, verbose = TRUE, force.recalc = FALSE)

m <- neighbor.graphs$snn

# alphas to look at
alphas <- c(0, 0.5, 1)

# barcode matrix
nm <- build_barcode_matrix_fast(bt)

# return hybrid graph and cluster assignments over alphas
ml <- lapply(alphas, function(a){

  mm2 <- barcluster_model(alpha = a, beta = 1, m = m, nm = nm)

  d <- barcluster(irl, bt, alpha = a, beta = 1, res = 2)

  return(list(mm2, d))

})

names(ml) <- alphas

# loop over alpha levels
pl <- lapply(ml %>% seq_along, function(n){

  m <- ml[[n]]

  al <- names(ml)[n]

  # this will be table of edges to draw and weights
  wt <- m[[1]] %>% as.matrix %>% reshape2::melt() %>% as.data.table()

  # this is the point positions in space of nodes, even 10 by 10 grid
  dt1 <- dt %>% data.table::copy()

  startt <- dt1[, .(rn, x, y)] %>% data.table::copy() %>%
    setnames("rn", "start")

  endt <- dt1[, .(rn, x, y)] %>% data.table::copy() %>%
    setnames(c("rn", "x", "y"), c("end", "xend", "yend"))

  wt %>% setnames(c("Var1", "Var2"), c("start", "end"))

  wt <- merge(startt, wt, by = "start")

  wt <- merge(endt, wt, by = "end")

  d <- m[[2]] %>% data.table::copy() %>% .[, .(rn, Group)] %>%
     .[, rn := as.numeric(rn)]

  dt1 <- merge(dt1, d, by = "rn")

  # text base size
  bs <- 3

  # plot
  p <- ggplot(dt1, aes(x = x, y = y)) +
    geom_segment(data = wt[value > 0], aes(xend = xend, yend = yend, alpha = value), size = 0.1) +
    geom_point(size = bs, col = "white") +
    geom_point(size = bs, shape = 1, col = "black") +
    geom_text(aes(col = as.factor(Group), label = Barcode), size = (bs - 1), fontface = "bold") +
    scale_color_manual(values = BarCluster::c25) +
    theme_void() +
    theme(legend.position = "none", plot.title = element_text(size = 8, hjust = 0.5, face = "bold")) +
    ggtitle(paste("\u03B1", "=", al))

  return(p)

})

plot <- cowplot::plot_grid(plotlist = pl, nrow = 1, scale = 0.9)

ggsave(plot = plot, "../Figs/grid_plot.pdf", height = 2, width = 6)
