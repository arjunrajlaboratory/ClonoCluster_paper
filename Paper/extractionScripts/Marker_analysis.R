library(magrittr)
library(data.table)
library(ClonoCluster)
library(ggplot2)
library(WebGestaltR)

set.seed(42)

source("Paper/extractionScripts/Constants.R")

# set cores for reorg marker analysis
ncores <- 4

lapply(sn_v, function(sn){

  set.seed(42)

  # get and transpose gene counts
  gt <- data.table::fread(paste(dgdir, sn, "_genes.txt", sep = ""))

  rn <- gt[, rn]

  gt <- gt[, .SD, .SDcols = 2:ncol(gt)] %>% as.matrix %>% t

  colnames(gt) <- rn

  if(sn == "CJ") rownames(gt) %<>% stringr::str_replace("-[0-9]$", "")

  if(sn == "NJ") rownames(gt) <- rownames(gt) %>%
    stringr::str_replace("^S[0-9]_", "") %>%
    stringr::str_replace("-[0-9]$", "")

  # barcodes
  bt <- paste(dbdir, sn, "_barcodes.tsv", sep = "") %>% data.table::fread()

  # make sure only barcoded cells included
  bt <- bt[cellID %chin% rownames(gt)]

  gt <- gt[bt[, cellID], ]

  bt %>% data.table::setnames("cellID", "rn")

  al <- mt[sample_id == sn, c(alevel2, alevel)]

  # get hybrid cluster assignments
  dl <- data.table::fread(
    paste(pdir, sn, "_cluster_assignments.txt", sep = "")
  )

  dl <- dl[alpha %in% c(0, al, 1)]

  # wrapping ROC so it works with apply
  fast_comp <- function(v, x, y){

    a <- ROCR_wrap(x = v[x], y = v[y])

  }

  # get transcriptome markers
  tr <- dl[alpha == 0]

  pb <- utils::txtProgressBar(min = 0,
    max = length(tr[, Group %>% unique %>% length]),
    style = 3)

  # loop to find transcriptome markers
  tr_auc <- lapply(tr[, Group %>% unique], function(go){

    utils::setTxtProgressBar(pb, which(go == tr[, Group %>% unique]))

    ing <- tr[Group == go, rn]

    if (length(ing) < 10) return(NULL)

    outg <- tr[Group != go, rn]

    auc <- apply(gt, 2, function(x) fast_comp(x, ing, outg))

    do <- data.table::as.data.table(auc, keep.rownames = TRUE)

    do[, tr := go]

    return(do)

  }) %>% data.table::rbindlist()

  close(pb)

  # loop to test for reorganization markers and hybrid cluster markers
  spl <- parallel::mclapply(dl[, alpha %>% unique], function(a){

    print("alpha")
    print(a)

    # get hybrid clusters at alpha
    lin <- dl[alpha == a]

    print("getting lin clusters")

    pb <- utils::txtProgressBar(min = 0,
      max = length(lin[, Group %>% unique %>% sort %>% length]),
      style = 3)

    # loop over hybrid clusters at this level and get markers for clusters
    lin_auc <- lapply(lin[, Group %>% unique %>% sort], function(go){

      utils::setTxtProgressBar(pb, which(go == lin[, Group %>% unique %>% sort]))

      print(go)

      # get cells in the group
      ing <- lin[Group == go, rn]

      if (length(ing) < 10) return(NULL)

      # get all other cells
      outg <- lin[Group != go, rn]

      auc <- apply(gt, 2, function(x) fast_comp(x, ing, outg))

      do <- data.table::as.data.table(auc, keep.rownames = TRUE)

      do[, lin := go]

      return(do)

    }) %>% data.table::rbindlist()

    close(pb)

    lin_auc[, alpha := a]

    pb <- utils::txtProgressBar(min = 0,
      max = length(lin[, Group %>% unique %>% length]),
      style = 3)

    # loop over hybrid clusters at this level to get reorg markers
    splits <- lapply(lin[, Group %>% unique], function(go){

      utils::setTxtProgressBar(pb, which(go == lin[, Group %>% unique]))

      print("lineage")
      print(go)

      # get hybrid cluster members
      d <- lin[Group == go]

      # subset transcriptome cluster assignments on the hybrid to get contributing cells
      tri <- tr[rn %chin% d[, rn]]

      # loop over contributing transcriptome clusters
      inner <- lapply(tri[, Group %>% unique %>% sort], function(g){

        print(g)

        # transcriptome cluster cells that go to hybrid cluster of interest
        ing <- tri[Group == g, rn]

        if (length(ing) < 10){

          print("too small")

          return(NULL)
        }

        # all other non-contributing cells in the transcriptome cluster
        outg <- tr[Group == g & !rn %chin% ing, rn]

        # return just straight transcriptome cluster markers for alpha == 0
        # when entire cluster is eaten by another cluster
        if (length(outg) == 0){

          d <- tr_auc[tr == g]

          d[, lin := go]

          return(d)

        }

        # get AUC
        auc <- apply(gt, 2, function(x) fast_comp(x, ing, outg))

        d <- data.table::as.data.table(auc, keep.rownames = TRUE)

        d[, tr := g]

        d[, lin := go]

        return(d)

      }) %>% data.table::rbindlist()

      if (nrow(inner) == 0) return(NULL)

      return(inner)

    }) %>% data.table::rbindlist()

    close(pb)

    splits[, alpha := a]

    splits[, beta := 0.1]

    return(list(splits, lin_auc))

  }, mc.cores = ncores)

  split_l <- lapply(spl %>% seq_along, function(n){

    return(spl[n][[1]][[1]])

  }) %>% data.table::rbindlist(use.names = TRUE)

  lin_l <- lapply(spl %>% seq_along, function(n){

    return(spl[n][[1]][[2]])

  }) %>% data.table::rbindlist(use.names = TRUE)

  # table of reorg marker tests
  split_l %>% data.table::fwrite(paste(pdir, sn, "_auc_barcodes.txt", sep = ""), sep = "\t")

  # table of hybrid cluster markers
  lin_l %>% data.table::fwrite(paste(pdir, sn, "_auc_linclust.txt", sep = ""), sep = "\t")

})

# analysis of reorganization marker overrepresentation
fl <- list.files(pdir, pattern = "auc_barcodes.txt", full.names = TRUE)

# read reorg marker tables for each sample
gst <- lapply(sn_v, function(i){

  f <- fl[fl %like% i]

  d <- data.table::fread(f)

  d[, sample_id := i]

  d[, lab := paste(sample_id, alpha, sep = "_")]

  return(d)

}) %>% data.table::rbindlist()

# GO search space
GOs <- c("geneontology_Cellular_Component_noRedundant",
        "geneontology_Molecular_Function_noRedundant")

# get tempdir
dn <- tempdir()

# loop over samples
fl_ora <- lapply(gst[, sample_id %>% unique], function(si){

  g <- gst[sample_id == si]

  # only evaluate markers with AUC > 0.8
  gl <- g[auc > 80, rn %>% unique, by = c("alpha", "lin", "tr")]

  # split into hybrid cluster contributing transcriptome cluster pairs
  gl %<>% split(by = c("alpha", "lin", "tr"))

  # loop over pairs
  dl <- lapply(gl %>% seq_along, function(n){

    d <- gl[[n]]

    org <- "hsapiens"

    if (si %like% "Klein") org <- "mmusculus"

    gs <- try(WebGestaltR::WebGestaltR(enrichMethod = "ORA", organism = org, enrichDatabase = GOs,
            enrichDatabaseFile = NULL, enrichDatabaseType = NULL, enrichDatabaseDescriptionFile = NULL,
            interestGeneFile = NULL, interestGene = d[, V1], interestGeneType = "genesymbol",
            collapseMethod = "mean", referenceGeneFile = NULL, referenceGene = NULL,
            referenceGeneType = "genesymbol", referenceSet = "genome_protein-coding", minNum = 10,
            maxNum = 500, sigMethod = "fdr", fdrMethod = "BH", fdrThr = 1,
            topThr = 10, reportNum = 20, perNum = 1000, gseaP = 1, isOutput = TRUE,
            outputDirectory = dn, projectName = paste(si, "_ORA", sep = ""), dagColor = "continuous",
            saveRawGseaResult = FALSE, gseaPlotFormat = c("png", "svg"),
            setCoverNum = 10, networkConstructionMethod = NULL, neighborNum = 10,
            highlightType = "Seeds", highlightSeedNum = 10, nThreads = 1,
            cache = NULL, hostName = "http://www.webgestalt.org/"))

    # handle error of empty return if no significantly diff gene sets
    if(class(gs)[1] == "try-error"){

      warning(paste("Try-error at alpha", d[, alpha %>% unique], "iter", d[, lin %>% unique]))

      return(NULL)

    }

    if (is.null(gs) || nrow(gs) == 0) return(NULL)

    # curate output columns
    gs %<>% as.data.table %>% .[, .SD, .SDcols = c(1:2, 7:9)]

    gs[, alpha := d[, alpha %>% unique]]

    gs[, lin := d[, lin %>% unique]]

    gs[, tr := d[, tr %>% unique]]

    return(gs)

  }) %>% data.table::rbindlist()

  dl[, sample_id := si]

  dl %>% data.table::fwrite(paste(pdir, si, "_rearrangement_ORA.txt", sep = ""), sep = "\t")

  return(paste(pdir, si, "_rearrangement_ORA.txt"))

}) %>% unlist
