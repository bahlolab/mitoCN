# Define mtDNA_CN function
mtDNA_CN <- function(res_file, region){
  
  # Read in depth files for autosomes, mitochondrial DNA, chrX, and chrY
  chrA <- read.table(paste0(res_file,"_at"), header = TRUE, sep = "\t")
  chrM <- read.table(paste0(res_file,"_mt"), header = TRUE, sep = "\t")
  chrX <- read.table(paste0(res_file,"_x"), header = TRUE, sep = "\t")
  chrY <- read.table(paste0(res_file,"_y"), header = TRUE, sep = "\t")

  # Define GC content bin boundaries: 6 bins from 33% to 58%, spaced by 5%
  g <- seq(33, 58, 5)
  len.g <- length(g)

  # Initialize counters
  mi <- ai <- xi <- yi <- c()        # bin counts
  MB <- AB <- XB <- YB <- c()        # read depth sums
  binM <- binA <- binX <- binY <- list()

  # Sanity check: autosome coverage must be complete
  if(all(is.na(chrA$RC))){
    stop("Coverage calculation in autosome wasn't completed.")
  } else if(any(is.na(chrA$RC))){
    chrA <- chrA[!is.na(chrA$RC),]
  }

  # Sanity check: mitochondrial coverage must be complete
  if(any(is.na(chrM$RC))){
    stop("Coverage calculation in chrM wasn't completed.")
  }

  # Total autosomal read count and mean probability (inverse of bin count)
  N <- sum(chrA$RC)
  p <- 1 / nrow(chrA)

  # Loop through GC bins and accumulate statistics
  for(ll in 1:len.g){
    # Bin chrM, chrX, chrY, and autosomes based on GC content
    binM[[ll]] <- chrM[chrM$GC >= g[ll] - 3 & chrM$GC < g[ll] + 2,]
    mi <- c(mi, nrow(binM[[ll]]))
    MB <- c(MB, sum(binM[[ll]]$RC))

    binX[[ll]] <- chrX[chrX$GC >= g[ll] - 3 & chrX$GC < g[ll] + 2,]
    xi <- c(xi, nrow(binX[[ll]]))
    XB <- c(XB, sum(binX[[ll]]$RC))

    binY[[ll]] <- chrY[chrY$GC >= g[ll] - 3 & chrY$GC < g[ll] + 2,]
    yi <- c(yi, nrow(binY[[ll]]))
    YB <- c(YB, sum(binY[[ll]]$RC))

    binA[[ll]] <- chrA[chrA$GC >= g[ll] - 3 & chrA$GC < g[ll] + 2,]
    ai <- c(ai, nrow(binA[[ll]]))
    AB <- c(AB, sum(binA[[ll]]$RC))
  }

  # Total bin sizes and depths
  m <- sum(mi); a <- sum(ai); x <- sum(xi); y <- sum(yi)
  M <- sum(MB); X <- sum(XB); Y <- sum(YB)

  # Combined read depth for each (mt + autosome, X + autosome, Y + autosome)
  t.m <- MB + AB
  t.x <- XB + AB
  t.y <- YB + AB

  # --- Estimate mtDNA-CN (mu.m) using iterative method (EM-like) ---
  b0.m <- rep(0.5, len.g)
  mu.m <- 1000
  mu0.m <- 0

  while(abs(mu.m - mu0.m) > 0.0001){
    mu0.m <- mu.m
    mu.m <- M / (N * p * sum(mi * b0.m))
    b.m <- t.m / (N * p * (mu.m * mi + 2 * ai))
    b0.m <- b.m
  }

  # --- Estimate chrX coverage (mu.x) ---
  b0.x <- rep(0.5, len.g)
  mu.x <- 1000
  mu0.x <- 0

  while(abs(mu.x - mu0.x) > 0.0001){
    mu0.x <- mu.x
    mu.x <- X / (N * p * sum(xi * b0.x))
    b.x <- t.x / (N * p * (mu.x * xi + 2 * ai))
    b0.x <- b.x
  }

  # --- Estimate chrY coverage (mu.y) ---
  b0.y <- rep(0.5, len.g)
  mu.y <- 1000
  mu0.y <- 0

  while(abs(mu.y - mu0.y) > 0.0001){
    mu0.y <- mu.y
    mu.y <- Y / (N * p * sum(yi * b0.y))
    b.y <- t.y / (N * p * (mu.y * yi + 2 * ai))
    b0.y <- b.y
  }

  # --- Estimate variance of mtDNA-CN ---
  var.m <- 1 / (N * p * (sum(mi * b.m) / mu.m - sum((mi^2 * b.m) / (mu.m * mi + 2 * ai))))

  # --- Prepare output ---
  res <- as.data.frame(matrix(
    c(region, round(c(N * p, sum(mi * b.m), sum(ai) / sum(22 * mi), mu.m, var.m, mu.x, mu.y), 2)),
    nrow = 1
  ))
  colnames(res) <- c("region", "Np", "m_beta", "k", "mt", "var.mt", "chrX", "chrY")

  # Write output table to a text file
  write.table(
    res,
    paste0(res_file, "_", region, ".mitoCN.txt"),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
res_file <- args[1]  # Sample name (prefix)
region <- args[2]    # Region name (e.g., "k500")

# Run mtDNA_CN function
mtDNA_CN(res_file, region)
