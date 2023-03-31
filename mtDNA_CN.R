
mtDNA_CN <- function(res_file, region){
  
  chrA <- read.table(paste0(res_file,"_",region), header = T, sep = "\t")
  chrM <- read.table(paste0(res_file,"_chrM"), header = T, sep = "\t")
  chrX <- read.table(paste0(res_file,"_chrX"), header = T, sep = "\t")
  chrY <- read.table(paste0(res_file,"_chrY"), header = T, sep = "\t")

  # 6 bins: 30% - 60%
  # 5%: K=6, 
  g <- seq(33,58,5)
  len.g <- length(g)
  mi <- ai <- xi <- yi <- c()
  MB <- AB <- XB <- YB <- c()
  binM <- binA <- binX <- binY <- list()
  
  if(all(is.na(chrA$RC))){
    stop("Coverage calculation in autosome wasn't completed.")
  }else if(any(is.na(chrA$RC))){
    chrA <- chrA[!(is.na(chrA$RC)),]
  }
  
  if(any(is.na(chrM$RC))){
    stop("Coverage calculation in chrM wasn't completed.")
  }
  
  N <- sum(chrA$RC) 
  p <- 1/(nrow(chrA))

  for(ll in 1:len.g){
    binM[[ll]] <- chrM[chrM$GC >= g[ll]-3 & chrM$GC < g[ll]+2,]
    mi <- c(mi,nrow(binM[[ll]]))
    MB <- c(MB,sum(binM[[ll]]$RC))

    binX[[ll]] <- chrX[chrX$GC >= g[ll]-3 & chrX$GC < g[ll]+2,]
    xi <- c(xi,nrow(binX[[ll]]))
    XB <- c(XB,sum(binX[[ll]]$RC))

    binY[[ll]] <- chrY[chrY$GC >= g[ll]-3 & chrY$GC < g[ll]+2,]
    yi <- c(yi,nrow(binY[[ll]]))
    YB <- c(YB,sum(binY[[ll]]$RC))
    
    binA[[ll]] <- chrA[(chrA$GC >= g[ll]-3) & (chrA$GC < g[ll]+2), ]
    #binA[[ll]] <- binA[[ll]][!(is.na(binA[[ll]][,1])),]
    ai <- c(ai,nrow(binA[[ll]]))
    AB <- c(AB,sum(binA[[ll]]$RC))
  }
  m <- sum(mi); a <- sum(ai); x <- sum(xi); y <- sum(yi)
  M <- sum(MB); X <- sum(XB); Y <- sum(YB)
  t.m <- MB+AB; t.x <- XB+AB; t.y <- YB+AB
  
  b0.m <- rep(0.5,len.g)
  mu.m <- 1000 
  mu0.m <- 0
  while(abs(mu.m-mu0.m)>0.0001){
    mu0.m <- mu.m
    mu.m <- M/(N*p*sum(mi*b0.m))
    b.m <- t.m/(N*p*(mu.m*mi+2*ai))
    b0.m <- b.m
  }

  b0.x <- rep(0.5,len.g)
  mu.x <- 1000 
  mu0.x <- 0
  while(abs(mu.x-mu0.x)>0.0001){
    mu0.x <- mu.x
    mu.x <- X/(N*p*sum(xi*b0.x))
    b.x <- t.x/(N*p*(mu.x*xi+2*ai))
    b0.x <- b.x
  }

  b0.y <- rep(0.5,len.g)
  mu.y <- 1000 
  mu0.y <- 0
  while(abs(mu.y-mu0.y)>0.0001){
    mu0.y <- mu.y
    mu.y <- Y/(N*p*sum(yi*b0.y))
    b.y <- t.y/(N*p*(mu.y*yi+2*ai))
    b0.y <- b.y
  }
  var.m <- 1/(N*p*(sum(mi*b.m)/mu.m-sum((mi^2*b.m)/(mu.m*mi+2*ai))))
  res <- as.data.frame(matrix(c(region,round(c(N*p, sum(mi*b.m), sum(ai)/sum(22*mi), mu.m, var.m, mu.x, mu.y),2)),nrow = 1))
  colnames(res) <- c("region","Np","m_beta","k","mt","var.mt","chrX","chrY")
  write.table(res,paste0(res_file,"_",region,".mitoCN.txt"),row.names = F, quote = F, sep = "\t")
  
}

args <- commandArgs(trailingOnly = TRUE)
res_file <- args[1]
region <- args[2]
mtDNA_CN(res_file,region)
