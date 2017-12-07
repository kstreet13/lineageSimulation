###################################################
### Lineage Inference Simulation Study
### Author: Kelly Street
###################################################
{
  # library(splatter)
  source('code/helper_functions.R')
  # samples <- system('ls data/trapnell14/salmon/', intern=TRUE)
  # hsmmCounts <- sapply(samples, function(samp){
  #   fname <- paste0('data/trapnell14/salmon/',samp,'/quant.sf')
  #   tab <- read.table(fname, header = TRUE)
  #   return(as.numeric(tab$NumReads))
  # })
  # rownames(hsmmCounts) <- as.character(read.table(
  #   paste0('data/trapnell14/salmon/',samples[1],'/quant.sf'), 
  #   header = TRUE)$Name)
  # m <- quantile(hsmmCounts[hsmmCounts > 0], probs = .5)
  # # filter genes based on robust expression in at least 10% of cells
  # gfilt <- apply(hsmmCounts,1,function(g){
  #   sum(g > m) >= .10 * ncol(hsmmCounts)
  # })
  # # filter cells: no more than 70% zero counts
  # # there seems to be a break around 70% separating two groups of cells
  # sfilt <- apply(hsmmCounts[filt,],2,function(s){
  #   mean(s == 0) <= .7
  # })
  # parHSMM <- splatEstimate(round(hsmmCounts[gfilt,sfilt]))
  load('data/parHSMM.RData')
} # load splatter and get parameters from HSMM data


cells.range <- seq(20,120, by=20)
deProb.range <- (1:5)/10
reps <- 10
it.params <- cbind(rep(cells.range,each=reps*length(deProb.range)), rep(deProb.range,each=reps,times=length(cells.range)))
nIts <- nrow(it.params)
smallestOKcluster <- 5

library(SingleCellExperiment)
library(slingshot)
library(destiny)
library(fastICA)
library(mclust)
# library(monocle)
# library(TSCAN)
library(parallel)

# i = 40
simResults_sling52 <- mclapply(1:nIts, function(i){
  print(paste('Iteration',i))
  {
    gpcells <- it.params[i,1]
    deProb <- it.params[i,2]
    splat.it.params <- splatter::setParams(parHSMM, update = list(batchCells = 11*gpcells, group.prob = c(1,2,1,1,1,1,2,1,1)/11,
                                                        path.from = c(0,1,2,2,1,5,5,6,6), path.nonlinearProb = .5,
                                                        de.prob = deProb, path.length = c(1,2,1,1,1,1,2,1,1)*100, 
                                                        de.facLoc = .25, de.facScale = 1, dropout.present = TRUE,
                                                        seed = i))
    sim <- splatter::splatSimulate(splat.it.params, method='paths')
    truth.vec <- sim$Step + 100*(sim$Group %in% c('Path2','Path5')) + 200*(sim$Group %in% c('Path6','Path7')) + 300*(sim$Group %in% c('Path3','Path4','Path8','Path9'))
    sim$Step <- truth.vec
    start.cell <- colnames(sim)[which.min(sim$Step)]
    
    assays(sim)$normcounts <- FQnorm(assays(sim)$counts)
    
    true.lineages <- list(Lineage1 = c('Path1','Path2','Path3'),
                          Lineage2 = c('Path1','Path2','Path4'),
                          Lineage3 = c('Path1','Path5','Path7'),
                          Lineage4 = c('Path1','Path5','Path6','Path8'),
                          Lineage5 = c('Path1','Path5','Path6','Path9'))
    
    truth <- matrix(rep(truth.vec,5), ncol=5)
    for(l in 1:5){
      truth[! sim$Group %in% true.lineages[[l]], l] <- NA
    }
    rownames(truth) <- colnames(sim)
    colnames(truth) <- names(true.lineages)
  }# simulate, normalize data
  out <- NULL
  
  pca <- prcomp(t(log1p(assays(sim)$normcounts)))
  # 4-D ICA
  p <- 4
  ica <- fastICA::fastICA(pca$x, n.comp = p, method = 'C')$S
  rownames(ica) <- rownames(pca$x)
  X <- ica
  X.orig <- X
  mc <- Mclust(X)
  cl <- mc$classification
  minClusSize <- min(table(cl))
  while(minClusSize < smallestOKcluster){
    clTab <- table(cl)
    smallClus <- names(clTab)[clTab < smallestOKcluster]
    X <- X[! cl %in% smallClus, ]
    mc <- Mclust(X)
    cl <- mc$classification
    minClusSize <- min(table(cl))
    if(nrow(X) < .5*nrow(X.orig)){ # prevent dropping too many cells
      X <- X.orig
      mc <- Mclust(X)
      cl <- mc$classification
      minClusSize <- smallestOKcluster + 1
    }
  }
  clus1 <- as.character(cl[which.max(names(cl)==start.cell)])
  pstIC <- tryCatch({
    slingshot(X, cl, start.clus = clus1)
  }, error = function(e){ e })
  if('error' %in% class(pstIC)){
    out <- cbind(out,
                 c(nlins = NA, rep(NA,ncol(truth))))
  }else{
    pstIC <- pseudotime(pstIC)
    out <- cbind(out,
                 c(nlins = ncol(pstIC), agreement.pair(truth, pstIC)))
  }
  colnames(out)[ncol(out)] <- paste0('slingIC',p,'_MCbest')
  
  # 8-D Diffusion Map
  dm <- destiny::DiffusionMap(t(log1p(assays(sim)$normcounts)))@eigenvectors
  rownames(dm) <- colnames(sim)
  p <- 8
  X <- dm[,1:p]
  X.orig <- X
  mc <- Mclust(X)
  cl <- mc$classification
  minClusSize <- min(table(cl))
  while(minClusSize < smallestOKcluster){
    clTab <- table(cl)
    smallClus <- names(clTab)[clTab < smallestOKcluster]
    X <- X[! cl %in% smallClus, ]
    mc <- Mclust(X)
    cl <- mc$classification
    minClusSize <- min(table(cl))
    if(nrow(X) < .5*nrow(X.orig)){ # prevent dropping too many cells
      X <- X.orig
      mc <- Mclust(X)
      cl <- mc$classification
      minClusSize <- smallestOKcluster + 1
    }
  }
  clus1 <- as.character(cl[which.max(names(cl)==start.cell)])
  pstDM <- tryCatch({
    slingshot(X, cl, start.clus = clus1)
  }, error = function(e){ e })
  if('error' %in% class(pstDM)){
    out <- cbind(out,
                 c(nlins = NA, rep(NA,ncol(truth))))
  }else{
    pstDM <- pseudotime(pstDM)
    out <- cbind(out,
                 c(nlins = ncol(pstDM), agreement.pair(truth, pstDM)))
  }
  colnames(out)[ncol(out)] <- paste0('slingDM',p,'_MCbest')
  
  return(out)
}, mc.cores = 16)

save(simResults_sling52, file='temp_sling52.RData')

# format output
nlins <- t(sapply(simResults_sling52, function(res){
  res[1,]
}))
taus <- lapply(1:ncol(simResults_sling52[[1]]),function(i){
  t(sapply(simResults_sling52, function(res){
    res[-1, i]
  }))
})
names(taus) <- colnames(simResults_sling52[[1]])

save(nlins,taus, file='results_sling52.RData')

