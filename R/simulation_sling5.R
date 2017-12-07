###################################################
### Lineage Inference Simulation Study
### Slingshot - other methods' best cases
### Five Lineages
###################################################

# setup
NCORES <- 1
source('R/helper_functions.R')
load('data/parHSMM.RData')
library(SingleCellExperiment)
library(parallel)
library(slingshot)
library(destiny)
library(fastICA)
library(mclust)

# five-lineage case setup
cells.range <- seq(20,120, by=20)
deProb.range <- (1:5)/10
reps <- 10
it.params <- cbind(rep(cells.range,each=reps*length(deProb.range)), rep(deProb.range,each=reps,times=length(cells.range)))
nIts <- nrow(it.params)
smallestOKcluster <- 5

simResults_sling5 <- mclapply(1:nIts, function(i){
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
  
  # 4-D PCA - to ~match TSCAN (hybrid method is more "matched")
  pca <- prcomp(t(log1p(assays(sim)$normcounts)))
  p <- 4
  X <- pca$x[,1:p]
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
  pstPC <- tryCatch({
      slingshot(X, cl, start.clus = clus1)
  }, error = function(e){ e })
  if('error' %in% class(pstPC)){
      out <- cbind(out,
                   c(nlins = NA, rep(NA,ncol(truth))))
  }else{
      pstPC <- pseudotime(pstPC)
      out <- cbind(out,
                   c(nlins = ncol(pstPC), agreement.pair(truth, pstPC)))
  }
  colnames(out)[ncol(out)] <- paste0('slingPC',p,'_MCbest')  
  
  # 4-D ICA - to match Monocle
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
  
  # 8-D Diffusion Map - to match DPT
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
}, mc.cores = NCORES)

# format output
nlins <- t(sapply(simResults_sling5, function(res){
  res[1,]
}))
taus <- lapply(1:ncol(simResults_sling5[[1]]),function(i){
  t(sapply(simResults_sling5, function(res){
    res[-1, i]
  }))
})
names(taus) <- colnames(simResults_sling5[[1]])

save(nlins,taus, file='results_sling5.RData')

