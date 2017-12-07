###################################################
### Lineage Inference Simulation Study
### Cluster MST
### Two Lineages
###################################################

# setup
NCORES <- 1
source('R/helper_functions.R')
load('data/parHSMM.RData')
library(SingleCellExperiment)
library(parallel)
library(mclust)

# two-lineage case setup
cells.range <- seq(40,500, by=20)
deProb.range <- (1:5)/10
reps <- 10
it.params <- cbind(rep(cells.range,each=reps*length(deProb.range)), rep(deProb.range,each=reps,times=length(cells.range)))
nIts <- nrow(it.params)
smallestOKcluster <- 5

simResults_mstClus <- mclapply(1:nIts, function(i){
  print(paste('Iteration',i))
  {
    ncells <- it.params[i,1]
    deProb <- it.params[i,2]
    splat.it.params <- setParams(parHSMM, update=list(batchCells=3*ncells, group.prob=c(1,1,1)/3, 
                                                      path.nonlinearProb=.5, de.prob=deProb, 
                                                      de.facLoc=.25, de.facScale=1,
                                                      path.length=1000, path.from=c(0,1,1), 
                                                      dropout.present = TRUE, seed = i))
    sim <- splatSimulate(splat.it.params, method='paths')
    sim$Step[sim$Group %in% c('Path2','Path3')] <- sim$Step[sim$Group %in% c('Path2','Path3')] + 1000
    assays(sim)$normcounts <- FQnorm(assays(sim)$counts)
    start.cell <- colnames(sim)[which.min(sim$Step)]
    truth <- cbind(sim$Step,sim$Step)
    truth[sim$Group == 'Path3', 1] <- NA
    truth[sim$Group == 'Path2', 2] <- NA
    rownames(truth) <- colnames(sim)
    colnames(truth) <- c('Lineage1','Lineage2')
  }# simulate, normalize data
  out <- NULL
  
  # fix dimensionality reduction
  X <- prcomp(t(log1p(assays(sim)$normcounts)))$x[,1:3]
  X.orig <- X
  
  {
    for(K in 3:10){
      X <- X.orig
      mc <- Mclust(X, G = K)
      cl <- mc$classification
      minClusSize <- min(table(cl))
      while(minClusSize < smallestOKcluster){
        clTab <- table(cl)
        smallClus <- names(clTab)[clTab < smallestOKcluster]
        X <- X[! cl %in% smallClus, ]
        mc <- Mclust(X, G = K)
        cl <- mc$classification
        minClusSize <- min(table(cl))
        if(nrow(X) < .5*nrow(X.orig)){ # prevent dropping too many cells
          X <- X.orig
          mc <- Mclust(X, G = K)
          cl <- mc$classification
          minClusSize <- smallestOKcluster + 1
        }
      }
      clus1 <- as.character(cl[which.max(names(cl)==start.cell)])
      
      mst.kMC <- tryCatch({
        mst_pst(X, cl, clus1)
      }, error = function(e){ e })
      if('error' %in% class(mst.kMC)){
        out <- cbind(out,
                     c(nlins = NA, rep(NA,ncol(truth))))
      }else{
        out <- cbind(out,
                     c(nlins = ncol(mst.kMC), agreement.pair(truth, mst.kMC)))
      }
      colnames(out)[ncol(out)] <- paste0('MST_MC',K)
    }
  } # MC_k
  #############
  {
    for(K in 3:10){
      X <- X.orig
      km <- kmeans(X, centers = K, nstart = 20)
      cl <- km$cluster
      minClusSize <- min(table(cl))
      while(minClusSize < smallestOKcluster){
        clTab <- table(cl)
        smallClus <- names(clTab)[clTab < smallestOKcluster]
        X <- X[! cl %in% smallClus, ]
        km <- kmeans(X, centers = K, nstart = 20)
        cl <- km$cluster
        minClusSize <- min(table(cl))
        if(nrow(X) < .5*nrow(X.orig)){ # prevent dropping too many cells
          X <- X.orig
          km <- kmeans(X, centers = K, nstart = 20)
          cl <- km$cluster
          minClusSize <- smallestOKcluster + 1
        }
      }
      clus1 <- as.character(cl[which.max(names(cl)==start.cell)])
      
      mst.kKM <- tryCatch({
        mst_pst(X, cl, clus1)
      }, error = function(e){ e })
      if('error' %in% class(mst.kMC)){
        out <- cbind(out,
                     c(nlins = NA, rep(NA,ncol(truth))))
      }else{
        out <- cbind(out,
                     c(nlins = ncol(mst.kKM), agreement.pair(truth, mst.kKM)))
      }
      colnames(out)[ncol(out)] <- paste0('MST_KM',K)
    }
  } # KM_k
  #############
  {
    for(K in 3:10){
      X <- X.orig
      d <- dist(X)
      h <- hclust(d)
      cl <- cutree(h, k = K)
      minClusSize <- min(table(cl))
      while(minClusSize < smallestOKcluster){
        clTab <- table(cl)
        smallClus <- names(clTab)[clTab < smallestOKcluster]
        X <- X[! cl %in% smallClus, ]
        d <- dist(X)
        h <- hclust(d)
        cl <- cutree(h, k = K)
        minClusSize <- min(table(cl))
        if(nrow(X) < .5*nrow(X.orig)){ # prevent dropping too many cells
          X <- X.orig
          d <- dist(X)
          h <- hclust(d)
          cl <- cutree(h, k = K)
          minClusSize <- smallestOKcluster + 1
        }
      }
      clus1 <- as.character(cl[which.max(names(cl)==start.cell)])
      
      mst.kHC <- tryCatch({
        mst_pst(X, cl, clus1)
      }, error = function(e){ e })
      if('error' %in% class(mst.kHC)){
        out <- cbind(out,
                     c(nlins = NA, rep(NA,ncol(truth))))
      }else{
        out <- cbind(out,
                     c(nlins = ncol(mst.kHC), agreement.pair(truth, mst.kHC)))
      }
      colnames(out)[ncol(out)] <- paste0('MST_HC',K)
    }
  } # HC_k
  
  return(out)
}, mc.cores = NCORES)

# format output
nlins <- t(sapply(simResults_mstClus, function(res){
  res[1,]
}))
taus <- lapply(1:ncol(simResults_mstClus[[1]]),function(i){
  t(sapply(simResults_mstClus, function(res){
    res[-1, i]
  }))
})
names(taus) <- colnames(simResults_mstClus[[1]])

save(nlins,taus, file='results_mstClus.RData')
