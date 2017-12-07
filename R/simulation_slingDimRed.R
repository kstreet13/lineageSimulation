###################################################
### Lineage Inference Simulation Study
### Slingshot - robustness to clustering
### Two Lineages
###################################################

# setup
NCORES <- 1
source('R/helper_functions.R')
load('data/parHSMM.RData')
library(SingleCellExperiment)
library(parallel)
library(slingshot)
library(mclust)
library(destiny)
library(fastICA)
library(Rtsne)

# two-lineage case setup
cells.range <- seq(40,500, by=20)
deProb.range <- (1:5)/10
reps <- 10
it.params <- cbind(rep(cells.range,each=reps*length(deProb.range)), rep(deProb.range,each=reps,times=length(cells.range)))
nIts <- nrow(it.params)
smallestOKcluster <- 5

simResults_slingDimRed <- mclapply(1:nIts, function(i){
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
  out <- list()
  
  # PCA
  pca <- prcomp(t(log1p(assays(sim)$normcounts)))
  for(p in 2:8){
    sds <- newSlingshotDataSet(pca$x[,1:p])
    out[[length(out)+1]] <- sds
    names(out)[length(out)] <- paste0('pca',p)
  }
  # ICA
  for(p in 2:8){
    # ICA is slow and large matrices can cause errors in the underlying
    # Fortran code, so we do ICA on the full principal components matrix
    ica <- fastICA::fastICA(pca$x, n.comp = p, method = 'C')$S
    rownames(ica) <- rownames(pca$x)
    sds <- newSlingshotDataSet(ica)
    out[[length(out)+1]] <- sds
    names(out)[length(out)] <- paste0('ica',p)
  }
  # Diffusion Map
  dm <- destiny::DiffusionMap(t(log1p(assays(sim)$normcounts)))@eigenvectors
  rownames(dm) <- rownames(pca$x)
  for(p in 2:8){
    sds <- newSlingshotDataSet(dm[,1:p])
    out[[length(out)+1]] <- sds
    names(out)[length(out)] <- paste0('dm',p)
  }
  # tSNE
  for(p in 2:8){
    # tSNE is similarly slow, so we use the matrix of principle components
    tsne <- Rtsne::Rtsne(pca$x, dims = p)$Y
    rownames(tsne) <- rownames(pca$x)
    sds <- newSlingshotDataSet(tsne)
    out[[length(out)+1]] <- sds
    names(out)[length(out)] <- paste0('tsne',p)
  }
  
  #######################################  
  
  # clustering 1: bestMC
  # Gaussian mixture modelling with best K chosen by BIC
  res <- sapply(out,function(sds){
    X <- reducedDim(sds)
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
    
    
    pstDR <- tryCatch({
      slingshot(X, cl, start.clus = clus1)
    }, error = function(e){ e })
    
    if('error' %in% class(pstDR)){
      ret <- c(nlins = NA, rep(NA,ncol(truth)))
    }else{
      pstDR <- pseudotime(pstDR)
      #rownames(pstDR) <- colnames(sim)
      ret <- c(nlins = ncol(pstDR), agreement.pair(truth, pstDR))
    }
    return(ret)
  })
  
  return(res)
}, mc.cores = NCORES)

# format output
nlins <- t(sapply(simResults_slingDimRed, function(res){
  res[1,]
}))
taus <- lapply(1:ncol(simResults_slingDimRed[[1]]),function(i){
  t(sapply(simResults_slingDimRed, function(res){
    res[-1, i]
  }))
})
names(taus) <- colnames(simResults_slingDimRed[[1]])

save(nlins,taus, file='results_slingDimRed.RData')
