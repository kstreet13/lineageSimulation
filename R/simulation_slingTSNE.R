###################################################
### Lineage Inference Simulation Study
### Author: Kelly Street
###################################################
{
  library(splatter)
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

cells.range <- seq(40,500, by=20)
deProb.range <- (1:5)/10
reps <- 10
it.params <- cbind(rep(cells.range,each=reps*length(deProb.range)), rep(deProb.range,each=reps,times=length(cells.range)))
nIts <- nrow(it.params)
smallestOKcluster <- 5

library(mclust)
library(cluster)
library(slingshot)
library(parallel)
#library(destiny)
#library(fastICA)
library(Rtsne)
simResults_slingTSNE <- mclapply(1:nIts, function(i){
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
  pca <- prcomp(t(log1p(assays(sim)$normcounts)))$x
  # tSNE
  for(p in 2:8){
    # tSNE is slow, so we use the matrix of principle components
    tsne <- Rtsne::Rtsne(pca, dims = p)$Y
    rownames(tsne) <- rownames(pca)
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
}, mc.cores = 16)

save(simResults_slingTSNE, file='temp_slingTSNE.RData')

# format output
nlins <- t(sapply(simResults_slingTSNE, function(res){
  res[1,]
}))
taus <- lapply(1:ncol(simResults_slingTSNE[[1]]),function(i){
  t(sapply(simResults_slingTSNE, function(res){
    res[-1, i]
  }))
})
names(taus) <- colnames(simResults_slingTSNE[[1]])

save(nlins,taus, file='results_slingTSNE.RData')
