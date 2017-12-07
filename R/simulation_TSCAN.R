###################################################
### Lineage Inference Simulation Study
### TSCAN
### Two Lineages
###################################################

# setup
source('R/helper_functions.R')
load('data/parHSMM.RData')
library(SingleCellExperiment)
library(parallel)
library(TSCAN)
library(slingshot)

# two-lineage case setup
cells.range <- seq(40,500, by=20)
deProb.range <- (1:5)/10
reps <- 10
it.params <- cbind(rep(cells.range,each=reps*length(deProb.range)), rep(deProb.range,each=reps,times=length(cells.range)))
nIts <- nrow(it.params)
smallestOKcluster <- 5

simResults_TSCAN <- mclapply(1:nIts, function(i){
  print(paste('Iteration',i))
  {
    ncells <- it.params[i,1]
    deProb <- it.params[i,2]
    splat.it.params <- splatter::setParams(parHSMM, update=list(batchCells=3*ncells, group.prob=c(1,1,1)/3, 
                                                                path.nonlinearProb=.5, de.prob=deProb, 
                                                                de.facLoc=.25, de.facScale=1,
                                                                path.length=1000, path.from=c(0,1,1), 
                                                                dropout.present = TRUE, seed = i))
    sim <- splatter::splatSimulate(splat.it.params, method='paths')
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
  
  ### Slingshot w/TSCAN dimensionality reduction and clustering
  
  tscanSling <- tryCatch({
    tscan_sling_pst(log1p(assays(sim)$normcounts), 
                    smallestOKcluster = 5, preproc = FALSE,
                    start.cell = colnames(sim)[which.min(sim$Step)])
  }, error = function(e){ e })
  if('error' %in% class(tscanSling)){
    out <- cbind(out,
                 c(nlins = NA, rep(NA,ncol(truth))))
  }else{
    out <- cbind(out,
                 c(nlins = ncol(tscanSling), agreement.pair(truth, tscanSling)))
  }
  colnames(out)[ncol(out)] <- 'tscanSling'
  
  ### TSCAN (FQ normalization, no preprocessing)
  tscanNorm <- tryCatch({
    tscan_pst(log1p(assays(sim)$normcounts), smallestOKcluster = 5, 
              preproc = FALSE,
              start.cell = colnames(sim)[which.min(sim$Step)])
  }, error = function(e){ e })
  if('error' %in% class(tscanNorm)){
    out <- cbind(out,
                 c(nlins = NA, rep(NA,ncol(truth))))
  }else{
    out <- cbind(out,
                 c(nlins = ncol(tscanNorm), agreement.pair(truth, tscanNorm)))
  }
  colnames(out)[ncol(out)] <- 'tscanNorm'
  
  
  ### TSCAN (no normalization, preprocessing)
  tscanPrep <- tryCatch({
    tscan_pst(assays(sim)$counts, smallestOKcluster = 5, 
              preproc = TRUE,
              start.cell = colnames(sim)[which.min(sim$Step)])
  }, error = function(e){ e })
  if('error' %in% class(tscanPrep)){
    out <- cbind(out,
                 c(nlins = NA, rep(NA,ncol(truth))))
  }else{
    out <- cbind(out,
                 c(nlins = ncol(tscanPrep), agreement.pair(truth, tscanPrep)))
  }
  colnames(out)[ncol(out)] <- 'tscanPrep'
  
  ### TSCAN (FQ normalization, preprocessing)
  tscanNormPrep <- tryCatch({
    tscan_pst(assays(sim)$normcounts, smallestOKcluster = 5, 
              preproc = TRUE,
              start.cell = colnames(sim)[which.min(sim$Step)])
  }, error = function(e){ e })
  if('error' %in% class(tscanNormPrep)){
    out <- cbind(out,
                 c(nlins = NA, rep(NA,ncol(truth))))
  }else{
    out <- cbind(out,
                 c(nlins = ncol(tscanNormPrep), agreement.pair(truth, tscanNormPrep)))
  }
  colnames(out)[ncol(out)] <- 'tscanNormPrep'
  
  return(out)
}, mc.cores = 1)

# format output
nlins <- t(sapply(simResults_TSCAN, function(res){
  res[1,]
}))
taus <- lapply(1:ncol(simResults_TSCAN[[1]]),function(i){
  t(sapply(simResults_TSCAN, function(res){
    res[-1, i]
  }))
})
names(taus) <- colnames(simResults_TSCAN[[1]])

save(nlins,taus, file='results_TSCAN.RData')
