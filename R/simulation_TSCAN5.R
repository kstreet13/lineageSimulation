###################################################
### Lineage Inference Simulation Study
### TSCAN
### Five Lineages
###################################################

# setup
NCORES <- 1
source('R/helper_functions.R')
load('data/parHSMM.RData')
library(SingleCellExperiment)
library(parallel)
library(TSCAN)
library(slingshot)

# five-lineage case setup
cells.range <- seq(20,120, by=20)
deProb.range <- (1:5)/10
reps <- 10
it.params <- cbind(rep(cells.range,each=reps*length(deProb.range)), rep(deProb.range,each=reps,times=length(cells.range)))
nIts <- nrow(it.params)
smallestOKcluster <- 5

simResults_tscan5 <- mclapply(1:nIts, function(i){
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
}, mc.cores = NCORES)

# format output
nlins <- t(sapply(simResults_tscan5, function(res){
  res[1,]
}))
taus <- lapply(1:ncol(simResults_tscan5[[1]]),function(i){
  t(sapply(simResults_tscan5, function(res){
    res[-1, i]
  }))
})
names(taus) <- colnames(simResults_tscan5[[1]])

save(nlins,taus, file='results_tscan5.RData')
