###################################################
### Lineage Inference Simulation Study
### Monocle 2
### Two Lineages
###################################################

# setup
NCORES <- 1
source('R/helper_functions.R')
load('data/parHSMM.RData')
library(SingleCellExperiment)
# library(destiny)
library(monocle)
# library(TSCAN)
library(parallel)

# two-lineage case setup
cells.range <- seq(40,500, by=20)
deProb.range <- (1:5)/10
reps <- 10
it.params <- cbind(rep(cells.range,each=reps*length(deProb.range)), rep(deProb.range,each=reps,times=length(cells.range)))
nIts <- nrow(it.params)
smallestOKcluster <- 5

simResults_mon2_redo <- mclapply(1:nIts, function(i){
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
  
  vars <- apply(log1p(assays(sim)$normcounts),1,var)
  means <- rowMeans(log1p(assays(sim)$normcounts))
  ind <- (vars >= sort(vars,decreasing = TRUE)[5000]) | (means >= sort(means,decreasing = TRUE)[5000])
  
  for(D in 2:5){
    mon2_5k <- tryCatch({
      mon2_pst(sim[ind,], ndim = D)
    }, error = function(e){ e })
    if('error' %in% class(mon2_5k)){
      out <- cbind(out,
                   c(nlins = NA, rep(NA,ncol(truth))))
    }else{
      out <- cbind(out,
                   c(nlins = ncol(mon2_5k), agreement.pair(truth, mon2_5k)))
    }
    colnames(out)[ncol(out)] <- paste0('mon2_',D,'D5k')
  }
  for(D in 2:5){
    mon2_pc <- tryCatch({
      mon2_pst(sim, ndim = D, PCgenes = TRUE)
    }, error = function(e){ e })
    if('error' %in% class(mon2_pc)){
      out <- cbind(out,
                   c(nlins = NA, rep(NA,ncol(truth))))
    }else{
      out <- cbind(out,
                   c(nlins = ncol(mon2_pc), agreement.pair(truth, mon2_pc)))
    }
    colnames(out)[ncol(out)] <- paste0('mon2_',D,'Dpc')
  }

  return(out)
}, mc.cores = NCORES)

# format output
nlins <- t(sapply(simResults_mon2, function(res){
  res[1,]
}))
taus <- lapply(1:ncol(simResults_mon2[[1]]),function(i){
  t(sapply(simResults_mon2, function(res){
    res[-1, i]
  }))
})
names(taus) <- colnames(simResults_mon2[[1]])

save(nlins,taus, file='results_mon2.RData')
