###################################################
### Lineage Inference Simulation Study
### DPT - highest-level branching event info
### Two Lineages
###################################################

# setup
NCORES <- 1
source('R/helper_functions.R')
load('data/parHSMM.RData')
library(SingleCellExperiment)
library(destiny)
# library(monocle)
# library(TSCAN)
library(parallel)

# two-lineage case setup
cells.range <- seq(40,500, by=20)
deProb.range <- (1:5)/10
reps <- 10
it.params <- cbind(rep(cells.range,each=reps*length(deProb.range)), rep(deProb.range,each=reps,times=length(cells.range)))
nIts <- nrow(it.params)
smallestOKcluster <- 5

simResults_dpt <- mclapply(1:nIts, function(i){
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
  
  start.cell <- colnames(sim)[which.min(sim$Step)]
  start.ind <- which(colnames(sim)==start.cell)
  
  dm <- DiffusionMap(t(log1p(assays(sim)$normcounts)))
  dpt <- DPT(dm, tips = start.ind)
  
  bm <- dpt@branch

  return(table(bm[,1], useNA = "always"))  
}, mc.cores = NCORES)

# format output
dpt.na <- t(simplify2array(simResults_dpt))

save(dpt.na, file='results_dptNA.RData')
