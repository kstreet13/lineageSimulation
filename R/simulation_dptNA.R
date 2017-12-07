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

cells.range <- seq(40,500, by=20)
deProb.range <- (1:5)/10
reps <- 10
it.params <- cbind(rep(cells.range,each=reps*length(deProb.range)), rep(deProb.range,each=reps,times=length(cells.range)))
nIts <- nrow(it.params)
smallestOKcluster <- 5

library(SingleCellExperiment)
library(destiny)
# library(monocle)
# library(TSCAN)
library(parallel)

# i = 40
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
}, mc.cores = 16)

save(simResults_dpt, file='temp_dptNA.RData')

# format output
dpt.na <- t(simplify2array(simResults_dpt))

save(dpt.na, file='results_dptNA.RData')
