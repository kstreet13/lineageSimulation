###################################################
### Lineage Inference Simulation Study
### Author: Kelly Street
###################################################
{
  # library(splatter)
  source('code/helper_functions.R')
  # samples <- system('ls data/trapnell14/salmon/', intern=TRUE)
  # txCounts <- sapply(samples, function(samp){
  #   fname <- paste0('data/trapnell14/salmon/',samp,'/quant.sf')
  #   tab <- read.table(fname, header = TRUE)
  #   return(as.numeric(tab$NumReads))
  # })
  # rownames(txCounts) <- as.character(read.table(
  #   paste0('data/trapnell14/salmon/',samples[1],'/quant.sf'),
  #   header = TRUE)$Name)
  # rownames(txCounts) <- sapply(rownames(txCounts),function(n){
  #   unlist(strsplit(n,split='[.]'))[1]
  # })
  # 
  # library(EnsDb.Hsapiens.v86)
  # map <- transcripts(EnsDb.Hsapiens.v86, columns = c('tx_name','gene_name'),
  #                    return.type = 'data.frame')
  # 
  # txSub <- txCounts[rownames(txCounts) %in% map$tx_name,]
  # gene_names <- map$gene_name[match(rownames(txSub), map$tx_id)]
  # by_gene <- by(txSub, gene_names, colSums)
  # hsmmCounts <- t(simplify2array(by_gene))
  # 
  # m <- quantile(hsmmCounts[hsmmCounts > 0], probs = .5)
  # # filter genes based on robust expression in at least 10% of cells
  # gfilt <- apply(hsmmCounts,1,function(g){
  #   sum(g > m) >= .1 * ncol(hsmmCounts)
  # })
  # # filter cells: no more than 70% zero counts
  # sfilt <- apply(hsmmCounts[gfilt,],2,function(s){
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
# library(destiny)
library(monocle)
# library(TSCAN)
library(parallel)

simResults_mon1 <- mclapply(1:nIts, function(i){
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
    mon1_5k <- tryCatch({
      mon1_pst(sim[ind,], ndim = D, num_paths = 2)
    }, error = function(e){ e })
    if('error' %in% class(mon1_5k)){
      out <- cbind(out,
                   c(nlins = NA, rep(NA,ncol(truth))))
    }else{
      out <- cbind(out,
                   c(nlins = ncol(mon1_5k), agreement.pair(truth, mon1_5k)))
    }
    colnames(out)[ncol(out)] <- paste0('mon1_',D,'D5k')
  }
  for(D in 2:5){
    mon1_pc <- tryCatch({
      mon1_pst(sim, ndim = D, PCgenes = TRUE, num_paths = 2)
    }, error = function(e){ e })
    if('error' %in% class(mon1_pc)){
      out <- cbind(out,
                   c(nlins = NA, rep(NA,ncol(truth))))
    }else{
      out <- cbind(out,
                   c(nlins = ncol(mon1_pc), agreement.pair(truth, mon1_pc)))
    }
    colnames(out)[ncol(out)] <- paste0('mon1_',D,'Dpc')
  }

  return(out)
}, mc.cores = 16)

save(simResults_mon1, file='temp_mon1.RData')

# format output
nlins <- t(sapply(simResults_mon1, function(res){
  res[1,]
}))
taus <- lapply(1:ncol(simResults_mon1[[1]]),function(i){
  t(sapply(simResults_mon1, function(res){
    res[-1, i]
  }))
})
names(taus) <- colnames(simResults_mon1[[1]])

save(nlins,taus, file='results_mon1.RData')
