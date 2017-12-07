###################################################
### Lineage Inference Simulation Study
### Author: Kelly Street
### 
### This script will:
### 1) Read in the raw data from Trapnell, et al. (2014), as downloaded from
###    the conquer repository. These must be saved in the data/ directory.
### 2) Condense the data from transript-level to gene-level counts
### 3) Filter out the least informative genes and cells
### 4) Learn simulation parameters from the data
###################################################

# Read in and format the data
samples <- system('ls data/trapnell14/salmon/', intern=TRUE)

txCounts <- sapply(samples, function(samp){
  fname <- paste0('data/trapnell14/salmon/',samp,'/quant.sf')
  tab <- read.table(fname, header = TRUE)
  return(as.numeric(tab$NumReads))
})

rownames(txCounts) <- as.character(read.table(
  paste0('data/trapnell14/salmon/',samples[1],'/quant.sf'),
  header = TRUE)$Name)
rownames(txCounts) <- sapply(rownames(txCounts),function(n){
  unlist(strsplit(n,split='[.]'))[1]
})

# summarize to gene-level counts
library(EnsDb.Hsapiens.v86)
map <- transcripts(EnsDb.Hsapiens.v86, columns = c('tx_name','gene_name'),
                   return.type = 'data.frame')
txSub <- txCounts[rownames(txCounts) %in% map$tx_name,]
gene_names <- map$gene_name[match(rownames(txSub), map$tx_id)]
by_gene <- by(txSub, gene_names, colSums)
hsmmCounts <- t(simplify2array(by_gene))

# filter genes based on robust expression in at least 10% of cells
m <- quantile(hsmmCounts[hsmmCounts > 0], probs = .5)
gfilt <- apply(hsmmCounts,1,function(g){
  sum(g > m) >= .1 * ncol(hsmmCounts)
})
# filter cells: no more than 70% zero counts
sfilt <- apply(hsmmCounts[gfilt,],2,function(s){
  mean(s == 0) <= .7
})

# learn simulation parameters
parHSMM <- splatEstimate(round(hsmmCounts[gfilt,sfilt]))

save(parHSMM, file='data/parHSMM.RData')
