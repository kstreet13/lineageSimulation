
### Pseudotime functions

slingshot_pst <- function(redX, clus.labels, ...){
  require(slingshot)
  l <- get_lineages(redX, clus.labels, ...)
  c <- get_curves(redX, clus.labels, l)
  out <- sapply(c,function(i){i$pseudotime})
  return(out)
}

tscan_pst <- function(counts, preproc = FALSE, clusNum = 2:9,
                      smallestOKcluster = 1, start.cell){
  require(TSCAN)
  require(igraph)
  if(preproc){
    counts <- preprocess(counts)
  }
  counts <- counts[apply(counts,1,var) > 0, ]
  counts.orig <- counts
  clust <- exprmclust(counts, clusternum = clusNum)
  
  minClusSize <- min(table(clust$clusterid))
  while(minClusSize < smallestOKcluster){
    clTab <- table(clust$clusterid)
    smallClus <- names(clTab)[clTab < smallestOKcluster]
    counts <- counts[, ! clust$clusterid %in% smallClus]
    counts <- counts[apply(counts,1,var) > 0, ]
    clust <- exprmclust(counts, clusternum = clusNum)
    minClusSize <- min(table(clust$clusterid))
    if(ncol(counts) < .5*ncol(counts.orig)){
      counts <- counts.orig
      clust <- exprmclust(counts, clusternum = clusNum)
      minClusSize <- smallestOKcluster + 1
    }
  }
  
  #start.cell <- colnames(sim)[which.min(sim$Step)]
  clus1 <- clust$clusterid[names(clust$clusterid)==start.cell]
  
  # define lineages as in Slingshot, ordered list of clusters
  lineages <- list()
  g <- clust$MSTtree
  tree <- as.matrix(as_adj(g))
  degree <- rowSums(tree)
  ends <- rownames(tree)[degree == 1 & rownames(tree) != clus1]
  paths <- shortest_paths(g, from = clus1, to = ends, mode = 'out', output = 'vpath')$vpath
  for(p in paths){
    lineages[[length(lineages)+1]] <- names(p)
  }
  
  pst <- sapply(lineages,function(l){
    ord <- TSCANorder(clust, MSTorder = as.numeric(l), orderonly = FALSE)
    return(ord$Pseudotime[match(colnames(counts.orig), ord$sample_name)])
  })
  rownames(pst) <- colnames(counts.orig)
  colnames(pst) <- paste0('Lineage',1:length(lineages))
  
  return(pst)
}

tscan_sling_pst <- function(counts, preproc = FALSE, clusNum = 2:9, 
                            smallestOKcluster = 1, start.cell){
  require(TSCAN)
  require(igraph)
  require(slingshot)
  if(preproc){
    counts <- preprocess(counts)
  }
  counts <- counts[apply(counts,1,var) > 0, ]
  counts.orig <- counts
  clust <- exprmclust(counts, clusternum = clusNum)
  
  minClusSize <- min(table(clust$clusterid))
  while(minClusSize < smallestOKcluster){
    print(ncol(counts))
    clTab <- table(clust$clusterid)
    smallClus <- names(clTab)[clTab < smallestOKcluster]
    counts <- counts[, ! clust$clusterid %in% smallClus]
    counts <- counts[apply(counts,1,var) > 0, ]
    clust <- exprmclust(counts, clusternum = clusNum)
    minClusSize <- min(table(clust$clusterid))
    if(ncol(counts) < .5*ncol(counts.orig)){
      counts <- counts.orig
      clust <- exprmclust(counts, clusternum = clusNum)
      minClusSize <- smallestOKcluster + 1
    }
  }
  
  #start.cell <- colnames(sim)[which.min(sim$Step)]
  clus1 <- clust$clusterid[names(clust$clusterid)==start.cell]
  
  sds <- slingshot(clust$pcareduceres, clust$clusterid, start.clus = clus1)
  pst <- as.matrix(pseudotime(sds))
  pst <- pst[match(colnames(counts.orig), rownames(pst)), ,drop = FALSE]
  rownames(pst) <- colnames(counts.orig)
  return(pst)
}

mon2_pst <- function(sim, ndims=2, PCgenes = FALSE){
  require(monocle)
  require(igraph)
  start.cell <- colnames(sim)[which.min(sim$Step)]
  rowData(sim)$gene_short_name <- rownames(sim)
  cds <- newCellDataSet(cellData = assays(sim)$normcounts, phenoData = AnnotatedDataFrame(as.data.frame(colData(sim))), featureData = AnnotatedDataFrame(as.data.frame(rowData(sim),row.names=rownames(sim))))
  cds <- estimateSizeFactors(cds)
  if(PCgenes){
    ordering_genes <- ordering_genes_pca(cds, which.pcs = seq_len(ndims))
  }else{
    ordering_genes <- rownames(cds)
  }
  cds <- setOrderingFilter(cds, ordering_genes)
  if(ncol(sim)==100){
    ncenter <- 99
  }else{
    ncenter <- NULL
  }
  cds <- reduceDimension(cds, reduction_method = 'DDRTree', ncenter = ncenter, max_components = ndims)
  cds <- orderCells(cds)
  state1 <- as.character(cds$State[which.max(colnames(cds)==start.cell)])
  if(length(cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell)==0){
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- rownames(pData(cds))[which.min(cds$Pseudotime)]
  }
  cds <- orderCells(cds, root_state = state1)
  if(length(cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell)==0){
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- rownames(pData(cds))[which.min(cds$Pseudotime)]
  }
  
  state1 <- as.character(cds$State[which.max(colnames(cds)==cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell)])
  nBranches <- length(cds@auxOrderingData$DDRTree$branch_points)
  
  if(nBranches > 0){
    nonStartStates <- unique(cds$State)
    nonStartStates <- nonStartStates[nonStartStates != state1]
    
    linInds <- lapply(nonStartStates,function(leafState){
      index <- tryCatch(
        getLineageID(cds,leafState),
        error=function(e) e
      )
      if(inherits(index, "error")) return(NULL)
      return(index)
    })
    linInds <- do.call("cbind", linInds)
    # #remove column if any other column is a subset
    # keep <- sapply(1:ncol(linInds),function(i){
    #   all(nonsubset <- sapply(1:ncol(linInds),function(j){
    #     if(i==j){
    #       return(TRUE)
    #     }
    #     if(all(linInds[linInds[,j],i])){
    #       return(FALSE)
    #     }else{
    #       return(TRUE)
    #     }
    #   }))
    # })
    # linInds <- linInds[,keep]
    pseudotime <- matrix(rep(cds$Pseudotime,ncol(linInds)), ncol = ncol(linInds))
    pseudotime[!linInds] <- NA
    ncellsLin <- apply(pseudotime,2,function(x){sum(!is.na(x))})
    pseudotime <- pseudotime[,order(ncellsLin, decreasing = TRUE)]
    lins <- ncol(pseudotime)
  }else{
    pseudotime <- matrix(cds$Pseudotime, ncol = 1)
    lins <- 1
  }
  return(pseudotime)
}

mon1_pst <- function(sim, ndim = 2, num_paths = NULL, PCgenes = FALSE){
  require(monocle)
  require(igraph)
  start.cell <- colnames(sim)[which.min(sim$Step)]
  rowData(sim)$gene_short_name <- rownames(sim)
  cds <- newCellDataSet(cellData = assays(sim)$normcounts, phenoData = AnnotatedDataFrame(as.data.frame(colData(sim))), featureData = AnnotatedDataFrame(as.data.frame(rowData(sim),row.names=rownames(sim))))
  cds <- estimateSizeFactors(cds)
  
  if(PCgenes){
    ordering_genes <- ordering_genes_pca(cds, which.pcs = seq_len(ndim))
  }else{
    ordering_genes <- rownames(cds)
  }
  cds <- setOrderingFilter(cds, ordering_genes)
  cds <- estimateSizeFactors(cds)
  
  cds <- reduceDimension(cds, max_components = ndim, reduction_method = "ICA")
  
  cds <- orderCells(cds, num_paths = num_paths)
  
  state1 <- as.character(cds$State[which.max(colnames(cds)==start.cell)])

  cds <- orderCells(cds, num_paths = num_paths, root_state = state1)

  # formatting output as n x L matrix
  # (states = clusters)
  clusLabels <- as.character(cds$State)
  clusters <- sort(unique(clusLabels))
  cellCon <- get.adjacency(cds@minSpanningTree)
  #cluster connectivity matrix
  clusterCon <- sapply(clusters,function(c1){
    sapply(clusters,function(c2){
      if(c1==c2){
        return(0)
      }
      ind1 <- clusLabels == c1
      ind2 <- clusLabels == c2
      return(max(cellCon[ind1,ind2]))
    })
  })
  g <- graph.adjacency(clusterCon, mode="undirected")
  start <- clusLabels[which.min(cds$Pseudotime)]
  ends <- clusters[(length(clusters)-num_paths+1):length(clusters)]
  paths <- shortest_paths(g, from = start, to = ends, mode = 'out', output = 'vpath')$vpath
  # define lineages as ordered sets of clusters
  lineages <- lapply(paths,function(p){
    lin <- names(p)
    # each lineage should have only one terminal state
    # mst geometry can be a little weird, can't rely on 
    # terminal states to have only one connection
    terminal <- lin[length(lin)]
    lin <- lin[! lin %in% ends[ends != terminal]]
    return(lin)
  })
  # extract pseudotime for each lineage
  pst <- sapply(lineages,function(lin){
    ind <- clusLabels %in% lin
    pst <- cds$Pseudotime
    pst[! ind] <- NA
    return(pst)
  })
  rownames(pst) <- colnames(cds)
  colnames(pst) <- paste0('Lineage',1:ncol(pst))
  return(pst)
}

mst_pst <- function(X, clusterLabels, clus1){
    require(slingshot)
    sds <- getLineages(X, clusterLabels, start.clus = clus1)
    lineages <- sds@lineages
    
    L <- length(grep("Lineage",names(lineages))) # number of lineages
    clusters <- unique(clusterLabels)
    d <- dim(X); n <- d[1]; p <- d[2]
    nclus <- length(clusters)
    centers <- t(sapply(clusters,function(clID){
        x.sub <- X[clusterLabels == clID, ,drop = FALSE]
        return(colMeans(x.sub))
    }))
    if(p == 1){
        centers <- t(centers)
    }
    rownames(centers) <- clusters
    W <- sapply(seq_len(L),function(l){
        as.numeric(clusterLabels %in% lineages[[l]])
    }) # weighting matrix
    rownames(W) <- rownames(X)
    colnames(W) <- names(lineages)[seq_len(L)]
    W.orig <- W
    D <- W; D[,] <- NA
    
    # determine curve hierarchy
    C <- as.matrix(sapply(lineages[seq_len(L)], function(lin) {
        sapply(clusters, function(clID) {
            as.numeric(clID %in% lin)
        })
    }))
    rownames(C) <- clusters
    segmnts <- unique(C[rowSums(C)>1,,drop = FALSE])
    segmnts <- segmnts[order(rowSums(segmnts),decreasing = FALSE), ,
                       drop = FALSE]
    avg.order <- list()
    for(i in seq_len(nrow(segmnts))){
        idx <- segmnts[i,] == 1
        avg.order[[i]] <- colnames(segmnts)[idx]
        new.col <- rowMeans(segmnts[,idx, drop = FALSE])
        segmnts <- cbind(segmnts[, !idx, drop = FALSE],new.col)
        colnames(segmnts)[ncol(segmnts)] <- paste('average',i,sep='')
    }
    
    # initial curves are piecewise linear paths through the tree
    pcurves <- list()
    for(l in seq_len(L)){
        idx <- W[,l] > 0
        clus.sub <- clusterLabels[idx]
        line.initial <- centers[clusters %in% lineages[[l]], , 
                                drop = FALSE]
        line.initial <- line.initial[match(lineages[[l]],
                                           rownames(line.initial)),  ,
                                     drop = FALSE]
        K <- nrow(line.initial)
        # special case: single-cluster lineage
        if(K == 1){
            pca <- prcomp(X[idx, ,drop = FALSE])
            ctr <- line.initial
            line.initial <- rbind(ctr - 10*pca$sdev[1] * 
                                      pca$rotation[,1], ctr, 
                                  ctr + 10*pca$sdev[1] *
                                      pca$rotation[,1])
            curve <- get.lam(X[idx, ,drop = FALSE], s = line.initial,
                             stretch = 9999)
            # do this twice because all points should have projections
            # on all lineages, but only those points on the lineage
            # should extend it
            pcurve <- get.lam(X, s = curve$s[curve$tag,], stretch=0)
            pcurve$lambda <- pcurve$lambda - min(pcurve$lambda, 
                                                 na.rm=TRUE)
            # ^ force pseudotime to start at 0
            pcurve$w <- W[,l]
            pcurves[[l]] <- pcurve
            D[,l] <- abs(pcurve$dist)
            next
        }
        
        curve <- get.lam(X[idx, ,drop = FALSE], s = line.initial,
                         stretch = 9999)
        
        pcurve <- get.lam(X, s = curve$s[curve$tag, ,drop=FALSE], 
                          stretch=0)
        # force pseudotime to start at 0
        pcurve$lambda <- pcurve$lambda - min(pcurve$lambda, 
                                             na.rm=TRUE) 
        pcurve$w <- W[,l]
        pcurves[[l]] <- pcurve
        D[,l] <- abs(pcurve$dist)
    }
    pst <- sapply(pcurves,function(x){
        x$lambda
    })
    pst[W==0] <- NA
    colnames(pst) <- paste0('Lineage',1:ncol(pst))
    rownames(pst) <- rownames(X)
    return(pst)
}

dpt_pst <- function(sim, num_paths = 'full'){
  num_paths <- as.character(num_paths)
  start.cell <- colnames(sim)[which.min(sim$Step)]
  start.ind <- which(colnames(sim)==start.cell)
  
  dm <- DiffusionMap(t(log1p(assays(sim)$normcounts)))
  dpt <- DPT(dm, tips = start.ind)
  
  pst.vec <- dpt[start.ind,]
  bm <- dpt@branch
  
  if(num_paths=='full'){
    if(sum(bm[,1] %in% 2:3) == 0){
      # one lineage
      pst <- as.matrix(pst.vec, ncol=1)
      rownames(pst) <- colnames(sim)
      colnames(pst) <- paste0('Lineage',1:ncol(pst))
      return(pst)
    }
    lineages <- list(c(1,2), c(1,3))
    ends <- 2:3
    if(max(bm,na.rm = TRUE)==3){
      done <- TRUE
    }else{
      done <- FALSE
    }
    while(! done){
      length.old <- length(lineages)
      ends <- sapply(lineages, function(lin){ lin[length(lin)] })
      ends <- ends[!is.na(ends)]
      endcols <- sapply(ends,function(end){
        which(sapply(1:ncol(bm),function(j){ end %in% bm[,j]}))
      })
      to.split <- ends[endcols < ncol(bm)]
      to.split <- to.split[order(table(bm)[to.split], decreasing = TRUE)]
      
      for(end in to.split){
        endcol <- which(sapply(1:ncol(bm),function(j){ end %in% bm[,j]}))
        
        new <- unique(bm[bm[,endcol]==end, endcol+1])
        new <- new[! is.na(new)]
        if(length(new)==0){
          next
        }
        lin <- lineages[[which(ends==end)]]
        lin <- lin[-length(lin)]
        lineages <- lineages[ends != end]
        new.ends <- new[new %% 3 != 1]
        new.root <- new[new %% 3 == 1]
        lineages[[length(lineages)+1]] <- c(lin, new.root, new.ends[1])
        lineages[[length(lineages)+1]] <- c(lin, new.root, new.ends[2])
        break
      }
      if(length(lineages)==length.old){
        done <- TRUE
      }
    }
    linInds <- lapply(lineages,function(lin){
      apply(bm,1,function(bm.i){
        any(bm.i %in% lin)
      })
    })
    pst <- sapply(linInds,function(ind){
      out <- pst.vec
      out[! ind] <- NA
      return(out)
    })
    rownames(pst) <- colnames(sim)
    colnames(pst) <- paste0('Lineage',1:ncol(pst))
    return(pst)
  }
  if(as.numeric(num_paths) == 1){
    pst <- as.matrix(pst.vec, ncol=1)
    rownames(pst) <- colnames(sim)
    colnames(pst) <- paste0('Lineage',1:ncol(pst))
    return(pst)
  }
  if(as.numeric(num_paths) %in% 2:20){
    lineages <- list(c(1,2), c(1,3))
    ends <- 2:3
    if(as.numeric(num_paths)==2 | max(bm,na.rm = TRUE)==3){
      done <- TRUE
    }else{
      done <- FALSE
    }
    while(! done){
      length.old <- length(lineages)
      ends <- sapply(lineages, function(lin){ lin[length(lin)] })
      ends <- ends[!is.na(ends)]
      endcols <- sapply(ends,function(end){
        which(sapply(1:ncol(bm),function(j){ end %in% bm[,j]}))
      })
      to.split <- ends[endcols < ncol(bm)]
      to.split <- to.split[order(table(bm)[to.split], decreasing = TRUE)]
      
      for(end in to.split){
        endcol <- which(sapply(1:ncol(bm),function(j){ end %in% bm[,j]}))
        
        new <- unique(bm[bm[,endcol]==end, endcol+1])
        new <- new[! is.na(new)]
        if(length(new)==0){
          next
        }
        lin <- lineages[[which(ends==end)]]
        lin <- lin[-length(lin)]
        lineages <- lineages[ends != end]
        new.ends <- new[new %% 3 != 1]
        new.root <- new[new %% 3 == 1]
        lineages[[length(lineages)+1]] <- c(lin, new.root, new.ends[1])
        lineages[[length(lineages)+1]] <- c(lin, new.root, new.ends[2])
        break
      }
      if(length(lineages)==length.old | length(lineages) == num_paths){
        done <- TRUE
      }
    }
    linInds <- lapply(lineages,function(lin){
      apply(bm,1,function(bm.i){
        any(bm.i %in% lin)
      })
    })
    pst <- sapply(linInds,function(ind){
      out <- pst.vec
      out[! ind] <- NA
      return(out)
    })
    rownames(pst) <- colnames(sim)
    colnames(pst) <- paste0('Lineage',1:ncol(pst))
    return(pst)
  }
  stop('Unrecognized num_paths argument.')
}

ordering_genes_pca <- function(cds, which.pcs = c(2,3)){
  exprs_filtered <- t(t(exprs(cds)/pData(cds)$Size_Factor))
  nz_genes <- which(exprs_filtered != 0)
  exprs_filtered[nz_genes] <- log(exprs_filtered[nz_genes] + 1)
  # Calculate the variance across genes without converting to a dense
  # matrix:
  expression_means <- Matrix::rowMeans(exprs_filtered)
  expression_vars <- Matrix::rowMeans((exprs_filtered - expression_means)^2)
  # Filter out genes that are constant across all cells:
  genes_to_keep <- expression_vars > 0
  exprs_filtered <- exprs_filtered[genes_to_keep,]
  expression_means <- expression_means[genes_to_keep]
  expression_vars <- expression_vars[genes_to_keep]
  # Here's how to take the top PCA loading genes
  pca_res <- prcomp(t(exprs_filtered), center=expression_means, scale.=sqrt(expression_vars))$rotation
  # Take the top 100 genes from components 2 and 3. Component
  # 1 usually is driven by technical noise.
  genes.list <- lapply(which.pcs, function(pc){
    names(sort(abs(pca_res[, pc]), decreasing = T))[1:100]
  })
  ordering_genes <- character()
  for(gl in genes.list){
    ordering_genes <- union(ordering_genes, gl)
  }
  return(ordering_genes)
}


### ACCURACY MEASURES

# s_{\pi_1\pi_2} from TSCAN paper
Spp <- function(pst1, pst2){
  pst1 <- as.numeric(pst1)
  pst2 <- as.numeric(pst2)
  # nicely formatted pst1, pst2
  spp <- function(p1,p2){
    A <- length(p1)
    cons <- 2/(A*(A-1))
    sigma <- sapply(seq_len(A),function(i){
      sapply(seq_len(A-i),function(k){
        j <- i+k
        if(any(is.na(c(p1[c(i,j)],p2[c(i,j)])))){
          return(0)
        }
        if((p1[j]-p1[i])*(p2[j]-p2[i]) >= 0){
          return(1)
        }
        return(0)
      })
    })
    return(cons*sum(unlist(sigma)))
  }
  if(is.null(names(pst1)) | is.null(names(pst2))){
    if(length(pst1)==length(pst2)){
      s <- spp(pst1,pst2)
    }else{
      stop('lengths must match or names must be provided')
    }
  }else{
    ns <- union(names(pst1),names(pst2))
    p1 <- pst1[match(ns,names(pst1))]
    p2 <- pst2[match(ns,names(pst2))]
    s <- spp(p1,p2)
  }
  return(s)
}

# modified Kendall's tau
kendall.mod <- function(pst1, pst2){
  # nicely formatted pst1, pst2
  spp <- function(p1,p2){
    A <- length(p1)
    d1 <- sign(outer(p1,p1,'-'))
    d1 <- d1[upper.tri(d1)]
    d2 <- sign(outer(p2,p2,'-'))
    d2 <- d2[upper.tri(d2)]
    concordance <- 2*(d1 == d2)-1
    concordance[is.na(concordance)] <- 0
    return(mean(concordance))
  }
  if(class(pst1)=='matrix'){
    n1 <- rownames(pst1)
    pst1 <- as.numeric(pst1)
    names(pst1) <- n1
  }
  if(class(pst2)=='matrix'){
    n2 <- rownames(pst2)
    pst2 <- as.numeric(pst2)
    names(pst2) <- n2
  }
  if(is.null(names(pst1)) | is.null(names(pst2))){
    if(length(pst1)==length(pst2)){
      keep <- !is.na(pst1) | !is.na(pst2)
      pst1 <- pst1[keep]
      pst2 <- pst2[keep]
      s <- spp(pst1,pst2)
    }else{
      stop('lengths must match or names must be provided')
    }
  }else{
    ns <- union(names(pst1),names(pst2))
    p1 <- pst1[match(ns,names(pst1))]
    p2 <- pst2[match(ns,names(pst2))]
    keep <- !is.na(p1) | !is.na(p2)
    p1 <- p1[keep]
    p2 <- p2[keep]
    s <- spp(p1,p2)
  }
  return(s)
}

# get agreement between a set of true lineages and inferred lineages,
# using the maximum modified Kendall's tau
agreement.pair <- function(truth, pst){
  truth <- as.matrix(truth)
  pst <- as.matrix(pst)
  
  kt.max <- apply(truth,2,function(truth.i){
    taus <- apply(pst,2,function(pst.j){
      names(truth.i) <- rownames(truth)
      names(pst.j) <- rownames(pst)
      return(kendall.mod(truth.i, pst.j))
    })
    return(max(taus,na.rm = TRUE))
  })
  return(kt.max)
}


### Other useful functions

# full-quantile normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank)
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}

# useful scaling function
scaleAB <- function(x,a=0,b=1){
    ((x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))*(b-a)+a
}

# pick best number of clusters based on silhouette width
pickBestK <- function(reducedDim, method = 'kmeans', reps = 10){
  require(cluster)
  if(method == 'kmeans'){
    clusterFun <- kmeans
  }
  if(method == 'pam'){
    clusterFun <- pam
  }
  asw <- sapply(2:15,function(k){
    asw.its <- sapply(seq_len(reps), function(i){
      cl <- clusterFun(reducedDim, k)
      if(is.null(cl$clustering)){
        cl$clustering <- cl$cluster
      }
      summary(silhouette(cl, dist = dist(reducedDim)))$avg.width
    })
    mean(asw.its)
  })
  return((2:15)[which.max(asw)])
}


# Extract an index denoting a single lineage in Monocle 2 output, corresponding
# to the given leaf state/cluster.
# Function based on code from buildBranchCellDataSet()
getLineageID <- function(cds, leaf_state){
  require(igraph)
  # checks
  if (is.null(pData(cds)$State) | is.null(pData(cds)$Pseudotime)) 
    stop("Please first order the cells in pseudotime using orderCells()")
  if (is.null(leaf_state)) 
    stop("Please specify the leaf_state to select subset of cells")
  # get mst
  if (cds@dim_reduce_type == "DDRTree") {
    pr_graph_cell_proj_mst <- minSpanningTree(cds)
  }else {
    pr_graph_cell_proj_mst <- cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree
  }
  # get root cell/state
  root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  if(length(root_cell)==0){
    root_cell <- rownames(pData(cds))[which.min(cds$Pseudotime)]
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
  }
  root_state <- pData(cds)[root_cell, ]$State
  pr_graph_root <- subset(pData(cds), State == root_state) # subset of pData corresponding to root cluster
  if (cds@dim_reduce_type == "DDRTree") {
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    # Y values are the principal points/graph
    root_cell_point_in_Y <- closest_vertex[row.names(pr_graph_root), ] # vertices in Y corresponding to root cluster 
  }else {
    root_cell_point_in_Y <- row.names(pr_graph_root)
  }
  root_cell <- names(which(degree(pr_graph_cell_proj_mst, v = root_cell_point_in_Y, 
                                  mode = "all") == 1, useNames = T))[1] # pick new root cell (a leaf)
  # above can sometimes return NA. As a backup: pick the original root cell (if its in the correct state)
  # or pick the cell from that state with the lowest pseudotime
  if(is.na(root_cell)){
    rc <- rownames(pData(cds))[which.min(cds$Pseudotime)]
    if(rc %in% names(degree(pr_graph_cell_proj_mst))){
      root_cell <- rc
    }else{
      root_cell <- names(degree(pr_graph_cell_proj_mst))[which.min(cds$Pseudotime[colnames(cds) %in% names(degree(pr_graph_cell_proj_mst))])]
    }
  }
  
  curr_cell <- subset(pData(cds), State == leaf_state)
  if (cds@dim_reduce_type == "DDRTree") {
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    curr_cell_point_in_Y <- closest_vertex[row.names(curr_cell), ] # vertices in Y corresopnding to leaf cluster
  }else {
    curr_cell_point_in_Y <- row.names(curr_cell)
  }
  curr_cell <- names(which(degree(pr_graph_cell_proj_mst, 
                                  v = curr_cell_point_in_Y, mode = "all") == 1, 
                           useNames = T))[1] # pick a leaf cell in leaf cluster
  path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, 
                                     curr_cell, root_cell)  # extract path from leaf cell to root cell
  path_to_ancestor <- names(unlist(path_to_ancestor$vpath)) # convert path to cell names
  # path_to_ancestor has extra cells for some reason (more than just those along path)
  
  if (cds@dim_reduce_type == "DDRTree") {
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    ancestor_cells_for_branch <- row.names(closest_vertex)[which(V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% 
                                                                   path_to_ancestor)]
    # ancestor_cells_for_branch has what we want
  }
  ancestor_cells_for_branch <- intersect(ancestor_cells_for_branch, 
                                         colnames(cds)) # sanity check?
  return(colnames(cds) %in% ancestor_cells_for_branch)
}


