

RunSlicerTraj <- function(Dat, m, kmin = 5, kmax = 50, by = 5, reg = 2){
  
  # m is instrinsic dimension of data
  # kmin is minimum number of nearest neighbors
  
  library(SLICER)
  
  Dat <- t(Dat)
  
  k = select_k(Dat, kmin, kmax, by)
  print(k)
  traj_lle = lle(Dat, m, k, reg = reg)$Y
  traj_graph = conn_knn_graph(traj_lle,5)
  ends = find_extreme_cells(traj_graph, traj_lle)
  
  l <- list()
  l$traj_lle <- traj_lle
  l$traj_graph <- traj_graph
  l$traj_graph2 <- as.matrix(as_adj(traj_graph))
  l$ends <- ends 
  
  return(l)
  
}

RunSlicerTrajK <- function(Dat, m, k = 5, reg = 2){
  
  # m is instrinsic dimension of data
  # kmin is minimum number of nearest neighbors
  
  library(SLICER)
  
  Dat <- t(Dat)
  
  traj_lle = lle(Dat, m, k, reg = reg)$Y
  
  if (m == 1){
    l <- length(traj_lle)
    traj_lle = rep(traj_lle,2)
    dim(traj_lle) <- c(l,2)
  }
  
  
  
  traj_graph = conn_knn_graph(traj_lle,5)
  ends = find_extreme_cells(traj_graph, traj_lle)
  
  l <- list()
  l$traj_lle <- traj_lle
  l$traj_graph <- traj_graph
  l$traj_graph2 <- as.matrix(as_adj(traj_graph))
  l$ends <- ends 
  
  return(l)
  
}

SlicerOrderCells <- function(traj_graph, start){
  
  #Given a starting cell, orders the other cells
  library(SLICER)
  
  cells_ordered = cell_order(traj_graph, start)
  branches = assign_branches(traj_graph,start)
  
  l <- list()
  
  l$cells_ordered <- cells_ordered
  l$branches <- branches 
  
  return(l)
  
}


WriteToMat <- function(l, fileName){
  
  library(R.matlab)
  writeMat(fileName, l = l)
  
  return(NULL)
}


RunMonocleDefault <- function(Dat, Labels, max_components = 2, meth = 'DDRTree'){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  Genes <- c(1:dim(Dat)[1])
  Genes <- data.frame(Genes)
  Labels <- data.frame(Labels)
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = Genes)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=negbinomial.size())
  
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  
  #HSMM <- setOrderingFilter(HSMM, c(1:length(Genes)))
  
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  HSMM <- orderCells(HSMM)
  plot_cell_trajectory(HSMM, color_by="Labels")
  
  return(HSMM)
  
}


RunMonocleTobit <- function(Dat, Labels, max_components=2, meth = 'DDRTree',C_by = NULL, 
                            gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  colnames(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  HSMM <- orderCells(HSMM)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}

Monocle.Plot <- function(MonRun, Labels, Discrete = T, 
                         titleVal = 'Lineage Scatter Plot'){
  
  library(ggplot2)
  library(reshape2)
  
  S <- MonRun@reducedDimS
  K <- MonRun@reducedDimK
  
  if(Discrete == T){
    Labels = factor(Labels)
  }
  
  df1 <- data.frame('x' = S[1,], 'y' = S[2,],
                    'Lab' = Labels)
  
  Pstime <- MonRun@phenoData@data$Pseudotime
  temp <- sort(Pstime, decreasing = F, index.return = T)
  
  
  #df2 <- data.frame('x2' = K[1,temp$ix], 'y2' = K[2,temp$ix])
  
  
  
  if(Discrete==T){
    g <- ggplot(df1, aes(x = x,y = y, color = Lab)) + geom_point(size=3)  + ggtitle(titleVal)
    #    geom_point(data = na.omit(dat), aes(x = x2,y = y2))  
  } else{
    g <- ggplot(df1, aes(x = x, y = y, color = Lab)) + 
      geom_point(size=3) + 
      scale_color_gradient(low = "black", high = "red") + ggtitle(titleVal)
    #      geom_line(data = df2, aes(x = x,y = y))  
  }
  
  return(g)
  
}


Mon.Plot.Genes <- function(MonRun, Dat, GeneNames, GeneID){
  
  In <- which(GeneNames == GeneID)
  Lab = as.numeric(Dat[In,])
  g <- Monocle.Plot(MonRun, Labels = Lab, Discrete = F, titleVal = GeneID)
  return(g)
  
}


Run.Lineage.BR <- function(Dat4, DatNorm4, FileName, BR, Sex, l = NULL){
  
  temp <- DatNorm4
  rownames(temp) <- NULL
  colnames(temp) <- NULL
  
  
  if (is.null(l)){
    MonRun <- RunMonocleTobit(temp, Dat4$AgeAtDeath)
    
    g <- Monocle.Plot(MonRun, Labels = Dat4$Tissue.Diagnosis, Discrete = T)
    t <- paste(FileName,BR,Sex,'Diag.png',sep = '')
    ggsave(filename=t, plot=g)
    
    
    g <- Monocle.Plot(MonRun, Labels = Dat4$Tissue.APOE4, Discrete = T)
    t <- paste(FileName,BR,Sex,'ApoE.png',sep = '')
    ggsave(filename=t, plot=g)
    
    
    g <- Monocle.Plot(MonRun, Labels = Dat4$AgeAtDeath, Discrete = F)
    t <- paste(FileName,BR,Sex,'Age.png',sep = '')
    ggsave(filename=t, plot=g)
    
  } else { 
    
    MonRun <- RunMonocleTobit(temp, Dat4[[l[3]]])
    
    g <- Monocle.Plot(MonRun, Labels = Dat4[[l[1]]], Discrete = T)
    t <- paste(FileName,BR,Sex,'Diag.png',sep = '')
    ggsave(filename=t, plot=g)
    
    
    g <- Monocle.Plot(MonRun, Labels = Dat4[[l[2]]], Discrete = T)
    t <- paste(FileName,BR,Sex,'ApoE.png',sep = '')
    ggsave(filename=t, plot=g)
    
    
    g <- Monocle.Plot(MonRun, Labels = Dat4[[l[3]]], Discrete = F)
    t <- paste(FileName,BR,Sex,'Age.png',sep = '')
    ggsave(filename=t, plot=g)
    
    
  }
  
  
  
}


Subset.Data.Sex.BR <- function(DatNorm2, Dat2, BR,Sex,Tdiag = T){
  
  if (Tdiag == T){
    In_BR <- grep(BR,Dat2$Tissue.Diagnosis)
    DatNorm3 <- DatNorm2[,In_BR]
    Dat3 <- Dat2[In_BR,]
  } else {
    DatNorm3 <- DatNorm2
    Dat3 <- Dat2
  }
  
  In_S <- which(Dat3$Sex == Sex)
  print(In_S)
  
  DatNorm4 <- DatNorm3[,In_S]
  Dat4 <- Dat3[In_S,]
  
  l <- list()
  l$Dat4 <- Dat4
  l$DatNorm4 <- DatNorm4
  
  return(l)
}


Normalize.Data <- function(Dat, Dat2, AMP_mods, DelChars = T){
  
  GeneNames <- Dat$ensembl_gene_id
  GeneNamesAD <- AMP_mods$GeneID
  
  Names <- colnames(Dat)
  
  if (DelChars ==T){
    for (i in 1:length(Names)){
      
      Names[i] <- substring(Names[i],2)
      
    }
    
  }
  
  
  colnames(Dat) <- Names
  cNames <- Dat2$SampleID
  l <- length(Names)
  
  #deleting columns not in the covariate list
  temp <- rep(T,l)
  for (i in 1:l){
    if (!(Names[i] %in% cNames)){
      temp[i] <- F
    }
  }
  
  In <- which(temp)
  #print(temp)
  Dat <- Dat[,In]
  
  #deleting extra rows in covariate list
  Names <- Names[In]
  l <- length(cNames)
  temp <- rep(T,l)
  for (i in 1:l){
    if (!(cNames[i] %in% Names)){
      temp[i] <- F
    }
  }
  In <- which(temp)
  Dat2 <- Dat2[In,]
  
  
  #Normalize all columns 
  source('MiscPreprocessing.R')
  
  DatNorm <- ColNorm(Dat)
  In <- which(GeneNames %in% GeneNamesAD)
  DatNorm2 <- DatNorm[In,]
  
  l <- list()
  l$Dat2 <- Dat2
  l$DatNorm <- DatNorm 
  l$DatNorm2 <- DatNorm2
  
  return(l)
  
}

Sliding.Window.Average.Cat <- function(PStime, X, wSize, t){
  
  l <- length(PStime)
  nB <- l-wSize + 1 
  pAvg <- rep(0,nB)
  fAvg <- rep(0,nB)
  
  In_P <- sort(PStime,index.return = T)
  In_P <- In_P$ix 
  PStime <- PStime[In_P]
  X <- X[In_P]
  #print(X)
  
  
  for (i in 1:nB){
    t2 <- i+wSize-1
    
    pAvg[i] <- mean(PStime[i:t2])
    fAvg[i] <- sum(X[i:t2])/wSize
    #print(fAvg[i])
    #print(sum(X[i:i+wSize-1]==t))
  }
  
  l2 <- list()
  l2$pAvg <- pAvg
  l2$fAvg <- fAvg
  l2$PStime <- PStime
  l2$X <- X
  
  return(l2)
}


Make.Gene.Symb <- function(GeneENSG){
  
  source('convertEnsemblToHgnc.R')
  GeneConv <- convertEnsemblToHgnc(GeneENSG)
  Symb <- as.character(c(1:length(GeneENSG)))
  
  for (i in 1:length(GeneENSG)){
    In <- which(GeneConv$ensembl_gene_id == GeneENSG[i])
    if (length(In)>0){
      Symb[i] <- GeneConv$external_gene_name[In]
    }
  }
  
  return(Symb)
  
}


ScaleMatrix <- function(Dat){
  
  for(i in 1:dim(Dat)[1]){
    Dat[i,] <- Dat[i,] - min(Dat[i,])
    Dat[i,] <- Dat[i,]/max(Dat[i,])
  }
  
  return(Dat)
  
}


