ColNorm <- function(Dat){
  
  M = max(colSums(Dat))
  l <- length(colnames(Dat))
  
  for( i in 1:l){
    
    Dat[,i] = Dat[,i]*(M/sum(Dat[,i]))
    
  }
  
  return(Dat)
}


RowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}


Get.High.Var.Genes <- function(Dat, Bins = 10, frac = .2){
  
  M = rowMeans(Dat)
  V = RowVar(Dat)
  
  BinLims <- seq(log(10/length(M)),max(M),Bins + 1)
  
  In = c()
  
  for (i in 1:Bins){
    In_red <- which(as.logical((M >= BinLims[i])*
                                 (M < BinLims[i+1])))
    temp <- V[In_red]/M[In_red]
    temp2 <- sort(temp, decreasing = T, index.return = T)
    I <- temp2$ix
    topX <- round(length(In_red)*frac)
    In <- c(In,In_red[I[1:topX]])
    
  }
  
  In <- In[!is.na(In)]
  In <- unique(In)
  return(In)
  
}


Rescale.Rows <- function(Dat){
  
  #s <- dim(Dat)
  
  #for (i in 1:s[1]){
  
  #  Dat[i,] <- (Dat[i,] - min(Dat[i,]))/max(Dat[i,] - min(Dat[i,])) 
  
  
  #}
  
  Dat <- t(apply(Dat, 1, function(x)(x-min(x))/(max(x)-min(x))))
  
  return(Dat)
}

Get.Module.Averages <- function(Dat, GeneNames, ModList, RescaleRows = T){
  
  library(Matrix)
  
  if (RescaleRows == T){
    Dat <- Rescale.Rows(Dat)
  }
  
  d <- dim(Dat)
  print(d)
  
  ModNames <- unique(ModList$Module)
  M <- matrix(rep(0,d[2]*length(ModNames)), nrow = length(ModNames), ncol = d[2])
  
  for( i in 1:length(ModNames)){
    In <- which(ModList$Module == ModNames[i])
    In2 <- which(GeneNames %in% ModList$GeneID[In])
    M[i,] <- colMeans(Dat[In2,])
    
  }
  
  colnames(M) <- colnames(Dat)
  rownames(M) <- ModNames
  
  M <- as.data.frame(M)
  
  return(M)
  
}

