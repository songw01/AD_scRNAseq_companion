load.data <- function (filename, gsub.from = "\\.", gsub.to = "-") 
{
    df <- read.delim(file = filename, sep = "\t", header = T)
    out <- as.matrix(df[, 2:ncol(df)])
    rownames(out) <- as.character(df[[1]])
    if (!is.null(gsub.from)) 
        colnames(out) <- gsub(gsub.from, gsub.to, colnames(out))
    return(out)
}
############################# FET pipeline
read.geneSet <- function(geneSet.file)
{
 gene.list <- readLines(geneSet.file)
 gene.list <- strsplit(gene.list,"\t")
 
 names(gene.list) <- sapply(gene.list,function(x) x[1])
 gene.list <- lapply(gene.list,function(x) x[2:length(x)])
 return(gene.list)
}

output.geneSet.file <- function(geneSet,outputfname)
{
 if (!is.list(geneSet)) stop("geneSet is not a list.")
 if (is.null(names(geneSet))) stop("names(geneSet) is not defined properly.")
  
 sink(outputfname)
 cat(paste(paste(names(geneSet),"\t",sapply(geneSet,function(x) paste(x,collapse = "\t")),sep = ""),collapse = "\n"))
 sink()

 return(0)
}

make.Pairwise.Tables <- function(geneSets1,geneSets2,background)
{
 mem1 <- do.call(cbind,lapply(geneSets1,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))
 mem2 <- do.call(cbind,lapply(geneSets2,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))

 d <- t(mem1) %*% mem2;
 b <- abs(t(mem1) %*% (mem2-1))
 c <- abs(t(mem1-1) %*% (mem2))
 a <- t(mem1-1) %*% (mem2-1);

 ij <- do.call(rbind,lapply(1:length(geneSets1),function(i,j) cbind(rep(i,length(j)),j),j = 1:length(geneSets2)))

 pairwise.tables <- lapply(1:nrow(ij),function(i,ij,a,b,c,d) as.table(matrix(c(a[ij[i,1],ij[i,2]],b[ij[i,1],ij[i,2]],c[ij[i,1],ij[i,2]],d[ij[i,1],ij[i,2]]),nrow = 2)),ij = ij,a = a,b = b,c = c,d = d)
 names(pairwise.tables) <- paste(names(geneSets1)[ij[,1]],"~vs~",names(geneSets2)[ij[,2]],sep = "")
 return(pairwise.tables)
}

do.FisherExactTest <- function(table.count,or = 1,alternative = "greater",N_bg = NULL)
{
 if (is.null(N_bg)) N_bg = sum(rowSums(table.count))
 
 out <- fisher.test(x = table.count,or = or,alternative = alternative)
 odds.ratio <- out$estimate
 p.value <- out$p.value;
 geneSet1.count <- rowSums(table.count)[2]
 geneSet2.count <- colSums(table.count)[2]
 expected.count <- geneSet1.count/N_bg * geneSet2.count
 overlap.count <- table.count[2,2];
 fold.change <- overlap.count/expected.count
 
 out <- c(N_bg,geneSet1.count,geneSet2.count,expected.count,overlap.count,fold.change,odds.ratio,p.value)
 names(out) <- c("Background","set1_size","set2_size","expected.overlap","actual.overlap","enrichment.foldchange","odds.ratio","FET_pvalue")
 return(out)
}

require(parallel)
require(doParallel)
require(foreach)
require(iterators)
perform.AllPairs.FET <- function(geneSets1,geneSets2,background,
or = 1,alternative = "greater",
adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
{
 pairwise.tables <- make.Pairwise.Tables(geneSets1,geneSets2,background)
 
 if (do.multicore)
 {
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
  fact <- do.call(c,lapply(1:n.cores,function(i,dn) rep(i,dn),dn = ceiling(length(pairwise.tables)/n.cores)))
  fact <- factor(fact[1:length(pairwise.tables)])
  split.tables <- lapply(split(1:length(pairwise.tables),fact),function(i,obj) obj[i],obj = pairwise.tables);rm(fact)
  output <- foreach (tbl = split.tables,.combine = 'c',.export = c("do.FisherExactTest")) %dopar% {
            out <- lapply(tbl,function(x,y,z) do.FisherExactTest(table.count = x,or = y,alternative = z),y = or,z = alternative) 
			return(out)
  }
 
 }else{
  output <- lapply(pairwise.tables,function(x,y,z) do.FisherExactTest(table.count = x,or = y,alternative = z),y = or,z = alternative)
 }
 # collect outputs
 output <- data.frame(set1_Name = gsub("~vs~(.*)$","",names(output)),set2_Name = gsub("^(.*)~vs~","",names(output)),as.data.frame(do.call(rbind,output)))
 int.element <- mapply(FUN = function(x,y) paste(intersect(x,y),collapse = ","),x = geneSets1[as.character(output[[1]])],y = geneSets2[as.character(output[[2]])],SIMPLIFY = F)
 int.element <- do.call(c,int.element)
 if (adjust.FET.pvalue) output <- cbind.data.frame(output,data.frame(corrected.FET.pvalue = p.adjust(output$FET_pvalue,"bonferroni"),intersects = int.element))
 
 return(output)
}

make.paired.Tables <- function(geneSets1,geneSets2,ij,background)
{
 
 #mem1 <- do.call(cbind,lapply(geneSets1,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))
 #mem2 <- do.call(cbind,lapply(geneSets2,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))

 pairwise.tables <- mapply(FUN = function(s1,s2,z) {
                                         v1 <- rep(0,length(z));v1[which(z %in% s1)] <- 1; 
										 v2 <- rep(0,length(z));v2[which(z %in% s2)] <- 1;
										 # n11,n12,n21,n22
										 as.table(matrix(c(sum(abs(v1 - 1) * abs(v2 - 1)),sum(abs(v2 - 1) * v1),sum(abs(v1 - 1) * v2),sum(v1 * v2)),nrow = 2))
                                },s1 = geneSets1[ij[,1]],s2 = geneSets2[ij[,2]],MoreArgs = list(z = background),SIMPLIFY = FALSE)
 names(pairwise.tables) <- paste(names(geneSets1)[ij[,1]],"~vs~",names(geneSets2)[ij[,2]],sep = "")
 return(pairwise.tables)
}


perform.ijPairs.FET <- function(geneSets1,geneSets2,ij,background,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
{
 require(doParallel)
 pairwise.tables <- make.paired.Tables(geneSets1,geneSets2,ij,background)
 
 if (do.multicore & getDoParWorkers() == 1)
 {
  export.func <- c("do.FisherExactTest","make.Pairwise.Tables","make.paired.Tables")
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  cat(paste("registered:",getDoParWorkers()," cores\n",sep = ""))
  fact <- do.call(c,lapply(1:n.cores,function(i,dn) rep(i,dn),dn = ceiling(length(pairwise.tables)/n.cores)))
  fact <- factor(fact[1:length(pairwise.tables)])
  split.tables <- lapply(split(1:length(pairwise.tables),fact),function(i,obj) obj[i],obj = pairwise.tables);rm(fact)
  cat("commence parallelized enrichment analyses\n")
  output <- foreach (tbl = split.tables,.combine = 'c',.export = export.func) %dopar% {
            out <- lapply(tbl,do.FisherExactTest) 
			return(out)
  }
  
  stopCluster(cl)
 
 }else{
  output <- lapply(pairwise.tables,do.FisherExactTest)
 }
 # collect outputs
 output <- data.frame(set1_Name = gsub("~vs~(.*)$","",names(output)),set2_Name = gsub("^(.*)~vs~","",names(output)),as.data.frame(do.call(rbind,output)))
 int.element <- mapply(FUN = function(x,y) paste(intersect(x,y),collapse = ","),x = geneSets1[as.character(output[[1]])],y = geneSets2[as.character(output[[2]])],SIMPLIFY = F)
 int.element <- do.call(c,int.element)
 if (adjust.FET.pvalue) output <- cbind.data.frame(output,data.frame(corrected.FET.pvalue = p.adjust(output$FET_pvalue,"bonferroni"),intersects = int.element))
 
 output <- output[order(output$FET_pvalue),]
 return(output)
}

############################## DEG pipeline

extract.DEG <- function(DEG.table,id.col,pval.col,pval.cutoff = 0.5,fc.col,fc.cutoff = 1.2,is.log = TRUE)
{
 if (is.character(pval.col)) pval.col <- which(colnames(DEG.table) == pval.col)
 if (is.character(fc.col)) fc.col <- which(colnames(DEG.table) == fc.col)
 if (is.character(id.col)) id.col <- which(colnames(DEG.table) == id.col)
 
 if (is.log) 
 {
  up.fc <- log2(fc.cutoff);
  dn.fc <- log2(1/fc.cutoff)
 }else{
  up.fc <- fc.cutoff;
  dn.fc <- 1/fc.cutoff
 }
 
 up.i <- which(DEG.table[[pval.col]] < pval.cutoff & DEG.table[[fc.col]] > up.fc)
 dn.i <- which(DEG.table[[pval.col]] < pval.cutoff & DEG.table[[fc.col]] < dn.fc)
 
 output <- list(UP = unique(as.character(DEG.table[[id.col]][up.i])),DN = unique(as.character(DEG.table[[id.col]][dn.i])))
 output <- c(output,list(ALL = Reduce("union",output)))
 
 names(output) <- paste(colnames(DEG.table)[pval.col],"_",pval.cutoff,"__fc_",fc.cutoff,"__",names(output),sep = "")
 
 return(output)
}


load.modules <- function(module.file,file.type)
{
 if (file.type == "MEGENA")
 {
  output <- read.geneSet(module.file)
 }
 
 if (file.type == "WGCNA")
 {
  module.table <- read.delim(file = module.file,sep = "\t",header = T)
  output <- lapply(split(as.character(module.table$Gene.symbol),factor(module.table$module)),unique)
 }
 
 return(output)
}

########## summary functions for FET

extract.SigTerms <- function(FET.file,pvalue.cutoff = 0.05,fold.cutoff = 2,pval.colname,fold.colname = "enrichment.foldchange")
{
 FET.table <- read.delim(file = FET.file,sep = "\t",header = T)
 pvalue.vec <- FET.table[[which(names(FET.table) == pval.colname)]]
 fold.vec <- FET.table[[which(names(FET.table) == fold.colname)]]
 
 ii <- which(pvalue.vec < pvalue.cutoff & fold.vec > fold.cutoff)
 
 sig.module <- factor(FET.table$set1_Name[ii])
 sig.term <- as.character(FET.table$set2_Name[ii])
 sig.pval <- pvalue.vec[ii]
 sig.fold <- fold.vec[ii]
 
 ## order by pvalues
 ii <- order(sig.pval);
 sig.module <- sig.module[ii];sig.term <- sig.term[ii];sig.pval <- sig.pval[ii];sig.fold <- sig.fold[ii]
 
 sig.str <- paste(sig.term,signif(sig.pval,4),signif(sig.fold,4),sep = "/");
 
 split.str <- split(sig.str,sig.module)
 
 output <- data.frame(Module = names(split.str),Significant.Term = sapply(split.str,function(x) paste(x,collapse = ",")))
 return(output)
}

format.FET.output <- function(FET.folder,pvalue.cutoff = 0.05,fold.cutoff = 2)
{
 cat(paste("##### formatting FET outputs from folder:",FET.folder,"\n",sep = ""))
 FET.table.files <- list.files(path = FET.folder,pattern = "_FET-Table.txt$",full.names = T);
 cat("Processing the outputs from the following files:\n");
 cat(paste(paste(FET.table.files,collapse = "\n"),"\n",sep = ""));
 
 names(FET.table.files) <- gsub("_FET-Table.txt","",gsub("(.*)/","",FET.table.files))
 
 SigTerms <- lapply(FET.table.files,extract.SigTerms,pval.colname = "FET_pvalue",pvalue.cutoff = pvalue.cutoff,fold.cutoff = fold.cutoff)
 AllPair.SigTerms <- lapply(FET.table.files,extract.SigTerms,pval.colname = "corrected.FET.pvalue",pvalue.cutoff = pvalue.cutoff,fold.cutoff = fold.cutoff)
 
 module.ticks <- sort(union(Reduce("union",lapply(SigTerms,function(x) as.character(x[[1]]))),Reduce("union",lapply(AllPair.SigTerms,function(x) as.character(x[[1]])))))
 
 output.matrix <- matrix(NA,nrow = length(module.ticks),ncol = length(FET.table.files)*2)
 
 sig.i <- seq(1,ncol(output.matrix)-1,2)
 allpair.i <- seq(2,ncol(output.matrix),2)
 
 for (i in 1:length(FET.table.files))
 {
  ii <- match(as.character(SigTerms[[i]][[1]]),module.ticks)
  if (length(ii) > 0) output.matrix[ii,sig.i[i]] <- as.character(SigTerms[[i]][[2]])
  rm(ii)

  ii <- match(as.character(AllPair.SigTerms[[i]][[1]]),module.ticks);
  if (length(ii) > 0) output.matrix[ii,allpair.i[i]] <- as.character(AllPair.SigTerms[[i]][[2]])
  rm(ii)
 }
 
 colname.vec <- rep("",ncol(output.matrix))
 colname.vec[sig.i] <- paste("SigTerms.",names(FET.table.files),sep = "")
 colname.vec[allpair.i] <- paste("AllPair.SigTerms.",names(FET.table.files),sep = "")
 
 colnames(output.matrix) <- colname.vec
 
 output.matrix <- data.frame(Module = module.ticks,as.data.frame(output.matrix))
 
 return(output.matrix)
}

############# 1st PC analysis
module.1stPC <- function(x,m)
{
 sub.m <- na.omit(rbind(m[x,]));
 colnames(sub.m) <- colnames(m)
 if (nrow(sub.m) > 1)
 {
  PCres = tryCatch(prcomp(x = sub.m,scale = T,center = T),error = function(e) return(NULL))
  out <- PCres$rotation[,1]
 }else{
  if (nrow(sub.m) == 1) {out <- sub.m}else{out <- NULL}
 }
 return(out)
}

trait.assoc <- function(trait.vec,pc.matrix)
{
 data.frame(module = rownames(pc.matrix),as.data.frame(t(apply(pc.matrix,1,function(x,y) do.call('c',cor.test(x,y,method = "spearman",use = "pairwise.complete.obs")[c("estimate","p.value")]),y = trait.vec))))
}

all.trait.assoc <- function(trait.df,pc.matrix)
{
 pc.res <- lapply(trait.df,function(x,y) trait.assoc(x,y),y = pc.matrix)
 rho.mat <- do.call('cbind',lapply(pc.res,function(x) x[[2]]));colnames(rho.mat) <- paste("spearman.rho",colnames(trait.df),sep = "__")
 p.mat <- do.call('cbind',lapply(pc.res,function(x) x[[3]]));colnames(p.mat) <- paste("spearman.rho.P",colnames(trait.df),sep = "__")
 data.frame(module = rownames(pc.matrix),as.data.frame(rho.mat),as.data.frame(p.mat))
}

##### update with permutation base.

count.FD <- function(rho,thresh.vec)
{
	 # combine threshold and correlation values for efficiency
	 n.rho <- length(rho)
	 label.vec <- c(rep(1,length(rho)),rep(0,length(thresh.vec)));# label correlation value and threshold value
	 rho <- c(rho,thresh.vec)# combine correlation and threhold values for efficiency.
	 # sort correlation values
	 i <- order(rho)
	 rho <- rho[i]
	 label.vec <- label.vec[i]
	 
	 # count permuted correlation values above thresholds
	 j <- which(label.vec == 0)
	 output <- cbind(rho[j],(n.rho - cumsum(label.vec)[j])/n.rho)
	 colnames(output) <- c("rho.cutoff","FDR")
	 return(output)
} 

trait.assoc.perm <- function(trait.vec,pc.matrix,do.permute = TRUE,n.perm = 100)
{
 output.df <- data.frame(module = rownames(pc.matrix),as.data.frame(t(apply(pc.matrix,1,function(x,y) do.call('c',cor.test(x,y,method = "spearman",use = "pairwise.complete.obs")[c("estimate","p.value")]),y = trait.vec))))
 
 # use correlation coefficient to calculate FDR
 if (do.permute)
 {
  
  thresh.vec <- abs(output.df$estimate.rho);
  names(thresh.vec) <- as.character(output.df[[1]])
  thresh.vec <- sort(thresh.vec)
  cat("-- commence permuting:");cat(n.perm);cat(" times\n");
  perm.ind <- lapply(1:n.perm,function(i,n) sample(1:n,n),n = ncol(pc.matrix))
  perm.rho <- lapply(perm.ind,function(i,tvec,mat) cor(tvec,t(mat[,i]),method = "spearman",use = "pairwise.complete.obs"),tvec = trait.vec,mat = pc.matrix) 
  cat("summarizing stats...\n")
  z <- rank(-thresh.vec);
  PR <- z/length(thresh.vec);
  count.out <- lapply(perm.rho,function(x,vec) count.FD(rho = x,thresh.vec = vec),vec = thresh.vec)
  FPR = Reduce("+",lapply(count.out,function(x) x[,2]))/n.perm; 
  FDR = FPR/PR;FDR[which(FPR == 0)] <- 0;FDR[which(FDR > 1)] <- 1;
  FDR.table <- data.frame(threshold = thresh.vec,FPR = FPR,PR = PR,FDR = FDR)
  
  FDR.table <- FDR.table[match(as.character(output.df[[1]]),rownames(FDR.table)),]
  output.df <- data.frame(output.df,BH.FDR = p.adjust(output.df$p.value,"BH"),perm.FDR = FDR.table$FDR)
 }
 
 return(output.df)
}

all.trait.permAssoc <- function(trait.df,pc.matrix)
{
 pc.res <- lapply(trait.df,function(x,y) trait.assoc.perm(x,y),y = pc.matrix)
 for (i in 1:length(pc.res)) colnames(pc.res[[i]])[-1] <- paste(colnames(trait.df)[i],colnames(pc.res[[i]])[-1],sep = "__");
 names(pc.res) <- NULL
 data.frame(module = rownames(pc.matrix),do.call('cbind.data.frame',lapply(pc.res,function(x) x[,-1])))
}
read.MSigDB <- function(MSigDB.file)
{
 gene.list <- readLines(MSigDB.file)
 gene.list <- strsplit(gene.list,"\t")
 
 names(gene.list) <- sapply(gene.list,function(x) x[1])
 gene.list <- lapply(gene.list,function(x) x[3:length(x)])
 return(gene.list)
}
perform.FET.output <- function(geneSets1,geneSets2,background,outputfname,adjust.FET.pvalue = T,pvalue.cutoff = 0.05,do.multicore = F,n.cores = NULL)
{
 if (is.null(names(geneSets1))) 
 {
  cat("#### No names were tagged to geneSets1 #####\nProviding names to geneSets1...\n")
  names(geneSets1) <- paste("A",1:length(geneSets1),sep = "")
 }
 
 if (is.null(names(geneSets2))) 
 {
  cat("#### No names were tagged to geneSets2 #####\nProviding names to geneSets1...\n")
  names(geneSets2) <- paste("B",1:length(geneSets2),sep = "")
 }
 cat("###### Performing Fisher Exact Test for over-representation\n")
 cat(paste("- Number of geneSets1:",length(geneSets1),"\n","- Number of geneSets2:",length(geneSets2),"\n",sep = ""))
 
 FET.table <- perform.AllPairs.FET(geneSets1 = geneSets1,geneSets2 = geneSets2,background = background,adjust.FET.pvalue = adjust.FET.pvalue,do.multicore = do.multicore,n.cores = n.cores)
 
 # output gene sets
 cat("###### Outputting files...\n- Output gene sets\n")
 geneSet.files <- paste(outputfname,"_geneSets",c(1,2),".txt",sep = "")
 cat(paste(geneSet.files,"\n",sep = ""))
 output.status <- output.geneSet.file(geneSet = geneSets1,outputfname = geneSet.files[1])
 output.status <- output.geneSet.file(geneSet = geneSets2,outputfname = geneSet.files[2])
 cat("- Output FET output table\n")
 write.table(FET.table,file = paste(outputfname,"_FET-Table.txt",sep = ""),sep = "\t",row.names = F,col.names = T,quote = F)
 cat("- Output significant results\n")
 sig.table <- FET.table[FET.table$corrected.FET.pvalue < pvalue.cutoff,];
 write.table(sig.table,file = paste(outputfname,"_SignificantFET-Table.txt",sep = ""),sep = "\t",row.names = F,col.names = T,quote = F)
 
 return(0) 
}



run.MSigDB <- function(module.geneSets,GeneSymbol,MSigDB.files,saveto,
annot.table = NULL,id.col = NULL,symbol.col = NULL,
do.multicore = F,n.cores = NULL)
{
 # replaces run.MSigDB.Enrichment() function 
 
 dir.create(saveto)
 
 if (is.null(MSigDB.files))
 {
  data(MSigDB.geneSets)
 }else{
  names(MSigDB.files) <- gsub("\\.gmt$","",gsub("(.*)/","",MSigDB.files))
  MSigDB.geneSets <- lapply(MSigDB.files,read.MSigDB)
 }
 
 if (is.null(GeneSymbol)) GeneSymbol <- Reduce("union",module.geneSets);
 
 if (!is.null(annot.table))
 {
  # if annotation is provided, project MSigDB gene symbols to array
  gene.vec <- as.character(annot.table[[symbol.col]]);
  names(gene.vec) <- as.character(annot.table[[id.col]])
  for (i in 1:length(MSigDB.geneSets))
  {
   MSigDB.geneSets[[i]] <- lapply(MSigDB.geneSets[[i]],function(x,y) {out <- names(y)[match(x,y)];out <- out[!is.na(out)];return(out)},y = gene.vec)
   MSigDB.geneSets[[i]] <- MSigDB.geneSets[[i]][sapply(MSigDB.geneSets[[i]],length) > 0]
  }
 }
 
 #### load up MSigDB gene Sets
 
 for (i in 1:length(MSigDB.geneSets))
 {
  output.status <- perform.FET.output(geneSets1 = module.geneSets,geneSets2 = MSigDB.geneSets[[i]],background = GeneSymbol,outputfname = paste(saveto,names(MSigDB.geneSets)[i],sep = "/"),do.multicore = do.multicore,n.cores = n.cores)
 }
 
 return(0)

}

get.FET.summary <- function(FET.table,id.col = "id",module.col = "module",pval.col = "corrected.FET.pvalue",
EFC.col = "enrichment.foldchange",pval.cutoff = 0.05)
{
	id.col <- which(colnames(FET.table) == id.col)
	module.col <- which(colnames(FET.table) == module.col)
	pval.col <- which(colnames(FET.table) == pval.col)
	EFC.col <- which(colnames(FET.table) == EFC.col)

	sig.row <- which(FET.table[[pval.col]] < pval.cutoff)
	
	output <- NULL

	if (length(sig.row) > 0)
	{
	 summary.vec <- paste(as.character(FET.table[[id.col]][sig.row]),format(signif(FET.table[[pval.col]][sig.row],3),scientific = TRUE),
	 format(signif(FET.table[[EFC.col]][sig.row],3),scientific = TRUE),sep = "/")
	 output <- sapply(split(summary.vec,factor(FET.table[[module.col]][sig.row])),function(x) paste(x,collapse = ","))
	}

	return(output)
}

get.FET.matrix <- function(FET.table,mod.names,sig.names,module.col = "module",id.col = "id")
{
	pval.matrix <- matrix(1,nrow = length(mod.names),ncol = length(sig.names));
	rownames(pval.matrix) <- mod.names;
	colnames(pval.matrix) <- paste("FET.P__",sig.names,sep = "");
	fc.matrix <- matrix(0,nrow = length(mod.names),ncol = length(sig.names));
	rownames(fc.matrix) <- mod.names;
	colnames(fc.matrix) <- paste("EFC__",sig.names,sep = "");

	row.ij <- cbind(match(FET.table[[which(colnames(FET.table) == module.col)]],mod.names),match(FET.table[[which(colnames(FET.table) == id.col)]],sig.names))
	for (n in 1:nrow(row.ij)) 
	{
		pval.matrix[row.ij[n,1],row.ij[n,2]] <- FET.table$FET_pvalue[n];
		fc.matrix[row.ij[n,1],row.ij[n,2]] <- FET.table$enrichment.foldchange[n];
	}
	
	output <- list(p.value = pval.matrix,EFC = fc.matrix)
	return(output)
}

#### ranking

module_rank <- function (X,thresh = 0.05)
{
    X = as.matrix(X)
    G = apply(X, 2, function(y,th) {
        z = rank(y)
        M = max(z)
		z[which(y > th)] <- M
        (M + 1 - z)/M
    },th = thresh)
    apply(G, 1, function(x) -sum(log(x)))
}


module_rank.v1 <- function (X)
{
    X = as.matrix(X)
    if (any(X) == 0)
    {
	 zero.i <- which(X == 0);
	 X[zero.i] <- 1e-320;
	 X[X == 1] <- 1 - 1e-320
    }
    apply(X, 1, function(x)-mean(-log10(x),na.rm = T))
}

module_rank.v2 <- function (X)
{
    X = as.matrix(X)
    if (any(X) == 0)
    {
	 zero.i <- which(X == 0);
	 X[zero.i] <- 1e-320;
	 X[X == 1] <- 1 - 1e-320
    }
    apply(X, 1, function(x) sum(log10(-log10(x),na.rm = T)))
}



generate.rankTable <- function(trait.file,sig.file,module.file,summary.file,exclude.var)
{
	modules <- load.modules(module.file = module.file,file.type = "MEGENA")
	module.ids <- names(modules);rm(modules)

	# process trait + signature for ranking and summary
	trait.Df <- read.delim(file = trait.file,sep = "\t",header = T)

	# construc a matrix of trait p-values/correlation
	trait.matrix <- as.matrix(trait.Df[,2:ncol(trait.Df)]);rownames(trait.matrix) <- as.character(trait.Df$module)
	trait.matrix <- trait.matrix[,which(colSums(is.na(trait.matrix)) < nrow(trait.matrix))]

	if (!is.na(sig.file))
	{
		sig.Df <- read.delim(file = sig.file,sep = "\t",header = T)
		
		# construct a matrix of p-values/EFC
		sig.matrix <- get.FET.matrix(FET.table = sig.Df,mod.names = module.ids,
		sig.names = levels(sig.Df$set2_Name),module.col = "set1_Name",id.col = "set2_Name")
		sig.matrix <- cbind(sig.matrix$p.value,sig.matrix$EFC)
	}else{
		sig.matrix <- matrix(0,nrow = length(module.ids),ncol = 0);rownames(sig.matrix) <- module.ids
	}

	# construct a matrix of all statistics
	all.stat <- matrix(NA,nrow = length(module.ids),ncol = (ncol(trait.matrix) + ncol(sig.matrix)))
	rownames(all.stat) <- module.ids;
	if (ncol(sig.matrix) > 0)
	{
		all.stat[match(rownames(sig.matrix),rownames(all.stat)),1:ncol(sig.matrix)] <- sig.matrix;
	}
	
	all.stat[match(rownames(trait.matrix),rownames(all.stat)),(ncol(sig.matrix)+1):ncol(all.stat)] <- trait.matrix;
	
	colnames(all.stat) <- c(colnames(sig.matrix),colnames(trait.matrix))

	# make NA entries in p-value columns to 1, to mitigate NA issues.
	pcol <- grep("P__",colnames(all.stat))
	all.stat[,pcol] <- apply(all.stat[,pcol],2,function(x) {out <- x;out[is.na(out)] <- 1;return(out)})
	
	# extract columns for ranking
	rank.col <- Reduce("union",lapply(paste("\\.P__",exclude.var,sep = ""),function(x,y) grep(x,y),y = colnames(all.stat)))

	# ranked file with signature summary
	module.score <- module_rank.v1(cbind(all.stat[,rank.col]))
	nominal.sigcount <- apply(cbind(all.stat[,rank.col]),1,function(x) sum(x < 0.05,na.rm = T))
	correct.sigcount <- apply(cbind(apply(cbind(all.stat[,rank.col]),2,function(x) p.adjust(x,"BH"))),1,function(x) sum(x < 0.05,na.rm = T))

	# summarize signature FET: nominal < 0.05
	if (!is.na(sig.file))
	{
		sig.summary <- rep(NA,length(module.ids))
		sig.str <- get.FET.summary(FET.table = sig.Df,id.col = "set2_Name",module.col = "set1_Name",pval.col = "FET_pvalue",EFC.col = "enrichment.foldchange",pval.cutoff = 0.05)
		sig.summary[match(names(sig.str),module.ids)] <- sig.str
		nominal.FET <- sig.summary;rm(sig.summary,sig.str)
		
		# summarize signature FET: corrected p-value < 0.05
		sig.summary <- rep(NA,length(module.ids))
		sig.str <- get.FET.summary(FET.table = sig.Df,id.col = "set2_Name",module.col = "set1_Name",pval.col = "corrected.FET.pvalue",EFC.col = "enrichment.foldchange",pval.cutoff = 0.05)
		sig.summary[match(names(sig.str),module.ids)] <- sig.str
		correct.FET <- sig.summary;rm(sig.summary,sig.str)
	
	}else{
		nominal.FET <- correct.FET <- rep(NA,length(module.ids))
	}

	
	# summarize trait correlation :nominal P < 0.05
	trait.summary <- rep(NA,length(module.ids))
	#trait.pval <- apply(trait.matrix[,grep("rho.P__",colnames(trait.matrix))],2,function(x) p.adjust(x,"BH"))
	trait.pval <- trait.matrix[,grep("rho.P__",colnames(trait.matrix))];colnames(trait.pval) <- gsub("^(.*)rho.P__","",colnames(trait.pval))
	trait.str <- apply(trait.pval,1,function(x) paste(paste(names(x)[which(x < 0.05)],format(signif(x[which(x < 0.05)],3),scientific = TRUE),sep = "/"),collapse = ","))
	trait.summary[match(names(trait.str),module.ids)] <- trait.str
	nominal.trait <- trait.summary;rm(trait.pval,trait.summary)
	
	# summarize trait correlation :nominal P < 0.05
	trait.summary <- rep(NA,length(module.ids))
	trait.pval <- apply(trait.matrix[,grep("rho.P__",colnames(trait.matrix))],2,function(x) p.adjust(x,"BH"));colnames(trait.pval) <- gsub("^(.*)rho.P__","",colnames(trait.pval))
	trait.str <- apply(trait.pval,1,function(x) paste(paste(names(x)[which(x < 0.05)],format(signif(x[which(x < 0.05)],3),scientific = TRUE),sep = "/"),collapse = ","))
	trait.summary[match(names(trait.str),module.ids)] <- trait.str
	correct.trait <- trait.summary
	rm(trait.pval,trait.summary)
	
	ranked.summary <- data.frame(modules = module.ids,nominal.count = nominal.sigcount,
	correct.count = correct.sigcount,FET.nominal.summary = nominal.FET,FET.correct.summary = correct.FET,
	trait.nominal.summary = nominal.trait,trait.correct.summary = correct.trait,
	module.rank = rank(module.score),module.rank.score = module.score,as.data.frame(all.stat))

	# combine with existing summary
	if (!is.null(summary.file))
	{
	 module.table <- read.delim(file = summary.file,sep = "\t",header = T)
	 ranked.summary <- combine.table(ranked.summary,module.table)
	}
	
	ranked.summary <- ranked.summary[order(as.numeric(as.character(ranked.summary$module.rank))),]

	return(ranked.summary)
}



generate.rankTable.v2 <- function(trait.file,sig.file,module.file,summary.file,exclude.var)
{
		
	modules <- load.modules(module.file = module.file,file.type = "MEGENA")
	module.ids <- names(modules);rm(modules)

	# process trait + signature for ranking and summary
	trait.Df <- read.delim(file = trait.file,sep = "\t",header = T)

	# construc a matrix of trait p-values/correlation
	trait.matrix <- as.matrix(trait.Df[,2:ncol(trait.Df)]);rownames(trait.matrix) <- as.character(trait.Df$module)
	trait.matrix <- trait.matrix[,which(colSums(is.na(trait.matrix)) < nrow(trait.matrix))]

	if (!is.na(sig.file))
	{
		sig.Df <- read.delim(file = sig.file,sep = "\t",header = T)
		
		# construct a matrix of p-values/EFC
		sig.matrix <- get.FET.matrix(FET.table = sig.Df,mod.names = module.ids,
		sig.names = levels(sig.Df$set2_Name),module.col = "set1_Name",id.col = "set2_Name")
		sig.matrix <- cbind(sig.matrix$p.value,sig.matrix$EFC)
	}else{
		sig.matrix <- matrix(0,nrow = length(module.ids),ncol = 0);rownames(sig.matrix) <- module.ids
	}

	# construct a matrix of all statistics
	all.stat <- matrix(NA,nrow = length(module.ids),ncol = (ncol(trait.matrix) + ncol(sig.matrix)))
	rownames(all.stat) <- module.ids;
	if (ncol(sig.matrix) > 0)
	{
		all.stat[match(rownames(sig.matrix),rownames(all.stat)),1:ncol(sig.matrix)] <- sig.matrix;
	}
	
	trait.row <- intersect(rownames(trait.matrix),rownames(all.stat))
	all.stat[match(trait.row,rownames(all.stat)),(ncol(sig.matrix)+1):ncol(all.stat)] <- trait.matrix[trait.row,];
	
	colnames(all.stat) <- c(colnames(sig.matrix),colnames(trait.matrix))

	# make NA entries in p-value columns to 1, to mitigate NA issues.
	pcol <- grep("P__",colnames(all.stat))
	all.stat[,pcol] <- apply(all.stat[,pcol],2,function(x) {out <- x;out[is.na(out)] <- 1;return(out)})
	
	# extract columns for ranking
	rank.col <- grep("FET\\.P|rho\\.P",colnames(all.stat))
	exclude.var <- c("Clinical.correct.BMI","birth.weight","length.of.gestation","Clinical.lnBaP")
	if (!is.null(exclude.var)) 
	{ 
	 rank.col <- setdiff(rank.col,Reduce("union",lapply(exclude.var,function(x,y) grep(x,y),y = colnames(all.stat))))
	}

	# ranked file with signature summary
	module.score <- module_rank.v1(cbind(all.stat[,rank.col]))
	nominal.sigcount <- apply(cbind(all.stat[,rank.col]),1,function(x) sum(x < 0.05,na.rm = T))
	correct.sigcount <- apply(cbind(apply(cbind(all.stat[,rank.col]),2,function(x) p.adjust(x,"BH"))),1,function(x) sum(x < 0.05,na.rm = T))

	# summarize signature FET: nominal < 0.05
	if (!is.na(sig.file))
	{
		sig.summary <- rep(NA,length(module.ids))
		sig.str <- get.FET.summary(FET.table = sig.Df,id.col = "set2_Name",module.col = "set1_Name",pval.col = "FET_pvalue",EFC.col = "enrichment.foldchange",pval.cutoff = 0.05)
		sig.id <- intersect(names(sig.str),module.ids)
		sig.summary[match(sig.id,module.ids)] <- sig.str[sig.id]
		nominal.FET <- sig.summary;rm(sig.summary,sig.str)
		
		# summarize signature FET: corrected p-value < 0.05
		sig.summary <- rep(NA,length(module.ids))
		sig.str <- get.FET.summary(FET.table = sig.Df,id.col = "set2_Name",module.col = "set1_Name",pval.col = "corrected.FET.pvalue",EFC.col = "enrichment.foldchange",pval.cutoff = 0.05)
		sig.id <- intersect(names(sig.str),module.ids)
		sig.summary[match(sig.id,module.ids)] <- sig.str[sig.id]
		correct.FET <- sig.summary;rm(sig.summary,sig.str)
	
	}else{
		nominal.FET <- correct.FET <- rep(NA,length(module.ids))
	}

	
	# summarize trait correlation :nominal P < 0.05
	trait.summary <- rep(NA,length(module.ids))
	#trait.pval <- apply(trait.matrix[,grep("rho.P__",colnames(trait.matrix))],2,function(x) p.adjust(x,"BH"))
	trait.pval <- trait.matrix[,grep("rho.P__",colnames(trait.matrix))];colnames(trait.pval) <- gsub("^(.*)rho.P__","",colnames(trait.pval))
	trait.str <- apply(trait.pval,1,function(x) paste(paste(names(x)[which(x < 0.05)],format(signif(x[which(x < 0.05)],3),scientific = TRUE),sep = "/"),collapse = ","))
	trait.id <- intersect(names(trait.str),module.ids)
	trait.summary[match(trait.id,module.ids)] <- trait.str[trait.id]
	nominal.trait <- trait.summary;rm(trait.pval,trait.summary)
	
	# summarize trait correlation :nominal P < 0.05
	trait.summary <- rep(NA,length(module.ids))
	trait.pval <- apply(trait.matrix[,grep("rho.P__",colnames(trait.matrix))],2,function(x) p.adjust(x,"BH"));colnames(trait.pval) <- gsub("^(.*)rho.P__","",colnames(trait.pval))
	trait.str <- apply(trait.pval,1,function(x) paste(paste(names(x)[which(x < 0.05)],format(signif(x[which(x < 0.05)],3),scientific = TRUE),sep = "/"),collapse = ","))
	trait.id <- intersect(names(trait.str),module.ids)
	trait.summary[match(trait.id,module.ids)] <- trait.str[trait.id]
	correct.trait <- trait.summary
	rm(trait.pval,trait.summary)
	
	ranked.summary <- data.frame(modules = module.ids,nominal.count = nominal.sigcount,
	correct.count = correct.sigcount,FET.nominal.summary = nominal.FET,FET.correct.summary = correct.FET,
	trait.nominal.summary = nominal.trait,trait.correct.summary = correct.trait,
	module.rank = rank(module.score),module.rank.score = module.score,as.data.frame(all.stat))

	# combine with existing summary
	if (!is.null(summary.file))
	{
	 module.table <- read.delim(file = summary.file,sep = "\t",header = T)
	 ranked.summary <- combine.table(ranked.summary,module.table)
	}
	
	ranked.summary <- ranked.summary[order(as.numeric(as.character(ranked.summary$module.rank))),]

	return(ranked.summary)
}



combine.table <- function(abl,bbl)
{
 common.id <- union(as.character(abl[[1]]),as.character(bbl[[1]]))
 abl.align <- do.call(cbind,lapply(abl[2:ncol(abl)],function(x,y,z) {
                                   names(x) <- z;
                                   vec <- rep(NA,length(y));names(vec) <- y;
								   if (is.numeric(x)) {vec[names(x)] <- x;}else{
								   vec[names(x)] <- as.character(x);}
								   
								   return(vec)
                          },y = common.id,z = as.character(abl[[1]])))
 bbl.align <- do.call(cbind,lapply(bbl[2:ncol(bbl)],function(x,y,z) {
                                   names(x) <- z;
                                   vec <- rep(NA,length(y));names(vec) <- y;
								   if (is.numeric(x)) {vec[names(x)] <- x;}else{
								   vec[names(x)] <- as.character(x);}
								   
								   return(vec)
                          },y = common.id,z = as.character(bbl[[1]])))
 out <- cbind.data.frame(data.frame(id = common.id),abl.align,bbl.align)
 return(out) 
}
