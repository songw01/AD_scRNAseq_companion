library(reshape)
overlap_modules_wt_categ <- function(mods,vec)
{
  tpe = setdiff(unique(vec),NA)
  mat = do.call('rbind',lapply(mods,function(x,y,t) {out = table(y[names(y) %in% x])/length(x);out[match(tpe,names(out))]},y = vec,t = tpe))
  mat[is.na(mat)] = 0;
  #rownames(mat) = gsub("^c1_","M",rownames(mat));
  colnames(mat) = tpe
  
  if (any(rowSums(mat,na.rm = TRUE) < 1))
  {
    mat = cbind(mat,1-rowSums(mat,na.rm = TRUE))
    colnames(mat)[ncol(mat)] = "unassigned"
  }
  return(mat)
}
get_proportion_plot <- function(modules,vec,cols = NULL)
{
  prmat = overlap_modules_wt_categ(mods = modules,vec = vec);
  prdf = melt(prmat);
  colnames(prdf) = c("X1","X2","value")
  
  prdf$X1 = factor(prdf$X1);
  if (!is.null(cols)) prdf$X2 = factor(prdf$X2,levels = names(cols))
  pobj = ggplot() + geom_bar(data = prdf,aes(x = X1,y = value,fill = X2),stat = "identity",position = "stack") 
  if (!is.null(cols)) pobj = pobj + scale_fill_manual(values = cols)
  pobj = pobj + theme_bw() + 
    guides(fill = guide_legend(title = "",ncol = 2)) + 
    theme(axis.title = element_blank(),axis.text = element_text(size = 15),legend.position = "bottom")
  return(pobj)
}

