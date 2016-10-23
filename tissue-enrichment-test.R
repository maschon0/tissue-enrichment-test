################
# USER OPTIONS #
#############################################################################################################################
filename <- 'nodine_2012_embryos' # Name of the .TSV file containing transcriptomes to be tested. Try also: 'reference_atlas'
samples <- 'all'                  # Either 'all', a vector of names, or a vector of numbers, indicating which samples in
                                  # 'filename.description' to test
make_heatmaps <- TRUE             # If TRUE, generates a heatmap of p-values for all samples. Requires package 'pheatmap'
row_breaks <- c(2,5)              # A vector containing row numbers where a gap should be place on the heatmap
col_breaks <- FALSE               # A vector containing column numbers where a gap should be place on the heatmap
tissue_colors=c('#609A33','#C0D84E','#C56DAC','#F9C9DE','#F26D7E','#2A89B6','#66C5E9') # For CDF plots; one color per tissue
names(tissue_colors)=c("EP","SUS","MCE","PEN","CZE","CZSC","GSC")                      # All tissues used in the test
#############################################################################################################################

trim <- function(x){
  gsub("^\\s+|\\s+$", "", x)
}

if(make_heatmaps==TRUE){
  install.packages(c('RColorBrewer','pheatmap'),lib = 'Rpackages',repos = "http://cran.us.r-project.org")
  library(RColorBrewer,lib.loc = 'Rpackages')
  library(pheatmap,lib.loc = 'Rpackages')
}


tissues <- scan("tissues.txt",'character',sep='\n',comment.char = '#')
tissues <- strsplit(tissues,'\t')
names(tissues) <- unlist(lapply(tissues,function(x)x[1]))
tissues <- unlist(lapply(tissues,function(x)x[2]))

timepoints <- scan("timepoints.txt",'character',sep='\n',comment.char = '#')
timepoints <- strsplit(timepoints,'\t')
names(timepoints) <- unlist(lapply(timepoints,function(x)x[1]))
timepoints <- unlist(lapply(timepoints,function(x)x[2]))

specific_genes <- scan("specific_genes.txt",'character',sep='\n')
specific_samples <- unlist(lapply(strsplit(specific_genes,'\t'),function(x)x[1]))
specific_genelist <- unlist(lapply(strsplit(specific_genes,'\t'),function(x)x[2]))
specifics=list()
for(i in 1:length(specific_genes)){
  timepoint=gsub('^([0-9]+)_(.+)$','\\1',specific_samples[i])
  tissue=gsub('^([0-9]+)_(.+)$','\\2',specific_samples[i])
  specifics[[timepoint]][[tissue]]=unlist(strsplit(specific_genelist[i],','))
}
rm(specific_genes,specific_samples,specific_genelist)

transcriptomes <- read.delim(paste(filename,'.tsv',sep=''),header = T,stringsAsFactors = F)
transcriptome_descriptions <- read.delim(paste(filename,'.description',sep=''),header = T,stringsAsFactors = F)
rownames(transcriptome_descriptions) <- transcriptome_descriptions[,1]
if(tolower(samples)=='all'){
  samples=rownames(transcriptome_descriptions)
}
if(is.numeric(samples)){
  samples=rownames(transcriptome_descriptions)[samples]
}

###################################################################################
perform_contamination_test=function(nm,columns,stages,gaps_row,gaps_col,output_name,make_heatmaps){
  outputdir <- 'results'
  KS_tissuetypes=function(seqdata,test_timepoint,random_reps,poolsize,low,high,stepby,colorset,width,legendlabel,axes){
    # Performs a one-sided KS Test on a set of genes and their expression values for each tissue-specific gene set.
    # Required inputs: 
    #          seqdata- a named numeric vector containing an expression value for each gene
    #          test_timepoint- developmental timepoint to compare to. Options in 'timepoints.txt'
    #          random_reps- how many random samples of genes to take as control group
    #          poolsize- the number of genes per random sample. Default is the length of each tissue-specific gene set
    #          low- minimum expression value. Uses the min of seqdata by default
    #          high- Maximum value for thresholding. Default is the 95th expression percentile
    #          stepby- How granular the threshold sampling should be. For expression percentiles, default is 100.
    
    set.seed(8675309)
    freqs=list()
    stepby=abs(stepby)
    thresh=seq(high,low,-stepby)
    if(thresh[length(thresh)]!=low)thresh=append(thresh,low)
    for(i in 1:length(specifics[[test_timepoint]])){
      tissue=names(specifics[[test_timepoint]])[i]
      x=specifics[[test_timepoint]][[i]]
      x=x[x%in%names(seqdata)]
      x=seqdata[x]
      exvec=vector(length=length(thresh))
      count=1
      for(threshold in thresh){
        exvec[count]=length(which(x>=threshold))/length(x)
        count=count+1
      }
      freqs[[tissue]]=exvec
    }
    RANDOM=list()
    if(poolsize=='samesize'){
      ps=length(specifics[[test_timepoint]][[tissue]])
    }else{
      ps=poolsize
    }
    for(z in 1:random_reps){
      random_names=sample(names(seqdata),ps)
      randoms=seqdata[random_names]
      exvec=vector(length=length(thresh))
      count=1
      for(threshold in thresh){
        exvec[count]=length(which(randoms>=threshold))/length(randoms)
        count=count+1
      }
      RANDOM[[z]]=exvec
    }
    par(mar=c(5,5,3,2),cex=1,lwd=2)
    plot(c(-10,-10),xlim=c(1,length(thresh)),ylim=c(0,1),xaxt="n",yaxt="n",xlab="Expression Threshold (Percentile)",ylab="Percent of Genes Expressed",bty="o")
    
    axis(1,c(1,seq(round(length(thresh)/4,0),length(thresh),length.out = 4)),c("",seq(75,0,length.out = 4)),lwd = 2)
    if(axes=='xy')axis(2,seq(0,1,.1),seq(0,100,10),las=2,lwd = 2,hadj = .7)
    random_matrix=matrix(unlist(RANDOM),nrow = length(RANDOM),byrow = T)
    sds=apply(random_matrix,2,function(x)sd(x))
    CI95=1.96*sds/sqrt(10)
    rmeans=colMeans(random_matrix)
    ks.pvals=list()
    for(name in names(freqs)){
      current=freqs[[name]]
      lines(current,lwd=width,col=tissue_colors[name])
      # x=ks.test(x = current[1:(length(thresh)-1)],y = rmeans[1:(length(thresh)-1)],alternative = "less",exact = T)$p.value
      dfs=which(diff(current)>0)
      x=ks.test(x = current[dfs][1:(length(dfs)-1)],y = rmeans[dfs][1:(length(dfs)-1)],alternative = "less",exact = T)$p.value
      if(x<.05){star="*"}else{star=""}
      if(x<.01){star="**"}
      if(x<.001){star="***"}
      ks.pvals[[name]]=paste(format(x,digits = 3,scientific = x<.05),star,sep="")
    }
    tissuetypes=names(ks.pvals)
    if(legendlabel){legend("topleft",legend=c("One-sided K-S Test",paste(tissuetypes,ks.pvals[trim(tissuetypes)],sep=": ")),bty="n",pch=15,pt.cex = 1.5,col=c("white",tissue_colors[trim(tissuetypes)]))}
    else{legend("topleft",legend=paste(tissuetypes,ks.pvals[trim(tissuetypes)],sep=" "),bty="n",pch=15,cex=.9,pt.cex = 1.2,col=tissue_colors[trim(tissuetypes)])}
    lines(rmeans,lwd=width,col="black")
    polygon(x=c(1:length(rmeans),length(rmeans):1),y=c(rmeans+sds,(rmeans-sds)[length(rmeans):1]),col=adjustcolor("black",alpha.f = .25),border = F)
    z=as.numeric(gsub("(.+?)[\\*]*","\\1",ks.pvals))
    names(z)=names(ks.pvals)
    return(z)
  }
  #####################################################################
  y=matrix(nrow=length(unique(unlist(lapply(specifics,names)))),ncol=length(columns))
  rownames(y)=unique(unlist(lapply(specifics,names)))
  for(i in 1:length(columns)){
    cat("Performing tissue-enrichment test on ",columns[i],'\n')
    pdf(paste(outputdir,'/',columns[i],".pdf",sep=""),4,6)
    seqdata=eval(parse(text = nm))[,columns[i]]
    seqdata=as.numeric(trim(seqdata))
    names(seqdata)=rownames(eval(parse(text = nm)))
    seqdata=seqdata[!(is.na(seqdata))]
    p= round(mean(unlist(lapply(specifics[[stages[i]]],length))),0)
    h= as.numeric(sort(seqdata,decreasing = T)[round(length(seqdata)*.05,0)])
    l= min(seqdata)
    z=KS_tissuetypes(seqdata,test_timepoint = stages[i]
                     ,random_reps = 100,poolsize = p
                     ,high = h,low = l,
                     stepby = min(diff(seq(h,l,length.out = 100))),width = 4,legendlabel = F,axes = 'xy',
                     colorset="")
    y[names(z),i]=as.numeric(z)
    title(main=columns[i],cex=.6)
    dev.off()
  }
  colnames(y)=columns
  cat('Storing p-values in ',paste(outputdir,paste(output_name,'p_values.txt',sep="_"),sep='/'),'\n',sep='')
  write.table(y,paste(outputdir,paste(output_name,'p_values.txt',sep="_"),sep='/'),sep="\t",quote=F)
  if(make_heatmaps){
    cat('Saving heatmap to ',paste(outputdir,paste(output_name,"_heatmap.pdf",sep=''),sep='/'),'\n',sep='')
    breaksList = c(0,-log10(0.05),-log10(0.01),3,4,5,10,15,20,50)
    pdf(paste(outputdir,paste(output_name,"_heatmap.pdf",sep=''),sep='/'),7,4,onefile = F)
    pheatmap(-log10(y[1:7,]),cluster_rows = F,cluster_cols = F,scale = 'none', breaks = breaksList,color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)),
             border_color = 'black',gaps_col = gaps_col,
             gaps_row = gaps_row,main="Tissue-Specific Gene Enrichment",
             labels_row = tissues[rownames(y)],
             labels_col = colnames(y),fontsize = 6)
    dev.off()
  }
}
#####################################################################################
  
perform_contamination_test(nm='transcriptomes',
                             columns=samples,
                             stages=transcriptome_descriptions[,2],
                             gaps_row=row_breaks,
                             gaps_col=col_breaks,
                             output_name=filename,
                             make_heatmaps=make_heatmaps)
