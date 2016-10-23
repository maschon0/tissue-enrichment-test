##########################
# STATISTICAL THRESHOLDS #
#############################################################################################################################
minimum_fold <- 8              # Smallest fold change to be considered 'specific' to a tissue
testing_correction <- 'fdr'    # P-value adjustment for multiple testing correction
maximum_pvalue <- 0.01         # Largest adjusted P-value to be considered significant in ANOVA
flanking <- TRUE               # If TRUE, buffers temporal changes by including adjacent developmental timepoints in analysis
#############################################################################################################################

cat('Loading reference atlas...\n')
reference_atlas <- read.delim("reference_atlas.tsv",header = T,stringsAsFactors = F)

tissues <- scan("tissues.txt",'character',sep='\n',comment.char = '#')
tissues <- strsplit(tissues,'\t')
names(tissues) <- unlist(lapply(tissues,function(x)x[1]))
tissues <- unlist(lapply(tissues,function(x)x[2]))
tissue_cols <- gsub('^(.+?)_?(.+)\\..*$','\\1',colnames(reference_atlas))
if(!all(tissue_cols%in%names(tissues))){
  cat("WARNING: Unrecognized tissue names! Ignoring samples:\n\t")
  ignr <- colnames(reference_atlas)[!tissue_cols%in%names(tissues)]
  cat(ignr,sep='\t')
  reference_atlas <- reference_atlas[,!colnames(reference_atlas)%in%ignr]
}
excluded_tissues <- names(tissues)[grep('//exclude',tissues)]
if(length(excluded_tissues)>0){
  for(i in 1:length(excluded_tissues)){
    tissue_cols <- gsub('^(.+?)_?(.+)\\..*$','\\1',colnames(reference_atlas))
    tissues_to_exclude <- grep(excluded_tissues[i],tissue_cols)
    cat(paste('Removed',length(tissues_to_exclude),gsub('//exclude','',tissues[excluded_tissues[i]]),'samples from reference atlas.\n',sep=' '))
    x=1:ncol(reference_atlas)
    reference_atlas <- reference_atlas[,x[!x%in%tissues_to_exclude]]
  }
}

timepoints <- scan("timepoints.txt",'character',sep='\n',comment.char = '#')
timepoints <- strsplit(timepoints,'\t')
names(timepoints) <- unlist(lapply(timepoints,function(x)x[1]))
timepoints <- unlist(lapply(timepoints,function(x)x[2]))
timepoint_cols <- gsub('^(.+?)_?(.+)\\..*$','\\2',colnames(reference_atlas))
if(!all(timepoint_cols%in%names(timepoints))){
  cat("WARNING: Unrecognized timepoint names! Ignoring samples:\n\t")
  ignr <- colnames(reference_atlas)[!timepoint_cols%in%names(timepoints)]
  cat(ignr,sep='\t')
  reference_atlas <- reference_atlas[,!colnames(reference_atlas)%in%ignr]
}
excluded_timepoints <- names(timepoints)[grep('//exclude',timepoints)]
if(length(excluded_timepoints)>0){
  for(i in 1:length(excluded_timepoints)){
    timepoint_cols <- gsub('^(.+?)_?(.+)\\..*$','\\2',colnames(reference_atlas))
    timepoints_to_exclude <- grep(excluded_timepoints[i],timepoint_cols)
    cat(paste('Removed',length(timepoints_to_exclude),gsub('//exclude','',timepoints[excluded_timepoints[i]]),'samples from reference atlas.\n',sep=' '))
    x=1:ncol(reference_atlas)
    reference_atlas <- reference_atlas[,x[!x%in%timepoints_to_exclude]]
  }
}

cat("Reference atlas prepared.\nNumber of samples detected at each stage:\n")
tissue_cols <- gsub('^(.+?)_?(.+)\\..*$','\\1',colnames(reference_atlas))
timepoint_cols <- gsub('^(.+?)_?(.+)\\..*$','\\2',colnames(reference_atlas))
for(i in names(timepoints)){
  cat(timepoints[i],'\n',sep='')
  tbl <- table(tissue_cols[grep(i,timepoint_cols)])
  if(length(tbl)==0){
    cat('\tNo samples')
  }else{
    cat(paste('\t',tissues[sort(names(tbl))],': ',tbl[sort(names(tbl))],'\n',sep=''))
  }
}


cat('\nIdentifying tissue-specific gene expression...\n')
specifics=list()
minimum_fold <- log2(minimum_fold)
for(stage in 1:length(timepoints)){
  specifics[[names(timepoints)[stage]]]=list()
  if(flanking){
    samples=which(timepoint_cols%in%as.character((stage-1):(stage+1)))
  }else{
    samples=which(timepoint_cols==as.character(stage))
  }
  subreference_atlas=reference_atlas[,samples]

  for(tissue in names(tissues)){
    tissuereps=grep(paste('^',tissue,sep=''),colnames(subreference_atlas),ignore.case = F,value = T)
    nottissuereps=colnames(subreference_atlas)[!colnames(subreference_atlas)%in%tissuereps]
    if(length(tissuereps)==0){
      next
    }
    cat(paste(tissuereps,collapse=','),'\nvs.\n',paste(nottissuereps,collapse=','),'\n',sep='')
    foldchanges=log2(rowMeans(subreference_atlas[,tissuereps])/rowMeans(subreference_atlas[,nottissuereps]))
    spec_fc=names(which(foldchanges>=minimum_fold))
    subsubreference_atlas=subreference_atlas[spec_fc,c(tissuereps,nottissuereps)]
    pvals=apply(subsubreference_atlas,1,function(x){
      dt=data.frame(values=as.numeric(x),group=c(rep('in',length(tissuereps)),rep('out',length(nottissuereps))),stringsAsFactors = F)
      dt.aov=aov(dt$values ~ dt$group)
      return(as.numeric(unlist(summary(dt.aov))['Pr(>F)1']))
    })
    pvals=p.adjust(pvals,method = testing_correction)
    pvalhits=names(pvals[pvals<maximum_pvalue])
    pvalhits=pvalhits[!pvalhits%in%grep(';',pvalhits,value=T)]
    cat(length(pvalhits),'\n\n')
    specifics[[names(timepoints)[stage]]][[tissue]]<-pvalhits
  }
}

cat('',file = 'specific_genes.txt',append = F)
for(i in names(timepoints)){
  for(j in names(specifics[[i]])){
    cat(paste(paste(i,j,sep='_'),paste(specifics[[i]][[j]],collapse=','),sep='\t'),'\n',sep='',file = 'specific_genes.txt',append=T)
  }
}
cat("Tissue-specific genes saved in 'specific_genes.txt'.")