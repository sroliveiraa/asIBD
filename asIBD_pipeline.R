


############################################################### Ancestry-specific IBD - Sandra Oliveira - 25 Nov 2021

# To run this script you will need:
# 1) local ancestry results from RFMix v2
# 2) IBD sharing results from refined IBD. Use the raw output, without the merging step where gaps are removed, because in that step the phase information is lost and you cannot know which of the two chromosomes is part of the IBD.

# In the publication "Genome wide variation in the Angolan Namib desert reveals unique Pre-Bantu ancestry" from Oliveira et al., the RFMix analysis was conducted with 3 source populations, coded as:
# East_Africa=0,	South_Africa=1,	and West_Africa=2



############################################################### Step 1 - Add ancestry assignment to each IBD

setwd("")

interp=function(startp, endp, sg, eg, newp){
  return((((newp-startp)*(eg-sg))/(endp-startp))+sg)
}

for (chr in 1:22){
  #read your input files
  rfmix=read.table(paste0("filename.chr",chr,".msp.tsv"), sep="\t", head=T, comment.char = "@", skip = 1)
  ibd=read.table(paste0("filename.chr",chr,".ibd"), sep="", head=F)

  #Below I add a few more columns to the IBD file containing the ancestry-specific info:
  #columns 10-12th store information on the ancestry assigned by RFMix to a particular IBD
  #columns 10-12th represent the proportion of East, South, and West African-related ancestries, respectively, in that particular IBD.
  #columns 10-12th must sum to 1
  #column 13th has a tag "U", "M", or "I", which classifies the IBD as having a UNIQUE ancestry, MULTIPLE ancestries, or INCONSISTENT ancestries.
  #In the case of MULTIPLE ancestries, both chromosomes still need to have a perfect ancestry match along the whole IBD.
  #In the case of INCONSISTENT ancestries, the ancestry assigned to the two haplotypes that are in IBD are not exactly the same
  #The inconsistencies should actually correspond to errors in the local ancestry assignment and they often occur at the beginning or end of an IBD
  #For these cases, I still track how much of the IBD is from each ancestry, you can think of it as the probability of the IBD being of ancestry 0, 1, or 2.
  #In a later step you can define a cutoff for what you accept as informative IBDs, and what you exclude for being too inconsistent.
  
  #Some examples of results for columns 10-13th
  # "0 1 0 U"       the IBD is completely of ancestry 1, in this case Southern Africa (column 11th)
  # "0 0.9 0.1 M"   the IBD is 90% of ancestry 1 and 10% of ancestry 2 (both ancestries where inherited together by these individuals) but consistent between the two chromosomes
  # "0.5 0.5 0 M"   the IBD is 50% of ancestry 0 and 50% of ancestry 1 (both ancestries where inherited together by these individuals) but consistent between the two chromosomes
  # "0 0.9 0.1 I"   the IBD is inconsistent between the two chromosomes but most of it (90%) is of ancestry 1. 90% of its length could thus be accounted for when reporting the sum total lengths per ancestry. See next step.
  # "0.5 0.5 0 I"   the IBD is inconsistent between the two chromosomes, with 50% ancestry assigned to pop 0 and 50% to pop 1. We definitely don't want to use these IBDs and should treat them as ambiguous.
  
  anc_0=rep(NA, length(ibd[,1]))
  anc_1=anc_0
  anc_2=anc_0
  tag=anc_0
  ibd=cbind(ibd, anc_0, anc_1, anc_2,tag)
  #length in cM
  len=apply(rfmix,1,function(x){x[5]-x[4]})
  half_len=len/2
  
  #haplotype 1 and 2 in the IBD files correspond to haplotype 0 and 1 in RFMix
  ibd[ibd[,2]==1,2]=0
  ibd[ibd[,4]==1,4]=0
  ibd[ibd[,2]==2,2]=1
  ibd[ibd[,4]==2,4]=1
  
  #assumes snps are ordered in RFMix file
  for(i in 1:length(ibd[,1])){
    #get ancestry assignment for the positions that span the IBD for each of the two chromosomes
    fc=rfmix[max(which(rfmix$spos<=ibd[i,6])):min(which(rfmix$epos>=ibd[i,7])),colnames(rfmix) %in% paste0(ibd[i,1],".",ibd[i,2])]
    sc=rfmix[max(which(rfmix$spos<=ibd[i,6])):min(which(rfmix$epos>=ibd[i,7])),colnames(rfmix) %in% paste0(ibd[i,3],".",ibd[i,4])]
    
    #Are the assignments matching for the two chromosomes?
    if(isTRUE(all.equal(fc,sc))){
      #Unique ancestry IBD
      #these are the cases in which the IBD is completely of one ancestry, I tag it as "U", in column 13th
      if(length(unique(c(fc,sc)))==1){
        ibd[i,colnames(ibd)==paste0("anc_",unique(c(fc,sc)))]=1
        ibd$tag[i]="U"
      }
      #Multiple ancestry IBD
      #these are the cases in which the IBD is made up of more than one ancestry but the ancestry is consistent (ie exactly the same) in both haplotypes, I tag it as "M", in column 13th
      else{
        mat=cbind(fc,sc)
        endmap=rfmix[min(which(rfmix$epos>=ibd[i,7])),5]
        startmap=rfmix[max(which(rfmix$spos<=ibd[i,6])),4]
        temp_half_len=half_len
        temp_len=len
        #if the start and end positions of the IBD do not match the RFMix crf windows, estimate the genetic positions (cM) via interpolation
        if(!ibd[i,7] %in% rfmix$epos){
          endmap=interp(rfmix[min(which(rfmix$epos>=ibd[i,7])),2], rfmix[min(which(rfmix$epos>=ibd[i,7])),3], rfmix[min(which(rfmix$epos>=ibd[i,7])),4], rfmix[min(which(rfmix$epos>=ibd[i,7])),5],ibd[i,7])
          temp_len[min(which(rfmix$epos>=ibd[i,7]))]=endmap-rfmix[min(which(rfmix$epos>=ibd[i,7])),4]
          temp_half_len[min(which(rfmix$epos>=ibd[i,7]))]=temp_len[min(which(rfmix$epos>=ibd[i,7]))]/2
        }
        if(!ibd[i,6] %in% rfmix$spos){
          startmap=interp(rfmix[max(which(rfmix$spos<=ibd[i,6])),2], rfmix[max(which(rfmix$spos<=ibd[i,6])),3], rfmix[max(which(rfmix$spos<=ibd[i,6])),4], rfmix[max(which(rfmix$spos<=ibd[i,6])),5],ibd[i,6])
          temp_len[max(which(rfmix$spos<=ibd[i,6]))]=rfmix[max(which(rfmix$spos<=ibd[i,6])),5]-startmap
          temp_half_len[max(which(rfmix$spos<=ibd[i,6]))]=temp_len[max(which(rfmix$spos<=ibd[i,6]))]/2
        }
        total_crf=endmap-startmap
        #estimate the proportion for two out of three ancestries. The third is what is left.
        for(k in 1:2){
          minor=as.numeric(names(sort(table(c(fc,sc)), decreasing=F)[k]))
          minor_i1=apply(mat,1, function(x){sum(x==c(minor,minor))})==1
          minor_i2=apply(mat,1, function(x){sum(x==c(minor,minor))})==2
          minor_crf=sum(temp_half_len[max(which(rfmix$spos<=ibd[i,6])):min(which(rfmix$epos>=ibd[i,7]))][minor_i1])+sum(temp_len[max(which(rfmix$spos<=ibd[i,6])):min(which(rfmix$epos>=ibd[i,7]))][minor_i2])
          ibd[i,colnames(ibd)==paste0("anc_",minor)]=minor_crf/total_crf
          ibd$tag[i]="M"
        }
      }
    }
    #Inconsistent ancestry IBD
    #If the two haplotypes that are in IBD are not exactly the same, I tag it as inconsistent ("I")
    #For these cases, I still track how much of the IBD is from each ancestry.
    else{
      mat=cbind(fc,sc)
      endmap=rfmix[min(which(rfmix$epos>=ibd[i,7])),5]
      startmap=rfmix[max(which(rfmix$spos<=ibd[i,6])),4]
      temp_half_len=half_len
      temp_len=len
      #if the start and end positions of the IBD do not match the RFMix crf windows, estimate the genetic positions (cM) via interpolation
      if(!ibd[i,7] %in% rfmix$epos){
        endmap=interp(rfmix[min(which(rfmix$epos>=ibd[i,7])),2], rfmix[min(which(rfmix$epos>=ibd[i,7])),3], rfmix[min(which(rfmix$epos>=ibd[i,7])),4], rfmix[min(which(rfmix$epos>=ibd[i,7])),5],ibd[i,7])
        temp_len[min(which(rfmix$epos>=ibd[i,7]))]=endmap-rfmix[min(which(rfmix$epos>=ibd[i,7])),4]
        temp_half_len[min(which(rfmix$epos>=ibd[i,7]))]=temp_len[min(which(rfmix$epos>=ibd[i,7]))]/2
      }
      if(!ibd[i,6] %in% rfmix$spos){
        startmap=interp(rfmix[max(which(rfmix$spos<=ibd[i,6])),2], rfmix[max(which(rfmix$spos<=ibd[i,6])),3], rfmix[max(which(rfmix$spos<=ibd[i,6])),4], rfmix[max(which(rfmix$spos<=ibd[i,6])),5],ibd[i,6])
        temp_len[max(which(rfmix$spos<=ibd[i,6]))]=rfmix[max(which(rfmix$spos<=ibd[i,6])),5]-startmap
        temp_half_len[max(which(rfmix$spos<=ibd[i,6]))]=temp_len[max(which(rfmix$spos<=ibd[i,6]))]/2
      }
      total_crf=endmap-startmap
      #estimate the proportion for two out of three ancestries. The third is what is left.
      for(k in 1:2){
        minor=as.numeric(names(sort(table(c(fc,sc)), decreasing=F)[k]))
        minor_i1=apply(mat,1, function(x){sum(x==c(minor,minor))})==1
        minor_i2=apply(mat,1, function(x){sum(x==c(minor,minor))})==2
        minor_crf=sum(temp_half_len[max(which(rfmix$spos<=ibd[i,6])):min(which(rfmix$epos>=ibd[i,7]))][minor_i1])+sum(temp_len[max(which(rfmix$spos<=ibd[i,6])):min(which(rfmix$epos>=ibd[i,7]))][minor_i2])
        ibd[i,colnames(ibd)==paste0("anc_",minor)]=minor_crf/total_crf
        ibd$tag[i]="I"
      }
    }
    #fill-in the table for the missing ancestry
    ibd[i,10:12][is.na(ibd[i,10:12])]=1-sum(ibd[i,10:12], na.rm = T)
  }
  write.table(ibd, paste0("filename.chr",chr,".ibd_la"), sep=" ", quote=F, col.names = F, row.names = F)
}



############################################################### Step 2 - If desired, group IBD's by category (in bash)

cat filename.chr*.ibd | awk '$9 > 1 && $9 <= 5' > filename.ibd.1to5cM
cat filename.chr*.ibd | awk '$9 > 5 && $9 <= 10' > filename.ibd.5to10cM
cat filename.chr*.ibd | awk '$9 > 10' > filename.ibd.over10cM

cat filename.chr*.ibd_la | awk '$9 > 1 && $9 <= 5' > filename.ibd_la.1to5cM
cat filename.chr*.ibd_la | awk '$9 > 5 && $9 <= 10' > filename.ibd_la.5to10cM
cat filename.chr*.ibd_la | awk '$9 > 10' > filename.ibd_la.over10cM



############################################################### Step 3 - compute pairwise matrices with the IBD summed lengths per ancestry

#define here the pop names
pop=""
#define here the number of chromosomes per pop (ie 2*number of individuals) in the same order
pop_count=""

#define paths to files (IBDs by category for standard IBD)
ibd1="filename.ibd.1to5cM"
ibd2="filename.ibd.5to10cM"
ibd3="filename.ibd.over10cM"

#define paths to files (IBD by category for ancestry-specific IBD)
asibd1="filename.ibd_la.1to5cM"
asibd2="filename.ibd_la.5to10cM"
asibd3="filename.ibd_la.over10cM"


#here is the function for the standard IBD (ie the without considering the local ancestry)
f_ibd_stats=function(x, pop_count, pop, tag){
  
  #read IBD's
  ibd=read.table(x, sep="", head=F)

  #add population pair column at the end
  pair=apply(ibd,1, function(x){ paste0(sort(c(strsplit(x[1],'_')[[1]][1], strsplit(x[3],'_')[[1]][1])), collapse = ":") })
  #create matrix with stats: the sum of IBD lengths (or number of IBD's) that are shared on average by a pair of individuals (haploid)
  #the matrix has all the pairwise comparisons
  number_ibd=matrix(NA, ncol=length(pop), nrow=length(pop))
  length_ibd=matrix(NA, ncol=length(pop), nrow=length(pop))
  colnames(number_ibd)= pop
  rownames(number_ibd)= pop
  colnames(length_ibd)= pop
  rownames(length_ibd)= pop
  for(i in 1:length(pop)){
    for(j in 1:length(pop)){
      if(i!=j){
        number_ibd[i,j]=length(ibd[which(ibd[,10] == paste0(sort(pop[c(i,j)]),collapse = ":")),9])/as.vector(pop_count[pop[i]]*pop_count[pop[j]])
        length_ibd[i,j]=sum(ibd[which(ibd[,10] == paste0(sort(pop[c(i,j)]),  collapse = ":")),9])/as.vector(pop_count[pop[i]]*pop_count[pop[j]])
      }
      else{
        number_ibd[i,j]=length(ibd[which(ibd[,10] == paste0(sort(pop[c(i,j)]),collapse = ":")),9])/as.vector(choose(pop_count[pop[i]],2))
        length_ibd[i,j]=sum(ibd[which(ibd[,10] == paste0(sort(pop[c(i,j)]),  collapse = ":")),9])/as.vector(choose(pop_count[pop[i]],2))
      }
    }
  }
  write.table(number_ibd, paste0(tag,"_number_ibd_allpops.csv"), sep=" ", col.names = T, row.names = T, quote = F)
  write.table(length_ibd, paste0(tag,"_length_ibd_allpops.csv"), sep=" ", col.names = T, row.names = T, quote = F)
  l=list()
  l[[1]]=number_ibd
  l[[2]]=length_ibd
  return(l)
  
}


#here is the function for the ancestry specific IBD
f_asibd_stats=function(x, pop_count, pop, tag, inc){
  
  #read IBD's
  ibd=read.table(x, sep="", head=F)
  
  #deal with inconsistencies
  #i have defined 3 different cutoffs to deal with inconsistencies: "strict", "moderate", and "all"
  #You can define the one you chose when you call this function below
  #strict: no inconsistencies
  if(inc=="strict"){
    ibd=ibd[ibd[,13] != "I",]
  }
  #moderate: take into account inconsistent IBDs that have more that 75% of one ancestry
  #the proportion of the total IBD length for each ancestry will be used for the sums
  if(inc=="moderate"){
    m=apply(ibd[,10:12], 1, function(x){sort(x, decreasing = T)[1]})
    ibd=rbind(ibd[ibd[,13] != "I",], ibd[intersect(which(ibd[,13] == "I"),which(m>0.75)),])
  }
  #if you chose "all" then no IBDs will be removed (not recommended)
  
  #add population pair column at the end
  pair=apply(ibd,1, function(x){ paste0(sort(c(strsplit(x[1],'_')[[1]][1], strsplit(x[3],'_')[[1]][1])), collapse = ":") })
  #create matrix with stats: the sum of IBD lengths (or number of IBD's) that are shared on average by a pair of individuals (haploid)
  #East_Africa=0	South_Africa=1	West_Africa=2
  number_ibd=matrix(NA, ncol=length(pop), nrow=length(pop))
  length_ibd=matrix(NA, ncol=length(pop), nrow=length(pop))
  colnames(number_ibd)= pop
  rownames(number_ibd)= pop
  colnames(length_ibd)= pop
  rownames(length_ibd)= pop
  number_ibd1=number_ibd
  length_ibd1=length_ibd
  number_ibd2=number_ibd
  length_ibd2=length_ibd
  for(i in 1:length(pop)){
    for(j in 1:length(pop)){
      if(i!=j){
        number_ibd[i,j]=sum(ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),collapse = ":")),10])/as.vector(pop_count[pop[i]]*pop_count[pop[j]])
        number_ibd1[i,j]=sum(ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),collapse = ":")),11])/as.vector(pop_count[pop[i]]*pop_count[pop[j]])
        number_ibd2[i,j]=sum(ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),collapse = ":")),12])/as.vector(pop_count[pop[i]]*pop_count[pop[j]])
        length_ibd[i,j]=sum(ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),  collapse = ":")),9]*ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),  collapse = ":")),10])/as.vector(pop_count[pop[i]]*pop_count[pop[j]])
        length_ibd1[i,j]=sum(ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),  collapse = ":")),9]*ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),  collapse = ":")),11])/as.vector(pop_count[pop[i]]*pop_count[pop[j]])
        length_ibd2[i,j]=sum(ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),  collapse = ":")),9]*ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),  collapse = ":")),12])/as.vector(pop_count[pop[i]]*pop_count[pop[j]])
      }
      else{
        number_ibd[i,j]=sum(ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),collapse = ":")),10])/as.vector(choose(pop_count[pop[i]],2))
        number_ibd1[i,j]=sum(ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),collapse = ":")),11])/as.vector(choose(pop_count[pop[i]],2))
        number_ibd2[i,j]=sum(ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),collapse = ":")),12])/as.vector(choose(pop_count[pop[i]],2))
        length_ibd[i,j]=sum(ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),  collapse = ":")),9]*ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),  collapse = ":")),10])/as.vector(choose(pop_count[pop[i]],2))
        length_ibd1[i,j]=sum(ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),  collapse = ":")),9]*ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),  collapse = ":")),11])/as.vector(choose(pop_count[pop[i]],2))
        length_ibd2[i,j]=sum(ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),  collapse = ":")),9]*ibd[which(ibd[,14] == paste0(sort(pop[c(i,j)]),  collapse = ":")),12])/as.vector(choose(pop_count[pop[i]],2))
      }
    }
  }
  #0, 1, and 2 in the filename indicate the ancestry
  write.table(number_ibd, paste0(tag,"_number0_ibd.csv"), sep=" ", col.names = T, row.names = T, quote = F)
  write.table(length_ibd, paste0(tag,"_length0_ibd.csv"), sep=" ", col.names = T, row.names = T, quote = F)
  write.table(number_ibd1, paste0(tag,"_number1_ibd.csv"), sep=" ", col.names = T, row.names = T, quote = F)
  write.table(length_ibd1, paste0(tag,"_length1_ibd.csv"), sep=" ", col.names = T, row.names = T, quote = F)
  write.table(number_ibd2, paste0(tag,"_number2_ibd.csv"), sep=" ", col.names = T, row.names = T, quote = F)
  write.table(length_ibd2, paste0(tag,"_length2_ibd.csv"), sep=" ", col.names = T, row.names = T, quote = F)
  l=list()
  l[[1]]=number_ibd
  l[[2]]=length_ibd
  l[[3]]=number_ibd1
  l[[4]]=length_ibd1
  l[[5]]=number_ibd2
  l[[6]]=length_ibd2
  return(l)
}


#here the function is called for the stardard IBD
res1=f_ibd_stats(ibd1, pop_count, pop, "stats_1to5cM")
res2=f_ibd_stats(ibd2, pop_count, pop, "stats_5to10cM")
res3=f_ibd_stats(ibd3, pop_count, pop, "stats_over10cM")

#here the function is called for the ancestry specific
#the last parameter defines the cutoff ("strict", "moderate", "all")
#the second last is just the name for the output file
res1=f_asibd_stats(asibd1, pop_count, pop, "stats_1to5cM_moderate", "moderate")
res2=f_asibd_stats(asibd2, pop_count, pop, "stats_5to10cM_moderate", "moderate")
res3=f_asibd_stats(asibd3, pop_count, pop, "stats_over10cM_moderate", "moderate")


#Tip for plotting the results
#Read the three matrices with the ancestry-specific IBD results and the one with the standard IBDs, and plot the results as stacked bar plots 
#3 cols for the asIBDs and grey for the standard IBDs minus the asIBDs. The grey would then show the proportion of IBDs you excluded because of ambiguities
