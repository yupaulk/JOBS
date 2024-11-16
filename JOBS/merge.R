.libPaths("/gpfs/group/dxl46/default/private/lidawang/software/Rpackages")
args = commandArgs(trailingOnly=TRUE)

options(stringsAsFactors=F)
library(dplyr)
library(data.table)
library(tidyr)

#args<-c("22")
chr<-as.numeric(args[1])


print(chr)
cell_type<-c("B_IN","B_Mem","Plasma","CD4_ET","CD4_NC","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","DC1","Mono_C","Mono_NC","NK","NK_R")

df<-as.data.frame(fread(paste0("/storage/group/dxl46/default/private/lidawang/sc_twas/expresso/eqtlgen_map_new/mapA.beta.chr",chr,".txt"),head=T))

df_s<-as.data.frame(fread(paste0("/storage/group/dxl46/default/private/lidawang/sc_twas/expresso/eqtlgen_map_new/mapA.beta.se.chr",chr,".txt"),head=T))


x_old<-df[,3:16]

sum1<-rowSums(x_old)
df$Whole_Blood<-sign(sum1)*abs(df$Whole_Blood)


length(sum1)
length(which(sum1!=0))

df<-df[which(sum1!=0),]

df_s<-df_s[which(sum1!=0),]


dim(df)
head(df)


tmp<-as.data.frame(df[,1])
colnames(tmp)<-"pair"
tmp$gene<-sapply(strsplit(df$pair, "-"), function(x){as.character(x[1])})
tmp$snp<-sapply(strsplit(df$pair, "-"), function(x){as.character(x[2])})
tmp$chr<-sapply(strsplit(tmp$snp, "_"), function(x){as.character(x[1])})
tmp$pos<-sapply(strsplit(tmp$snp, "_"), function(x){as.character(x[2])})
tmp$ref<-sapply(strsplit(tmp$snp, "_"), function(x){as.character(x[3])})
tmp$alt<-sapply(strsplit(tmp$snp, "_"), function(x){as.character(x[4])})


df[df==0]<-NA
df_s[df_s==0]<-NA
colnames(df)[2:ncol(df)]<-paste0("BETA_",colnames(df)[2:ncol(df)])
colnames(df_s)[2:ncol(df)]<-paste0("SE_",colnames(df_s)[2:ncol(df)])

df_out<-merge(tmp,df,by="pair")
df_out<-merge(df_out,df_s,by="pair")
dim(df_out)


fwrite(df_out,paste0("/storage/group/dxl46/default/private/datasets_refs/OneK1K_sumstat/eQTLGen/MERGE_eQTLGen_OneK1K_chr",chr,".txt"),col.names = T,row.names = F,sep="\t",quote=F)


