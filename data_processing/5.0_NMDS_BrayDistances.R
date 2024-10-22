library(writexl)
library(readxl)
library(tidyr)
library(dplyr)
library(edgeR)
library(ggplot2)
library(vegan)

metadata=read_excel("Supplementary_Data_3.xlsx",sheet="Metatranscriptomes")
metadata=metadata%>%
  select(Sample,Treatment,Time)

# read in MAG database annotations (access from zenodo repository: https://doi.org/10.5281/zenodo.13936221)
annotations=read.delim("reactorEMERGE_annotations.tsv",header=TRUE,sep="\t")

# read in gene metaT data (access from zenodo repository: https://doi.org/10.5281/zenodo.13937409)
getmm=read.delim("geTMMs_ge5_1X_REV_26samples.txt",header=TRUE,sep="\t")
#make tidy
tidy_getmm=getmm%>%gather(-gene,key="sample",value="geTMM")%>%filter(geTMM>0)

#add annotations to genes
tidy_getmm_annos=left_join(tidy_getmm,annotations,by=c("gene"="X"))


##################
### MAG level NMDS 
##################
##make table with features as columns and samples as rows
wide_mag_abund=read_excel("Supplementary_Data_3.xlsx",sheet="MAG_metatranscriptome")
wide_mag_abund=as.data.frame(wide_mag_abund)
row.names(wide_mag_abund)=wide_mag_abund[,1]
wide_mag_abund=wide_mag_abund[,-1]
wide_mag_abund=t(wide_mag_abund)

##NMDS on MAG metaT
set.seed(13)
NMDS_Bray_data_mag <-metaMDS(wide_mag_abund, distance = "bray", k=2,
                             autotransform = FALSE, noshare = 0.1, trace = 1)
bray_dist_mag = metaMDSdist(wide_mag_abund, distance = "bray", k=2,
                            autotransform = FALSE, noshare = 0.1, trace = 1)
NMDS_Bray_data_mag$stress
# stress 0.1705884
ord.data = as.data.frame(vegan::scores(NMDS_Bray_data_mag, display="sites"))

# adding metadata
ord.data2=ord.data
ord.data2$sample=row.names(ord.data2)
ord.data2=ord.data2%>%left_join(.,metadata,by=c("sample"="Sample"))
ord.data2$group=paste(ord.data2$Treatment,ord.data2$Time,sep="_")


ord.data2=ord.data2%>%mutate(day=paste("day",Time,sep=""))
a=ord.data2%>%
  ggplot(aes(x=NMDS1,y=NMDS2))+
  geom_point(aes(color=Treatment,shape=day),size=3,alpha=0.8)+
  scale_color_manual(values=c("#fe9929","#7bccc4"))+
  theme_classic()+
  labs(title="metaT_MAG")+
  xlim(-1.75,1.75)+ylim(-1,1)

# significant by group?
set.seed(13)
adonis2(bray_dist_mag ~ ord.data2$Treatment*ord.data2$day, dist="bray",perm=999) #0.001 yes

##################
### gene level NMDS 
##################
##make table with features as columns and samples as rows
#colnames(getmm)
tidy_getmm_wMags=tidy_getmm_annos%>%
  select(gene,fasta,sample,geTMM)%>%
  left_join(.,metadata,by=c("sample"="Sample"))

tidy_mag_geTMM=read_excel("Supplementary_Data_3.xlsx",sheet="MAG_metatranscriptome")%>%
  gather(-MAG,key="sample",value="abund")%>%
  filter(abund>0)

genes_mags22=tidy_mag_geTMM%>%
  select(MAG,sample)%>%
  left_join(.,tidy_getmm_wMags,by=c("MAG"="fasta","sample"))%>%
  select(-MAG,-Treatment,-Time)

genes_mags22=genes_mags22%>%distinct()

genes_mags22=genes_mags22%>%
  spread(key=sample,value=geTMM,fill=0)

genes_mags22=as.data.frame(genes_mags22)
t_getmm=t(genes_mags22[,-1])
t_getmm=as.data.frame(t_getmm)

##NMDS on gene MetaT
set.seed(13)
NMDS_Bray_data_gene <-metaMDS(t_getmm, distance = "bray", k=2,
                              autotransform = FALSE, noshare = 0.1, trace = 1)
bray_dist_gene = metaMDSdist(t_getmm, distance = "bray", k=2,
                             autotransform = FALSE, noshare = 0.1, trace = 1)
NMDS_Bray_data_gene$stress # stress 0.1765114
ord.data = as.data.frame(vegan::scores(NMDS_Bray_data_gene, display="sites"))

# adding metadata
ord.data2=ord.data
ord.data2$sample=row.names(ord.data2)
ord.data2=ord.data2%>%left_join(.,metadata,by=c("sample"="Sample"))
ord.data2$group=paste(ord.data2$Treatment,ord.data2$Time,sep="_")


ord.data2=ord.data2%>%mutate(day=paste("day",Time,sep=""))
b=ord.data2%>%
  ggplot(aes(x=NMDS1,y=NMDS2))+
  geom_point(aes(color=Treatment,shape=day),size=3,alpha=0.8)+
  scale_color_manual(values=c("#fe9929","#7bccc4"))+
  theme_classic()+
  labs(title="metaT_gene")+
  xlim(-1.75,1.75)+ylim(-1,1)

# significant by group?
set.seed(13)
adonis2(bray_dist_gene ~ ord.data2$Treatment*ord.data2$day, dist="bray",perm=999) #0.001 yes

##################
## gene fx metaT
##################
k=tidy_getmm_annos%>%
  select(gene,sample,geTMM,fasta,camper_id,ko_id,camper_rank,cazy_best_hit)%>%
  gather(camper_id,ko_id,cazy_best_hit,key="source",value="annotation")%>%
  filter(annotation!="")
k2=k%>%
  mutate(camper_rank=ifelse(annotation=="D00009","A",camper_rank))%>%
  mutate(camper_rank=ifelse(annotation=="D00001","A",camper_rank))%>%
  mutate(camper_rank=ifelse(annotation=="D00006","A",camper_rank))

k2=tidy_getmm_annos%>%filter(gene=="STM_0716_E_M_E034_A_bin.1_k121_848621_4")%>%
  mutate(camper_id="D00001")%>%
  mutate(camper_rank="A")%>%
  select(gene,sample,geTMM,fasta,camper_id,ko_id,cazy_best_hit,camper_rank)%>%
  gather(camper_id,ko_id,cazy_best_hit,key="source",value="annotation")%>%
  filter(annotation!="")%>%
  full_join(.,k2,by=c("gene","sample","geTMM","fasta","camper_rank","source","annotation"))
# updating rank of manually curated CAMPER hits
k3=k2%>%
  mutate(camper_rank=ifelse(gene=="STM_0716_E_M_E026_E030_E034_A_bin.12_k121_6755434_28","A",camper_rank))%>%
  mutate(camper_rank=ifelse(gene=="STM_0716_E_M_E034_A_bin.1_k121_311991_29","A",camper_rank))%>%
  mutate(camper_rank=ifelse(gene=="STM_0716_E_M_E034_A_bin.1_k121_1799737_13","A",camper_rank))%>%
  mutate(camper_rank=ifelse(gene=="STM_0716_E_M_E034_A_bin.1_k121_1457892_4","A",camper_rank))%>%
  mutate(camper_rank=ifelse(gene=="STM_0716_E_M_E034_A_bin.10_k121_505302_19","A",camper_rank))%>%
  mutate(camper_rank=ifelse(gene=="STM_0716_E_M_E034_A_bin.10_k121_505302_8","A",camper_rank))%>%
  mutate(camper_rank=ifelse(gene=="STM_0716_E_M_E034_A_bin.10_k121_955452_59","A",camper_rank))%>%
  mutate(camper_rank=ifelse(gene=="STM_0716_E_M_E034_A_bin.10_k121_955452_66","A",camper_rank))%>%
  mutate(camper_rank=ifelse(gene=="STM_0716_E_M_E034_A_bin.10_k121_955452_94","A",camper_rank))%>%
  mutate(camper_rank=ifelse(gene=="STM_0716_E_M_E034_A_bin.10_k121_1767847_3","A",camper_rank))%>%
  mutate(camper_rank=ifelse(gene=="STM_0716_E_M_E034_A_bin.10_k121_1879983_13","A",camper_rank))%>%
  mutate(camper_rank=ifelse(gene=="STM_0716_E_M_E034_A_bin.10_k121_1974006_42","A",camper_rank))%>%
  mutate(camper_rank=ifelse(gene=="STM_0716_E_M_E034_A_bin.10_k121_2595405_107","A",camper_rank))
k4=k3%>%
  filter(camper_rank!="B")
# get one annotation per gene- order camper > cazy > kegg
genes_best_anno=k4%>%
  spread(key=source,value=annotation)%>%
  mutate(best_anno=if_else(is.na(camper_id)==0,camper_id,
                           ifelse(is.na(cazy_best_hit)==0,cazy_best_hit,ko_id)))%>%
  ungroup()%>%
  select(gene,sample,best_anno,fasta,geTMM)

# update annotations from phylogenetic trees
fixed_annos=read_excel("Supplementary_Data_3.xlsx",sheet="curated_annotations")

genes_best_anno=genes_best_anno%>%
  left_join(.,fixed_annos,by="gene")%>%
  mutate(best_anno=ifelse(is.na(fx)==1,best_anno,fx))

##make table with features as columns and samples as rows
#colnames(getmm)
wide_fx=genes_best_anno%>%
  left_join(.,metadata,by=c("sample"="Sample"))%>%
  group_by(best_anno,sample)%>%
  summarise(sum=sum(geTMM))%>%
  spread(key=best_anno,value=sum,fill=0)
wide_fx=as.data.frame(wide_fx)
row.names(wide_fx)=wide_fx[,1]
wide_fx=wide_fx[,-1]

##NMDS on fx MetaT
set.seed(13)
NMDS_Bray_data <-metaMDS(wide_fx, distance = "bray", k=2,
                         autotransform = FALSE, noshare = 0.1, trace = 1)
bray_dist_fx = metaMDSdist(wide_fx, distance = "bray", k=2,
                           autotransform = FALSE, noshare = 0.1, trace = 1)
NMDS_Bray_data$stress # stress 0.08661316
ord.data = as.data.frame(vegan::scores(NMDS_Bray_data, display="sites"))

# adding metadata
ord.data2=ord.data
ord.data2$sample=row.names(ord.data2)
ord.data2=ord.data2%>%left_join(.,metadata,by=c("sample"="Sample"))
ord.data2$group=paste(ord.data2$Treatment,ord.data2$Time,sep="_")


ord.data2=ord.data2%>%mutate(day=paste("day",Time,sep=""))
c=ord.data2%>%
  ggplot(aes(x=NMDS1,y=NMDS2))+
  geom_point(aes(color=Treatment,shape=day),size=3,alpha=0.8)+
  scale_color_manual(values=c("#fe9929","#7bccc4"))+
  theme_classic()+
  labs(title="metaT_fx")+
  xlim(-1.75,1.75)+ylim(-1,1)

# significant by group?
set.seed(13)
adonis2(bray_dist_fx ~ ord.data2$Treatment*ord.data2$day, dist="bray",perm=999) #0.001 yes

##################
### 16S level NMDS 
##################
##read in 16S feature table 
data = read_excel("Supplementary_Data_2.xlsx",sheet="16S rRNA gene ASV table")
data=as.data.frame(data)
rownames(data)=data[,1]
data = data[,-1]

# make data with features as columns and samples as rows
t_data=t(data[,1:30])

##NMDS on ASV features
set.seed(13)
NMDS_Bray_data_16S <-metaMDS(t_data, distance = "bray", k=2,
                             autotransform = FALSE, noshare = 0.1, trace = 1,trymax = 100)
bray_dist_16S = metaMDSdist(t_data, distance = "bray", k=2,
                            autotransform = FALSE, noshare = 0.1, trace = 1)
NMDS_Bray_data_16S$stress
# stress = 0.1729862
ord.data = as.data.frame(vegan::scores(NMDS_Bray_data_16S, display="sites"))

# adding metadata
metadata2=read_excel("Supplementary_Data_2.xlsx",sheet="16S rRNA gene sequencing")
ord.data2=ord.data
ord.data2$sample=row.names(ord.data2)
ord.data2=ord.data2%>%left_join(.,metadata2,by=c("sample"="Sample"))
ord.data2$group=paste(ord.data2$Treatment,ord.data2$Time,sep="_")


ord.data2=ord.data2%>%mutate(day=paste("day",Time,sep=""))
# by treatment
d=ord.data2%>%
  ggplot(aes(x=NMDS1,y=NMDS2))+
  geom_point(aes(color=Treatment,shape=day),size=3,alpha=0.8)+
  scale_color_manual(values=c("#fe9929","#7bccc4"))+
  theme_classic()+
  labs(title="16S")+
  xlim(-1.75,1.75)+ylim(-1,1)
# significant by treatment? yes
set.seed(13)
adonis2(bray_dist_16S ~ ord.data2$Treatment*ord.data2$day, perm=999) #0.012 yes


##################
### metabolite level NMDS 
##################
##read in lcms feature table (ranks 1 and 2, data with features as columns and samples as rows)
metadata=read_excel("Supplementary_Data_5.xlsx",sheet="Metadata")

hilic=read_excel("Supplementary_Data_5.xlsx",sheet="HILIC LC-MS MS")
hilic=hilic%>%
  filter(Level=="Level 1"|Level=="Level 2")
hilic=hilic[,17:47]
hilic=hilic%>%
  gather(-metabolite,key="sample",value="value")

rp=read_excel("Supplementary_Data_5.xlsx",sheet="RP LC-MS MS")
rp=rp%>%
  filter(Rank=="Level 1"|Rank=="Level 2")
rp=rp[,18:48]
rp=rp%>%
  gather(-metabolite,key="sample",value="value")
metabs=full_join(rp,hilic,by=c("sample","value","metabolite"))

wide_metabs=metabs%>%
  spread(key=metabolite,value=value,fill=0)
wide_metabs=as.data.frame(wide_metabs)
rownames(wide_metabs)=wide_metabs[,1]
wide_metabs = wide_metabs[,-1]

##NMDS on metabolites
set.seed(13)
NMDS_Bray_data_lcms <-metaMDS(wide_metabs, distance = "bray", k=2,
                              autotransform = FALSE, noshare = 0.1, trace = 1,trymax = 100)
bray_dist_lcms = metaMDSdist(wide_metabs, distance = "bray", k=2,
                             autotransform = FALSE, noshare = 0.1, trace = 1)
NMDS_Bray_data_lcms$stress
# stress = 0.06024526
ord.data = as.data.frame(vegan::scores(NMDS_Bray_data_lcms, display="sites"))

ord.data2=ord.data
ord.data2$sample=row.names(ord.data2)
ord.data2=ord.data2%>%left_join(.,metadata,by=c("sample"))
ord.data2$group=paste(ord.data2$treatment,ord.data2$day,sep="_")


ord.data2=ord.data2%>%mutate(day=paste("day",day,sep=""))
# by treatment
e=ord.data2%>%
  ggplot(aes(x=NMDS1,y=NMDS2))+
  geom_point(aes(color=treatment,shape=day),size=3,alpha=0.8)+
  scale_color_manual(values=c("#fe9929","#7bccc4"))+
  theme_classic()+
  labs(title="Metabolites")+
  xlim(-1.9,1.9)+ylim(-1.2,1.2)
# significant by treatment? yes
set.seed(13)
adonis2(bray_dist_lcms ~ ord.data2$treatment*ord.data2$day, perm=999) #0.012 yes


##################
### plot all plots together
##################
library(cowplot)
plot_grid(d,e,"",b,a,c)
