library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)


#################
## plots in figure 3A
#################

metadata=read_excel("Supplementary_Data_5.xlsx",sheet="Metadata")

hilic=read_excel("Supplementary_Data_5.xlsx",sheet="HILIC LC-MS MS")
hilic=hilic%>%
  filter(`Figure 3`=="yes")
hilic=hilic[,17:47]
hilic=hilic%>%
  gather(-metabolite,key="sample",value="value")

rp=read_excel("Supplementary_Data_5.xlsx",sheet="RP LC-MS MS")
rp=rp%>%
  filter(`Figure 3`=="yes")
rp=rp[,18:48]
rp=rp%>%
  gather(-metabolite,key="sample",value="value")

nmr=read_excel("Supplementary_Data_5.xlsx",sheet="NMR")
nmr=nmr%>%
  filter(`Figure 3`=="yes")
nmr=nmr[,c(1,11:37)]
nmr=nmr%>%
  gather(-metabolite,key="sample",value="value")%>%
  mutate(value=ifelse(is.na(value)==1,0,value))

metabs=full_join(rp,hilic,by=c("sample","value","metabolite"))
metabs=full_join(metabs,nmr,by=c("sample","value","metabolite"))

metabs%>%
  left_join(.,metadata,by="sample")%>%
  filter(day>0)%>%
  ggplot()+
  geom_point(aes(x=day,y=value,color=treatment,fill=treatment))+
  geom_smooth(aes(x=day,y=value,color=treatment,fill=treatment))+
  facet_wrap(~metabolite,scales="free_y")+
  scale_color_manual(values=c("#fe9929","#7bccc4"))+
  scale_fill_manual(values=c("#fe9929","#7bccc4"))+
  theme_classic()

#################
## plots in figure 3B
#################

# read in gene metaT data (access from zenodo repository: 10.5281/zenodo.13375972)
getmm=read.delim("geTMMs_ge5_1X_REV_26samples.txt",header=TRUE,sep="\t")

# read in MAG database annotations (access from zenodo repository: 10.5281/zenodo.13375972)
annotations=read.delim("reactorEMERGE_annotations.tsv",header=TRUE,sep="\t")

# read in metadata
metadata=read_excel("Supplementary_Data_3.xlsx",sheet="Metatranscriptomes")
metadata=metadata%>%
  select(Sample,Treatment,Time)

#make tidy
tidy_getmm=getmm%>%
  gather(-gene,key="sample",value="geTMM")%>%
  filter(geTMM>0)

#add annotations to genes
tidy_getmm_annos=left_join(tidy_getmm,annotations,by=c("gene"="X"))
remove(annotations)

k=tidy_getmm_annos%>%
  select(gene,sample,geTMM,fasta,camper_id,ko_id,camper_rank)%>%
  gather(camper_id,ko_id,key="source",value="annotation")%>%
  filter(annotation!="")

# update annotations with manual curations
k2=k%>%
  mutate(camper_rank=ifelse(annotation=="D00009","A",camper_rank))%>%
  mutate(camper_rank=ifelse(annotation=="D00001","A",camper_rank))%>%
  mutate(camper_rank=ifelse(annotation=="D00006","A",camper_rank))
k2=tidy_getmm_annos%>%filter(gene=="STM_0716_E_M_E034_A_bin.1_k121_848621_4")%>%
  mutate(camper_id="D00001")%>%
  mutate(camper_rank="A")%>%
  select(gene,sample,geTMM,fasta,camper_id,ko_id,camper_rank)%>%
  gather(camper_id,ko_id,key="source",value="annotation")%>%
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
#k4=k3%>%
  left_join(.,camper,by="annotation",relationship = "many-to-many")


# plot dynamics of genes for figure
figure_genes=read_excel("Supplementary_Data_3.xlsx",sheet="catechin_genes")

figure_genes%>%
  left_join(k3,by="annotation")%>%
  left_join(.,metadata,by=c("sample"="Sample"))%>%
  filter(filter!="yes")%>%
  filter(geTMM>0)%>%
  ungroup()%>%
  group_by(sample,Treatment,Time,group)%>%
  summarise(group_sum=sum(geTMM))%>%
  spread(key=group,value=group_sum)%>%
  gather(-sample,-Treatment,-Time,key=group,value=group_sum)%>%
  mutate(group_sum=ifelse(is.na(group_sum)==1,0,group_sum))%>%
  ggplot()+
  geom_point(aes(x=Time,y=group_sum,color=Treatment))+
  geom_smooth(aes(x=Time,y=group_sum,color=Treatment,fill=Treatment))+
  facet_wrap(~group,scales="free_y",ncol = 2)+
  scale_color_manual(values=c("#fe9929","#7bccc4"))+
  scale_fill_manual(values=c("#fe9929","#7bccc4"))+
  theme_classic()
