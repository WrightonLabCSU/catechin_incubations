library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

# consider specific catechin-related genes
figure_genes=read_excel("Supplementary_Data_3.xlsx",sheet="catechin_genes")

# read in gene metaT data (access from zenodo repository: 10.5281/zenodo.13375972)
getmm=read.delim("geTMMs_ge5_1X_REV_26samples.txt",header=TRUE,sep="\t")

# read in MAG database annotations (access from zenodo repository: 10.5281/zenodo.13375972)
annotations=read.delim("reactorEMERGE_annotations.tsv",header=TRUE,sep="\t")

# read in MAG taxonomy
mags=read_excel("Supplementary_Data_2.xlsx",sheet="MAG Database")
mags=mags%>%
  select(MAG,`GTDB v2.3.0 r214`)

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


# plot dot plot of MAGs
figure_genes%>%
  left_join(k3,by="annotation")%>%
  left_join(.,metadata,by=c("sample"="Sample"))%>%
  filter(Treatment=="catechin")%>%
  filter(filter!="yes")%>%
  ungroup()%>%
  select(sample,Treatment,Time,group,path,fasta,geTMM,description,completion_percent)%>%
  distinct()%>%
  group_by(sample,Treatment,Time,group,path,fasta,geTMM)%>%
  summarise(completion=sum(completion_percent))%>%
  mutate(completion_2=ifelse(path=="carABCDE",ifelse(completion>=0.4,1,0),completion))%>%
  filter(completion_2>=0.2)%>%
  left_join(.,mags,by=c("fasta"="MAG"))%>%
  separate(`GTDB v2.3.0 r214`,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  mutate(genus=paste(d,p,c,o,f,g,sep=";"))%>%
  mutate(mag=paste(p,f,g,fasta,sep=";"))%>%
  ungroup()%>%
  select(mag,p,group,fasta)%>%
  distinct()%>%
  group_by(mag,p,group)%>%
  summarise(n_MAGs=n())%>%
  filter(is.na(p)==0)%>%
  spread(key=group,value=n_MAGs,fill=0)%>%
  mutate(count=carABCDE+CHI+CDH+FCR+hpaAB+hpaDEFGH+mhpABCD+PGR+pgthAB+PHY)%>%
  filter(count>=2)%>%select(-count)%>%
  gather(-mag,-p,key="group",value="n_MAGs")%>%
  filter(n_MAGs>0)%>%
  ggplot()+
  geom_point(aes(x=mag,y=ordered(group,levels=c("FCR","PHY","CHI","CDH","PGR","pgthAB","hpaAB","hpaDEFGH","mhpABCD","carABCDE")),size=n_MAGs>0,color=p))+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8),
        axis.text.y = element_text(size=5))
