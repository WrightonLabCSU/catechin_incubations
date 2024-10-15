library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)



###############
## plot fig. 6A
###############
hyd_type=read_excel("SOM_files/Supplementary_Data_8.xlsx",sheet="hydrogenase_key")
hyd_metat=read_excel("SOM_files/Supplementary_Data_8.xlsx",sheet="hydrogenase_expression")

# read in metadata
metadata=read_excel("SOM_files/Supplementary_Data_3.xlsx",sheet="Metatranscriptomes")
metadata=metadata%>%
  select(Sample,Treatment,Time)

# read in MAG taxonomy
mags=read_excel("SOM_files/Supplementary_Data_2.xlsx",sheet="MAG Database")
mags=mags%>%
  select(MAG,`GTDB v2.3.0 r214`)

hyd_metat=hyd_metat%>%
  gather(-gene,-hydrogenase_group,-MAG,key="sample",value="geTMM")%>%
  left_join(.,metadata,by=c("sample"="Sample"))

hyd_metat%>%
  left_join(.,mags,by=c("MAG"))%>%
  separate(`GTDB v2.3.0 r214`,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  filter(p!="p__Halobacteriota"|p!="p__Methanobacteriota"|p!="p__Thermoplasmatota")%>%
  left_join(.,hyd_type,by=c("hydrogenase_group"))%>%
  filter(is.na(direction)==0)%>%
  filter(direction!="sensing")%>%
  group_by(direction,sample,Treatment,Time)%>%
  summarise(sum=sum(geTMM))%>%
  ggplot()+
  geom_point(aes(x=Time,y=sum,color=Treatment))+
  geom_smooth(aes(x=Time,y=sum,color=Treatment,fill=Treatment))+
  facet_wrap(~direction,scales="free_y")+
  scale_color_manual(values=c("#fe9929","#7bccc4"))+
  scale_fill_manual(values=c("#fe9929","#7bccc4"))+
  theme_classic()

#stats on hydrogenase expression
hyd_data=hyd_metat%>%
  left_join(.,mags,by=c("MAG"))%>%
  separate(`GTDB v2.3.0 r214`,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  filter(p!="p__Halobacteriota"|p!="p__Methanobacteriota"|p!="p__Thermoplasmatota")%>%
  left_join(.,hyd_type,by=c("hydrogenase_group"))%>%
  filter(is.na(direction)==0)%>%
  filter(direction!="sensing")%>%
  group_by(direction,sample,Treatment,Time)%>%
  summarise(sum=sum(geTMM))

hyd_data%>%
  filter(direction=="bidirectional")%>%
  ungroup()%>%
  select(-sample)%>%
  filter(Time>0)%>%
  nest(data = c(Treatment,sum)) %>% 
  mutate(kruskal.raw = map(data, ~ kruskal.test(sum~Treatment, .,)),
         kruskal = map(kruskal.raw, broom::tidy)) %>%
  select(-data) %>%
  unnest(kruskal)%>%
  select(Time,p.value)

hyd_data%>%
  filter(direction=="bifurcating")%>%
  ungroup()%>%
  select(-sample)%>%
  filter(Time>0)%>%
  nest(data = c(Treatment,sum)) %>% 
  mutate(kruskal.raw = map(data, ~ kruskal.test(sum~Treatment, .,)),
         kruskal = map(kruskal.raw, broom::tidy)) %>%
  select(-data) %>%
  unnest(kruskal)%>%
  select(Time,p.value)

hyd_data%>%
  filter(direction=="producing")%>%
  ungroup()%>%
  select(-sample)%>%
  filter(Time>0)%>%
  nest(data = c(Treatment,sum)) %>% 
  mutate(kruskal.raw = map(data, ~ kruskal.test(sum~Treatment, .,)),
         kruskal = map(kruskal.raw, broom::tidy)) %>%
  select(-data) %>%
  unnest(kruskal)%>%
  select(Time,p.value)

hyd_data%>%
  filter(direction=="uptake")%>%
  ungroup()%>%
  select(-sample)%>%
  filter(Time>0)%>%
  nest(data = c(Treatment,sum)) %>% 
  mutate(kruskal.raw = map(data, ~ kruskal.test(sum~Treatment, .,)),
         kruskal = map(kruskal.raw, broom::tidy)) %>%
  select(-data) %>%
  unnest(kruskal)%>%
  select(Time,p.value)


###############
## plot fig. 6B
###############
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

mat1=figure_genes%>%
  left_join(k3,by="annotation")%>%
  left_join(.,metadata,by=c("sample"="Sample"))%>%
  filter(filter!="yes")%>%
  ungroup()%>%
  select(sample,Treatment,Time,group,path,completion_percent,fasta,geTMM,description)%>%
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
  select(sample,mag,p,group,fasta,geTMM)%>%
  group_by(sample,mag,group,fasta)%>%
  summarise(sum=sum(geTMM))%>%
  distinct()%>%
  full_join(.,hyd_metat,by=c("fasta"="MAG","sample"),relationship = "many-to-many")%>%
  filter(hydrogenase_group!="fefe_C1",hydrogenase_group!="fefe_C2",hydrogenase_group!="fefe_C3",hydrogenase_group!="nife_2c",hydrogenase_group!="nife_2b")%>%
  filter(is.na(group)==0)%>%
  filter(is.na(hydrogenase_group)==0)%>%
  group_by(sample,group,hydrogenase_group,Treatment,Time)%>%
  summarise(group_sum=sum(sum),hyd_sum=sum(geTMM))%>%
  ungroup()%>%
  select(sample,group,group_sum,hydrogenase_group,hyd_sum)%>%
  spread(key=group,value=group_sum,fill=0)%>%
  spread(key=hydrogenase_group,value=hyd_sum,fill=0)%>%
  gather(-sample,key="id",value="sum")%>%
  group_by(sample,id)%>%
  summarise(total=sum(sum))%>%
  spread(key=id,value=total)

mat1=as.data.frame(mat1)
row.names(mat1)=mat1[,1]
mat1=mat1[,-1]
mat1=as.matrix(mat1)

library(corrplot)
mat2=cor(mat1,method="spearman")
p_val=cor.mtest(mat1)
corrplot(mat2,type="upper",p.mat=p_val$p,sig.level=0.05, insig = "blank", order="hclust")
