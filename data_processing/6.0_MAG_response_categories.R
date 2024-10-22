library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)


# get metadata
metadata=read_excel("Supplementary_Data_3.xlsx",sheet="Metatranscriptomes")
metadata=metadata%>%
  select(Sample,Treatment,Time)

# read in gene metaT data (access from zenodo repository: https://doi.org/10.5281/zenodo.13937409)
getmm=read.delim("geTMMs_ge5_1X_REV_26samples.txt",header=TRUE,sep="\t")

# read in MAG database annotations (access from zenodo repository: https://doi.org/10.5281/zenodo.13936221)
annotations=read.delim("reactorEMERGE_annotations.tsv",header=TRUE,sep="\t")

#make tidy
tidy_getmm=getmm%>%
  gather(-gene,key="sample",value="geTMM")%>%
  filter(geTMM>0)

#add annotations to genes
tidy_getmm_annos=left_join(tidy_getmm,annotations,by=c("gene"="X"))

tidy_getmm_wMags=tidy_getmm_annos%>%
  select(gene,fasta,sample,geTMM)%>%
  left_join(.,metadata,by=c("sample"="Sample"))

tidy_mag_geTMM=read_excel("Supplementary_Data_3.xlsx",sheet="MAG_metatranscriptome")%>%
  gather(-MAG,key="sample",value="abund")%>%
  filter(abund>0)

genes_mags22=tidy_mag_geTMM%>%
  select(MAG,sample)%>%
  left_join(.,tidy_getmm_wMags,by=c("MAG"="fasta","sample"))

# get genes per timepoint
f7=genes_mags22%>%
  filter(Time==7)%>%
  group_by(Treatment,MAG,gene,Time)%>%
  summarise(ave=mean(geTMM))%>%
  spread(key=Treatment,value=ave)%>%
  mutate(catechin_pa=ifelse(is.na(catechin)==1,0,1),
                                       unamended_pa=ifelse(is.na(unamended)==1,0,1))%>%
  mutate(status=ifelse(catechin_pa==1,ifelse(unamended_pa==1,"both","cat"),"un"))%>%
  group_by(MAG,Time,status)%>%summarise(n=n())

f14=genes_mags22%>%
  filter(Time==14)%>%
  group_by(Treatment,MAG,gene,Time)%>%
  summarise(ave=mean(geTMM))%>%
  spread(key=Treatment,value=ave)%>%
  mutate(catechin_pa=ifelse(is.na(catechin)==1,0,1),
         unamended_pa=ifelse(is.na(unamended)==1,0,1))%>%
  mutate(status=ifelse(catechin_pa==1,ifelse(unamended_pa==1,"both","cat"),"un"))%>%
  group_by(MAG,Time,status)%>%summarise(n=n())

f21=genes_mags22%>%
  filter(Time==21)%>%
  group_by(Treatment,MAG,gene,Time)%>%
  summarise(ave=mean(geTMM))%>%
  spread(key=Treatment,value=ave)%>%
  mutate(catechin_pa=ifelse(is.na(catechin)==1,0,1),
         unamended_pa=ifelse(is.na(unamended)==1,0,1))%>%
  mutate(status=ifelse(catechin_pa==1,ifelse(unamended_pa==1,"both","cat"),"un"))%>%
  group_by(MAG,Time,status)%>%summarise(n=n())

f35=genes_mags22%>%
  filter(Time==35)%>%
  group_by(Treatment,MAG,gene,Time)%>%
  summarise(ave=mean(geTMM))%>%
  spread(key=Treatment,value=ave)%>%
  mutate(catechin_pa=ifelse(is.na(catechin)==1,0,1),
         unamended_pa=ifelse(is.na(unamended)==1,0,1))%>%
  mutate(status=ifelse(catechin_pa==1,ifelse(unamended_pa==1,"both","cat"),"un"))%>%
  group_by(MAG,Time,status)%>%summarise(n=n())


mag_gene_stats=rbind(f7,f14,f21,f35)
remove(f7)
remove(f14)
remove(f21)
remove(f35)

mag_gene_stats_sum=mag_gene_stats%>%
  group_by(MAG,Time)%>%summarise(total=sum(n))
mag_gene_stats=mag_gene_stats%>%
  left_join(.,mag_gene_stats_sum,by=c("MAG","Time"))%>%
  mutate(per=n/total)

expression_pattern=mag_gene_stats%>%
  select(-n,-total)%>%
  spread(key=status,value=per)%>%
  mutate(both=ifelse(is.na(both)==1,0,both),un=ifelse(is.na(un)==1,0,un),cat=ifelse(is.na(cat)==1,0,cat))%>%
  gather(both,cat,un,key="status",value="per")


# assigning categories
q=mag_gene_stats%>%
  select(-n,-total)%>%
  spread(key=status,value=per,fill=0)%>%
  mutate(x_axis=un-cat,y_axis=both)
quantile(q$x_axis)
ecdf(q$x_axis)
cat50=-0.25632277 #25%
un50=0.24571994 #75%
quantile(q$y_axis)
both50=0.2643335 #50%

mag_categories=mag_gene_stats%>%
  select(-n,-total)%>%
  spread(key=status,value=per)%>%
  mutate(both=ifelse(is.na(both)==1,0,both),
         un=ifelse(is.na(un)==1,0,un),
         cat=ifelse(is.na(cat)==1,0,cat))%>%
  mutate(x_axis=un-cat,y_axis=both)%>%
  mutate(category=ifelse(y_axis>=both50,
                         ifelse(x_axis>un50,"lost_fx",ifelse(x_axis<cat50,"gain_fx","resistant")),
                         ifelse(x_axis>un50,"un_only",ifelse(x_axis<cat50,"cat_only","responsive"))))


########
## plots for Fig S2
v=mag_categories%>%
  ggplot()+
  geom_histogram(aes(x=y_axis),bins=100)+
  geom_vline(xintercept=0.264)+
  theme_classic()

x=mag_categories%>%
  ggplot()+
  geom_histogram(aes(x=x_axis),bins=100)+
  geom_vline(xintercept=0.245)+
  geom_vline(xintercept=-0.256)+
  theme_classic()

plot_grid(v,x)


mag_categories%>%
  ggplot()+
  geom_point(aes(x=x_axis,y=y_axis,color=ordered(category,levels=c("cat_only","gain_fx","responsive","lost_fx","un_only","resistant"))),alpha=0.4)+
  #facet_wrap(~time)+
  theme_classic()+
  #scale_color_manual(values=c("#7bccc4","#fe9929"))+
  geom_hline(aes(yintercept=both50))+
  geom_vline(aes(xintercept=un50))+
  geom_vline(aes(xintercept=cat50))+
  xlab("%Genes Unique to Unamended-%Genes Unique to Catechin")+ylab("%Genes in Both")+
  scale_color_manual(values=c("#fe9929","#fdbf6f","#fb8072","#a8ddb5","#7bccc4","#7570b3"))+
  facet_wrap(~Time)
