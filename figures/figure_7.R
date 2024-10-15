library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)



mag_categories=read_excel("Supplementary_Data_4.xlsx",sheet="Classifications")
curated_annos=read_excel("Supplementary_Data_3.xlsx",sheet="curated_annotations")
product=read_excel("Supplementary_Data_3.xlsx",sheet="c_cycle_pathway_genes")


# isolate the non-methanogen genera that are sensitive or lost functions at days 21 and 35. must be at least 2 MAGs in the genus, with >=50% total active mags
tot_mags=mag_categories%>%
  select(MAG,GTDB)%>%
  distinct()%>%
  separate(GTDB,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  mutate(short_name=paste(p,f,g,sep=";"))%>%
  group_by(short_name)%>%
  summarise(total_active=n())

others=mag_categories%>%
  filter(time>14)%>%
  filter(classification=="lost_fx"|classification=="un_only")%>%
  select(MAG,GTDB)%>%
  distinct()%>%
  separate(GTDB,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  filter(p!="p__Halobacteriota")%>%
  filter(p!="p__Methanobacteriota")%>%
  filter(p!="p__Thermoplasmatota")%>%
  mutate(short_name=paste(p,f,g,sep=";"))%>%
  group_by(short_name)%>%
  summarise(n=n())%>%
  filter(n>=2)%>%
  left_join(.,tot_mags,by="short_name")%>%
  mutate(per_impacted=n/total_active)%>%
  filter(per_impacted>=0.5)

others=mag_categories%>%
  filter(time>14)%>%
  filter(classification=="lost_fx"|classification=="un_only")%>%
  select(MAG,GTDB)%>%
  distinct()%>%
  separate(GTDB,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  filter(p!="p__Halobacteriota")%>%
  filter(p!="p__Methanobacteriota")%>%
  filter(p!="p__Thermoplasmatota")%>%
  mutate(short_name=paste(p,f,g,sep=";"))%>%
  select(MAG,short_name)%>%
  right_join(.,others,by="short_name")




# read in gene metaT data (access from zenodo repository: 10.5281/zenodo.13375972)
getmm=read.delim("geTMMs_ge5_1X_REV_26samples.txt",header=TRUE,sep="\t")

# read in MAG database annotations (access from zenodo repository)
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

# make at MAG level
tidy_mag_geTMM=tidy_getmm_annos%>%
  group_by(sample,fasta)%>%
  summarise(n_genes=n(),mag_geTMM=sum(geTMM))
# deciding 22 because it is the 50% value
tidy_mag_geTMM=tidy_mag_geTMM%>%
  filter(n_genes>=22)
#get only genes from MAGs >=22 genes
tidy_getmm_ge22=tidy_mag_geTMM%>%
  left_join(.,tidy_getmm_annos,by=c("sample","fasta"))
tidy_getmm_ge22=tidy_getmm_ge22%>%
  select(sample,fasta,gene,geTMM)


wide_gen_fx=fx_ge2%>%
  mutate(best_id=ifelse(substr(best_id,1,2)=="GH","GH",
                        ifelse(substr(best_id,1,2)=="PL","PL",
                               ifelse(substr(best_id,1,2)=="CE","CE",
                                      ifelse(substr(best_id,1,2)=="AA","AA",best_id)))))%>%
  left_join(.,product,by="best_id",relationship = "many-to-many")%>%
  filter(is.na(category)==0)%>%
  ungroup()%>%
  group_by(short_name,best_id)%>%
  summarise(n=n())%>%
  mutate(n=ifelse(n>0,1,0))%>%
  spread(key=best_id,value=n,fill=0)%>%
  mutate(methanotroph=MMO_A+MMO_B+MMO_C)%>%
  mutate(sulfate=K11180+K11181)%>%
  mutate(nitrate=K00362+K00368+K04561)%>%
  mutate(prop=(K01847+K01849+K01848+K05606+K03416+K11264)/3)%>%
  mutate(prop=ifelse(prop>=0.5,1,0))%>%
  mutate(butyrate=K00634+K00929+K01034)%>%
  mutate(v1_acetate=(K00625+K13788+K00925)/3)%>%
  mutate(v1_acetate=ifelse(v1_acetate>0.5,1,0))%>%
  mutate(v2_acetate=(K01905+K22224+K24012)/3)%>%
  mutate(v1_acetate=ifelse(v2_acetate>0.5,1,0))%>%
  mutate(v3_acetate=(K01895))%>%
  mutate(v4_acetate=K18118)%>%
  mutate(acetate=v1_acetate+v2_acetate+v3_acetate+v4_acetate)%>%
  mutate(fructose=K00844+K00847)%>%
  mutate(fucose=K00879+K01628+K01818)%>%
  mutate(galactose=(K00849+K00965+K01784+K01785)/4)%>%
  mutate(galactose=ifelse(galactose>=0.5,1,0))%>%
  mutate(galacturonic=K00041+K00874+K01625+K01685+K16850)%>%
  mutate(mannose=(K00844+K01809)/2)%>%
  mutate(mannose=ifelse(mannose>=0.5,1,0))%>%
  mutate(sucrose=K00690+K01187)%>%
  mutate(v1_xylose=K00854+K01805)%>%
  mutate(v2_xylose=(K00008+K00854+K05351)/3)%>%
  mutate(v2_xylose=ifelse(v2_xylose>=0.5,1,0))%>%
  mutate(v3_xylose=K14273)%>%
  mutate(xylose=v1_xylose+v2_xylose+v3_xylose)%>%
  mutate(wl=(K00198+K05299+K00297+K00194+K00197+K15023+K14138+K15022+K25124+K01938+K01491+K01500)/12)%>%
  mutate(wl=ifelse(wl>=0.25,1,0))%>%
  mutate(syntroph=ifelse(short_name=="p__Desulfobacterota_G;f__Syntrophorhabdaceae;g__PNOF01",1,
                         ifelse(short_name=="p__Desulfobacterota_G;f__WCHB1-27;g__",1,
                                ifelse(short_name=="p__Desulfobacterota;f__Fen-1087;g__Fen-1087",1,
                                       ifelse(short_name=="p__Desulfobacterota;f__Smithellaceae;g__FEN-1160",1,
                                              ifelse(short_name=="p__Desulfobacterota;f__Smithellaceae;g__Smithella",1,
                                                     ifelse(short_name=="p__Desulfobacterota;f__UBA2185;g__Fen-1135",1,0)))))))%>%
  gather(-short_name,key="path",value="pres")%>%
  filter(substr(path,1,2)!="K0")%>%
  filter(substr(path,1,2)!="K1")%>%
  filter(substr(path,1,2)!="K2")%>%
  filter(substr(path,1,2)!="MM")%>%
  filter(substr(path,1,1)!="v")%>%
  filter(substr(path,1,1)!="A")%>%
  filter(substr(path,1,1)!="P")%>%
  filter(substr(path,1,1)!="C")%>%
  mutate(pres=ifelse(pres>0,1,0))%>%
  spread(key=path,value=pres,fill=0)
wide_gen_fx=as.data.frame(wide_gen_fx)
row.names(wide_gen_fx)=wide_gen_fx[,1]
wide_gen_fx=wide_gen_fx[,-1]
library(pheatmap)
pheatmap(wide_gen_fx,cluster_rows = TRUE,cluster_cols = FALSE)

