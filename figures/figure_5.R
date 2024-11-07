library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(BBmisc)
library(pheatmap)

# read in MAG-level metaT table
mag_geTMM=read_excel("Supplementary_Data_3.xlsx",sheet="MAG_metatranscriptome")
tidy_mag_geTMM=mag_geTMM%>%
  gather(-MAG,key="sample",value="mag_geTMM")

# read in MAG taxonomy
mags=read_excel("Supplementary_Data_2.xlsx",sheet="MAG Database")
mags=mags%>%
  select(MAG,`GTDB v2.3.0 r214`)

# read in metadata
metadata=read_excel("Supplementary_Data_3.xlsx",sheet="Metatranscriptomes")
metadata=metadata%>%
  select(Sample,Treatment,Time)

#####
## plot figure 5A
#####

mgen_max=tidy_mag_geTMM%>%
  left_join(.,metadata,by=c(c("sample"="Sample")))%>%
  left_join(.,mags,by="MAG")%>%
  separate(`GTDB v2.3.0 r214`,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  filter(p=="p__Halobacteriota"|p=="p__Methanobacteriota"|p=="p__Thermoplasmatota")%>%
  mutate(short_name=paste(f,g,sep=";"))%>%
  left_join(.,mgen_sum,by="sample")%>%
  group_by(short_name,MAG,sample)%>%
  summarise(sum=sum(mag_geTMM))%>%
  ungroup()%>%
  group_by(short_name,MAG)%>%
  summarise(max=max(sum))

wide_mgen_sum=tidy_mag_geTMM%>%
  left_join(.,metadata,by=c(c("sample"="Sample")))%>%
  left_join(.,mags,by="MAG")%>%
  separate(`GTDB v2.3.0 r214`,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  filter(p=="p__Halobacteriota"|p=="p__Methanobacteriota"|p=="p__Thermoplasmatota")%>%
  mutate(short_name=paste(f,g,sep=";"))%>%
  left_join(.,mgen_max,by=c("MAG","short_name"))%>%
  group_by(short_name,MAG,sample,Treatment,Time,max)%>%
  summarise(sum=sum(mag_geTMM),norm=sum(mag_geTMM)/max)%>%
  mutate(label=paste(short_name,MAG,sep="_"))%>%
  ungroup()%>%
  select(label,sample,norm)%>%
  spread(key=sample,value=norm,fill=0)
wide_mgen_sum=as.data.frame(wide_mgen_sum)
row.names(wide_mgen_sum)=wide_mgen_sum[,1]
wide_mgen_sum=wide_mgen_sum[,-1]
wide_mgen_sum=wide_mgen_sum[,c("STM_0716_E_M_E002","STM_0716_E_M_E003","STM_0716_E_M_E004","STM_0716_E_M_E121","STM_0716_E_M_E122","STM_0716_E_M_E123","STM_0716_E_M_E025","STM_0716_E_M_E027","STM_0716_E_M_E050","STM_0716_E_M_E051","STM_0716_E_M_E052","STM_0716_E_M_E062","STM_0716_E_M_E063","STM_0716_E_M_E064","STM_0716_E_M_E129","STM_0716_E_M_E130","STM_0716_E_M_E131","STM_0716_E_M_E033","STM_0716_E_M_E034","STM_0716_E_M_E035","STM_0716_E_M_E058","STM_0716_E_M_E059","STM_0716_E_M_E060","STM_0716_E_M_E070","STM_0716_E_M_E071","STM_0716_E_M_E072")]


pheatmap(wide_mgen_sum,scale="none",cluster_rows = FALSE,cluster_cols=FALSE,color = colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100))


# response category heatmap
mag_categories=read_excel("SOM_files/Supplementary_Data_4.xlsx",sheet="Classifications")

mag_categories%>%
  separate(GTDB,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  filter(p=="p__Halobacteriota"|p=="p__Methanobacteriota"|p=="p__Thermoplasmatota")%>%
  mutate(short_name=paste(f,g,sep=";"))%>%
  mutate(label=paste(short_name,MAG,sep="_"))%>%
  select(label,time,classification)%>%
  ggplot()+
  geom_tile(aes(x=time,y=label,fill=classification))+
  theme_classic()+
  scale_fill_manual(values=c("#fe9929","#fdbf6f","#a8ddb5","#7570b3","#fb8072","#7bccc4"))



#####
## plot figure 5B
#####

mgen_genes=read_excel("SOM_files/Supplementary_Data_7.xlsx",sheet="Summarized_Methanogen_MetaT")
mgen_genes=mgen_genes%>%
  gather(-MAG,-GTDB,-MAG_Potential_Curation,key="set",value="substrate")%>%
  separate(set,into=c("Treatment","Time"),sep="_")%>%
  mutate(Time=as.double(Time))%>%
  select(-GTDB)

tidy_mag_geTMM%>%
  left_join(.,metadata,by=c(c("sample"="Sample")))%>%
  left_join(.,mags,by="MAG")%>%
  separate(`GTDB v2.3.0 r214`,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  filter(p=="p__Halobacteriota"|p=="p__Methanobacteriota"|p=="p__Thermoplasmatota")%>%
  left_join(.,mgen_genes,by=c("MAG","Treatment","Time"),relationship = "many-to-many")%>%
  filter(substrate!="general_methanogenesis")%>%
  group_by(substrate,sample,Time,Treatment)%>%
  summarise(sum=sum(mag_geTMM))%>%
  filter(substrate!="x")%>%
  ggplot()+
  geom_point(aes(x=Time,y=sum,color=Treatment))+
  geom_smooth(aes(x=Time,y=sum,color=Treatment,fill=Treatment))+
  theme_classic()+
  facet_wrap(~substrate,scales="free_y")+
  scale_color_manual(values=c("#f8992c","#7accc4"))+
  scale_fill_manual(values=c("#f8992c","#7accc4"))

# stats for substrate use
substrate_data=tidy_mag_geTMM%>%
  left_join(.,metadata,by=c(c("sample"="Sample")))%>%
  left_join(.,mags,by="MAG")%>%
  separate(`GTDB v2.3.0 r214`,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  filter(p=="p__Halobacteriota"|p=="p__Methanobacteriota"|p=="p__Thermoplasmatota")%>%
  left_join(.,mgen_genes,by=c("MAG","Treatment","Time"),relationship = "many-to-many")%>%
  filter(substrate!="general_methanogenesis")%>%
  group_by(substrate,sample,Time,Treatment)%>%
  summarise(sum=sum(mag_geTMM))

substrate_data%>%
  filter(substrate=="Acetate")%>%
  ungroup()%>%
  select(-sample)%>%
  filter(Time>0)%>%
  nest(data = c(Treatment,sum)) %>% 
  mutate(kruskal.raw = map(data, ~ kruskal.test(sum~Treatment, .,)),
         kruskal = map(kruskal.raw, broom::tidy)) %>%
  select(-data) %>%
  unnest(kruskal)%>%
  select(Time,p.value)

substrate_data%>%
  filter(substrate=="formate")%>%
  ungroup()%>%
  select(-sample)%>%
  filter(Time>0)%>%
  nest(data = c(Treatment,sum)) %>% 
  mutate(kruskal.raw = map(data, ~ kruskal.test(sum~Treatment, .,)),
         kruskal = map(kruskal.raw, broom::tidy)) %>%
  select(-data) %>%
  unnest(kruskal)%>%
  select(Time,p.value)

substrate_data%>%
  filter(substrate=="Hydrogen/CO2")%>%
  ungroup()%>%
  select(-sample)%>%
  filter(Time>0)%>%
  nest(data = c(Treatment,sum)) %>% 
  mutate(kruskal.raw = map(data, ~ kruskal.test(sum~Treatment, .,)),
         kruskal = map(kruskal.raw, broom::tidy)) %>%
  select(-data) %>%
  unnest(kruskal)%>%
  select(Time,p.value)

substrate_data%>%
  filter(substrate=="MeOH")%>%
  ungroup()%>%
  select(-sample)%>%
  filter(Time>0)%>%
  nest(data = c(Treatment,sum)) %>% 
  mutate(kruskal.raw = map(data, ~ kruskal.test(sum~Treatment, .,)),
         kruskal = map(kruskal.raw, broom::tidy)) %>%
  select(-data) %>%
  unnest(kruskal)%>%
  select(Time,p.value)

substrate_data%>%
  filter(substrate=="MethylS")%>%
  ungroup()%>%
  select(-sample)%>%
  filter(Time>0)%>%
  nest(data = c(Treatment,sum)) %>% 
  mutate(kruskal.raw = map(data, ~ kruskal.test(sum~Treatment, .,)),
         kruskal = map(kruskal.raw, broom::tidy)) %>%
  select(-data) %>%
  unnest(kruskal)%>%
  select(Time,p.value)



#####
## plot figure 5C
#####
metadata=read_excel("SOM_files/Supplementary_Data_5.xlsx",sheet="Metadata")

nmr=read_excel("SOM_files/Supplementary_Data_5.xlsx",sheet="NMR")
nmr=nmr%>%
  filter(metabolite=="Acetate"|metabolite=="Formate"|metabolite=="Methanol")
nmr=nmr[,c(1,11:37)]
nmr=nmr%>%
  gather(-metabolite,key="sample",value="value")%>%
  mutate(value=ifelse(is.na(value)==1,0,value))

nmr%>%
  left_join(.,metadata,by="sample")%>%
  ggplot()+
  geom_point(aes(x=day,y=value,color=treatment,fill=treatment))+
  geom_smooth(aes(x=day,y=value,color=treatment,fill=treatment))+
  facet_wrap(~metabolite,scales="free_y")+
  scale_color_manual(values=c("#fe9929","#7bccc4"))+
  scale_fill_manual(values=c("#fe9929","#7bccc4"))+
  theme_classic()


# stats for substrate use
substrate_data=nmr%>%
  left_join(.,metadata,by="sample")

substrate_data%>%
  filter(metabolite=="Acetate")%>%
  ungroup()%>%
  select(-sample)%>%
  filter(day>0)%>%
  nest(data = c(treatment,value)) %>% 
  mutate(kruskal.raw = map(data, ~ kruskal.test(value~treatment, .,)),
         kruskal = map(kruskal.raw, broom::tidy)) %>%
  select(-data) %>%
  unnest(kruskal)%>%
  select(day,p.value)

substrate_data%>%
  filter(metabolite=="Methanol")%>%
  ungroup()%>%
  select(-sample)%>%
  filter(day>0)%>%
  nest(data = c(treatment,value)) %>% 
  mutate(kruskal.raw = map(data, ~ kruskal.test(value~treatment, .,)),
         kruskal = map(kruskal.raw, broom::tidy)) %>%
  select(-data) %>%
  unnest(kruskal)%>%
  select(day,p.value)

substrate_data%>%
  filter(metabolite=="Formate")%>%
  ungroup()%>%
  select(-sample)%>%
  filter(day>0)%>%
  nest(data = c(treatment,value)) %>% 
  mutate(kruskal.raw = map(data, ~ kruskal.test(value~treatment, .,)),
         kruskal = map(kruskal.raw, broom::tidy)) %>%
  select(-data) %>%
  unnest(kruskal)%>%
  select(day,p.value)
