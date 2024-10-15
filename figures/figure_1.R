library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

gas_data=read_excel("Supplementary_Data_1.xlsx")

# plot gas curves
gas_data%>%
  gather(`CO2 (ppm)`,`CH4 (ppm)`,key="gas",value="ppm")%>%
  filter(ppm!="NA")%>%
  mutate(ppm=as.numeric(ppm))%>%
  ggplot()+
  geom_point(aes(x=day,y=ppm,color=treatment))+
  geom_smooth(aes(x=day,y=ppm,color=treatment,fill=treatment))+
  facet_wrap(~gas,scales="free_y")+
  theme_classic()+
  scale_color_manual(values=c("#f8992c","#7accc4"))+
  scale_fill_manual(values=c("#f8992c","#7accc4"))


# kruskal-wallis test
gas_data%>%
  gather(`CO2 (ppm)`,`CH4 (ppm)`,key="gas",value="ppm")%>%
  filter(ppm!="NA")%>%
  mutate(ppm=as.numeric(ppm))%>%
  select(-Sample,-IGSN,-other_name,-site,-`peat (g)`,-pH)%>%
  nest(data = c(treatment,ppm)) %>% 
  mutate(kruskal_raw = map(data, ~ kruskal.test(ppm~treatment, .)),
         kruskal = map(kruskal_raw, broom::tidy)) %>%
  select(-data) %>%
  unnest(kruskal)%>%
  select(gas,day,p.value)
