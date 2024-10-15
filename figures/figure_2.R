library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)


#######################
# figure 2A
#######################
bray_distances=read_excel("Supplementary_Data_4.xlsx",sheet="Distances")

ggplot(bray_distances)+
  geom_boxplot(aes(x=ordered(as.character(day),levels=c("7","14","21","35")),y=distance,color=ordered(data,levels=c("16S","metaT_MAG","metaT_gene","metaT_fx","lcms"))))+
  geom_jitter(aes(x=ordered(as.character(day),levels=c("7","14","21","35")),y=distance,color=ordered(data,levels=c("16S","metaT_MAG","metaT_gene","metaT_fx","lcms"))),width = 0.3)+
  theme_classic()+
  xlab("Day")+ylab("Bray Distance")

#####
## testing for differences in distance at each timepoint
#####
# day 7
bray_stat=bray%>%filter(day==7)
# run ANOVA and store results
res_aov <- aov(distance ~ data,data = bray_stat)
# check for normality using the residuals
library(car)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals, id = FALSE)
# need the qqplot to be a straight line, normal distribution in histogram
shapiro.test(res_aov$residuals) # p<0.05 = not normal
#assumption 2:  are the variance equal?
leveneTest(distance ~ data,data = bray_stat) # p<0.05 = not equal variance

kruskal.test(distance~data, data=bray_stat) #sig p-value = 2.154e-06
# if no, stop here
# if yes, run Dunn's test
d=dunnTest(distance~data,
           data=bray_stat,
           method="bh")
# > #           Comparison     Z          P.unadj      P.adj
#   > # 1             16S - lcms -4.7428361 2.107467e-06 0.0000126448
#   > # 2       16S - metaT_gene -4.1164238 3.847965e-05 0.0001154389
#   > # 3      lcms - metaT_gene  0.6264123 5.310445e-01 0.5310445320
#   > # 4        16S - metaT_MAG -1.6555183 9.781941e-02 0.1173832964
#   > # 5       lcms - metaT_MAG  3.0873179 2.019715e-03 0.0040394301
#   > # 6 metaT_gene - metaT_MAG  2.4609055 1.385868e-02 0.0207880275
p7=c(d[["res"]][["P.unadj"]])

# day 14
bray_stat=bray%>%filter(day==14)
# run ANOVA and store results
res_aov <- aov(distance ~ data,data = bray_stat)
# check for normality using the residuals
library(car)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals, id = FALSE)
# need the qqplot to be a straight line, normal distribution in histogram
shapiro.test(res_aov$residuals) # p>0.05 = normal
#assumption 2:  are the variance equal?
leveneTest(distance ~ data,data = bray_stat) # p<0.05 = not equal variance

kruskal.test(distance~data, data=bray_stat) #sig p-value = 6.147e-05
# if no, stop here
# if yes, run Dunn's test
d=dunnTest(distance~data,
           data=bray_stat,
           method="bh")
# Comparison           Z      P.unadj       P.adj
# 1              16S - lcms -3.62424270 0.0002898094 0.002898094
# 2          16S - metaT_fx -0.03001501 0.9760550813 0.976055081
# 3         lcms - metaT_fx  3.21160620 0.0013199517 0.003299879
# 4        16S - metaT_gene -3.60180135 0.0003160198 0.001580099
# 5       lcms - metaT_gene -0.36018014 0.7187124291 0.798569366
# 6   metaT_fx - metaT_gene -3.26057991 0.0011118463 0.003706154
# 7         16S - metaT_MAG -0.69034526 0.4899770925 0.699967275
# 8        lcms - metaT_MAG  2.55127596 0.0107329306 0.017888218
# 9    metaT_fx - metaT_MAG -0.60279629 0.5466442173 0.683305272
# 10 metaT_gene - metaT_MAG  2.65778363 0.0078656374 0.015731275
p14=c(d[["res"]][["P.unadj"]])

# day 21
bray_stat=bray%>%filter(day==21)
# run ANOVA and store results
res_aov <- aov(distance ~ data,data = bray_stat)
# check for normality using the residuals
library(car)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals, id = FALSE)
# need the qqplot to be a straight line, normal distribution in histogram
shapiro.test(res_aov$residuals) # p>0.05 = normal
#assumption 2:  are the variance equal?
leveneTest(distance ~ data,data = bray_stat) # p<0.05 = not equal variance

kruskal.test(distance~data, data=bray_stat) #sig p-value = 8.82e-07
# if no, stop here
# if yes, run Dunn's test
d=dunnTest(distance~data,
           data=bray_stat,
           method="bh")
# Comparison           Z      P.unadj        P.adj
# 1              16S - lcms -4.86338318 1.153961e-06 5.769807e-06
# 2          16S - metaT_fx  0.03589213 9.713684e-01 9.713684e-01
# 3         lcms - metaT_fx  4.89927531 9.619078e-07 9.619078e-06
# 4        16S - metaT_gene -3.48153630 4.985462e-04 1.246366e-03
# 5       lcms - metaT_gene  1.38184688 1.670187e-01 1.855764e-01
# 6   metaT_fx - metaT_gene -3.51742843 4.357498e-04 1.452499e-03
# 7         16S - metaT_MAG -1.83049847 6.717544e-02 9.596491e-02
# 8        lcms - metaT_MAG  3.03288471 2.422281e-03 4.844562e-03
# 9    metaT_fx - metaT_MAG -1.86639059 6.198674e-02 1.033112e-01
# 10 metaT_gene - metaT_MAG  1.65103783 9.873085e-02 1.234136e-01
p21=c(d[["res"]][["P.unadj"]])

# day 35
bray_stat=bray%>%filter(day==35)
# run ANOVA and store results
res_aov <- aov(distance ~ data,data = bray_stat)
# check for normality using the residuals
library(car)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals, id = FALSE)
# need the qqplot to be a straight line, normal distribution in histogram
shapiro.test(res_aov$residuals) # p>0.05 = normal
#assumption 2:  are the variance equal?
leveneTest(distance ~ data,data = bray_stat) # p>0.05 = equal variance

kruskal.test(distance~data, data=bray_stat) #sig p-value = 1.251e-06
# if no, stop here
# if yes, run Dunn's test
d=dunnTest(distance~data,
           data=bray_stat,
           method="bh")
# Comparison          Z      P.unadj        P.adj
# 1              16S - lcms -4.2893917 1.791631e-05 8.958156e-05
# 2          16S - metaT_fx  1.1665710 2.433837e-01 2.704263e-01
# 3         lcms - metaT_fx  5.4559627 4.870822e-08 4.870822e-07
# 4        16S - metaT_gene -2.8177176 4.836633e-03 9.673266e-03
# 5       lcms - metaT_gene  1.4716742 1.411089e-01 1.763861e-01
# 6   metaT_fx - metaT_gene -3.9842886 6.768259e-05 2.256086e-04
# 7         16S - metaT_MAG -0.9691513 3.324697e-01 3.324697e-01
# 8        lcms - metaT_MAG  3.3202405 8.993994e-04 2.248499e-03
# 9    metaT_fx - metaT_MAG -2.1357222 3.270206e-02 5.450343e-02
# 10 metaT_gene - metaT_MAG  1.8485663 6.452046e-02 9.217209e-02
p35=c(d[["res"]][["P.unadj"]])

p=c(p7,p14,p21,p35)
p.adjust(p)


#                        day 7        day14         day21         day 35
#16S - lcms             #3.955169e-04 #8.984091e-03 #4.385053e-05 #6.270709e-04
#16S - metaT_fx         1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00
#lcms - metaT_fx        #2.913540e-04 #3.299879e-02 #3.751441e-05 #1.948329e-06
#16S - metaT_gene       #3.151502e-03 #9.480594e-03 #1.395929e-02 1.112426e-01
#lcms - metaT_gene      1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00
#metaT_fx - metaT_gene  #2.411123e-03 #2.890800e-02 #1.263675e-02 #2.301208e-03
#16S - metaT_MAG        1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00
#lcms - metaT_MAG       1.328210e-01  2.146586e-01  5.813475e-02  #2.428378e-02
#metaT_fx - metaT_MAG   1.000000e+00  1.000000e+00  1.000000e+00  5.886371e-01
#metaT_gene - metaT_MAG 4.727236e-01  1.651784e-01  1.000000e+00  1.000000e+00



#######################
# figure 2B
#######################
mag_categories=read_excel("Supplementary_Data_4.xlsx",sheet="Classifications")

mag_total=mag_categories%>%
  ungroup()%>%
  group_by(time)%>%
  summarise(tot=n())

mag_categories%>%
  ungroup()%>%
  group_by(time,classification)%>%
  summarise(n=n())%>%
  left_join(.,mag_total,by="time")%>%
  mutate(abund=n/tot)%>%
  ggplot()+
  geom_area(aes(x=time,y=abund,fill=ordered(classification,levels=c("cat_only","gain_fx","responsive","lost_fx","un_only","resistant"))))+
  scale_fill_manual(values=c("#fe9929","#fdbf6f","#fb8072","#a8ddb5","#7bccc4","#7570b3"))+
  theme_classic()
