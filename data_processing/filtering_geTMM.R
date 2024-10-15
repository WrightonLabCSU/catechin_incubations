library(writexl)
library(readxl)
library(tidyr)
library(dplyr)
library(edgeR)
library(ggplot2)

setwd("/Users/mcgivern.9/Desktop/Projects/EMERGE/2016_incubations/data/metaT/JGI_100M_2304MAGs/")

##################
## starting with reverse stranded
##################
# read in htseq table
# this table has rows with all zeros removed
counts = read.delim("htseq_2304MAGs_100M_97_REVSTRANDED_no0s.txt",sep="\t",header = FALSE)
colnames(counts)=c("gene","STM_0716_E_M_E002","STM_0716_E_M_E003","STM_0716_E_M_E004","STM_0716_E_M_E025","STM_0716_E_M_E027","STM_0716_E_M_E029","STM_0716_E_M_E030","STM_0716_E_M_E031","STM_0716_E_M_E033","STM_0716_E_M_E034","STM_0716_E_M_E035","STM_0716_E_M_E050","STM_0716_E_M_E051","STM_0716_E_M_E052","STM_0716_E_M_E054","STM_0716_E_M_E055","STM_0716_E_M_E056","STM_0716_E_M_E058","STM_0716_E_M_E059","STM_0716_E_M_E060","STM_0716_E_M_E062","STM_0716_E_M_E063","STM_0716_E_M_E064","STM_0716_E_M_E066","STM_0716_E_M_E067","STM_0716_E_M_E068","STM_0716_E_M_E070","STM_0716_E_M_E071","STM_0716_E_M_E072","STM_0716_E_M_E121","STM_0716_E_M_E122","STM_0716_E_M_E123","STM_0716_E_M_E125","STM_0716_E_M_E126","STM_0716_E_M_E127","STM_0716_E_M_E129","STM_0716_E_M_E130","STM_0716_E_M_E131")
# get lib size
lib=as.data.frame(colSums(counts[,2:39]))
lib$no_feature=t(counts[2022482,2:39])
lib$ambig=t(counts[2022483,2:39])
colnames(lib)=c("sum","no_feature","ambig")

#plot mapping stats
lib$mapped=(lib$sum-lib$no_feature-lib$ambig)
lib$sample=row.names(lib)
lib%>%filter(sample!="sum")%>%gather(-sample,-sum,key="status",value="counts")%>%
  ggplot()+
  geom_bar(aes(x=sample,y=counts,fill=status),stat="identity",position="stack")+
  theme_classic()+
  geom_vline(xintercept=50000000)+
  coord_flip()

#remove bottom two stats rows
counts=counts[1:2022481,]
row.names(counts)=counts[,1]
counts=counts[,-1]

# read in gene lengths
len=read.delim("gene_lengths.txt",sep="\t",header = FALSE)
colnames(len)=c("gene","length")

# moving on to filtering
counts$sum=rowSums(counts)
range(counts$sum) # range 1 to 975112

filtered_counts=counts %>%
  filter(sum>0)
range(filtered_counts$sum) # range 1 to 975112, because i already removed rows of all zero

#  making read counts less than 5 = 0
counts5 = filtered_counts
counts5[counts5==1|counts5==2|counts5==3|counts5==4] <-0

pa_counts=ifelse(counts5>0,1,0)
pa_counts=as.data.frame(pa_counts)
pa_counts=pa_counts[,-39] # get rid of old sum column
pa_counts$n=rowSums(pa_counts)

# just has to be in 1 sample
n2 = pa_counts %>%
  filter(n>=1)

x = row.names(n2)
x=as.data.frame(x)
x$n=n2$n
filtered_counts$name=row.names(filtered_counts)

filtered_1x = inner_join(filtered_counts,x,by=c("name"="x"))
colnames(filtered_1x)
filtered_1x = filtered_1x[,c(40,1:38)]
colnames(filtered_1x)

#write_xlsx(filtered_1x,"filtered_5counts_1sample_Rev.xlsx") # too big to write
write.csv(filtered_1x,"filtered_5counts_1sample_Rev.txt",sep = "\t")

# convert to geTMM
count=filtered_1x
count=left_join(count,len,by=c("name"="gene"))
count5=count%>%distinct()
rpk = ((count5[,2:39]*10^3)/count5$length)
row.names(rpk)=count5$name
group <- c(rep("A",ncol(rpk)))

rpk.norm <- DGEList(counts=rpk,group=group)
lib
lib_size=c(100000000,98819300,95400980,78691532,100000000,88540314,98819300,100000000,83168470,100000000,99235256,91608098,93570944,80222604,81306340,100000000,100000000,100000000,84557472,100000000,88283948,100000000,100000000,100000000,100000000,96316002,63972406,61512094,100000000,100000000,100000000,100000000,100000000,100000000,100000000,78973250,100000000,100000000)
lib_size=as.data.frame(lib_size)
lib_size$sample=c(colnames(filtered_1x[,2:39]))
row.names(lib_size)=lib_size$sample
lib_size=select(lib_size,-sample)

rpk.norm$samples$lib.size = lib_size[rownames(rpk.norm$samples),]
rpk.norm$samples$lib.size
rpk.norm$samples
rpk.norm <- calcNormFactors(rpk.norm)

input=rpk.norm$samples
input$sample=row.names(input)
input_getmm=cpm(rpk.norm)
input_getmm=as.data.frame(input_getmm)
input_getmm$gene=row.names(input_getmm)
input_getmm=input_getmm%>%
  gather(-gene,key="sample",value="geTMM")%>%
  filter(gene=="STM_0716_E_M_E069_A_bin.67_k121_1854211_13"|
           gene=="STM_0716_E_M_E058_A_bin.18_k121_1959993_17"|
           gene=="STM_0716_E_M_E026_E030_E034_A_bin.12_k121_7316787_117"|
           gene=="STM_0716_E_M_E058_A_bin.18_k121_584664_4"|
           gene=="20120600_E2D_38_c_000000010023_10"|
           gene=="STM_0716_E_M_E026_E030_E034_A_bin.12_k121_7316787_28")
input_getmm=input_getmm%>%mutate(type="input")

#
lib_size=c(3061262,3561254,2668554,4975801,4481615,2841390,3561254,3720585,6442429,6322189,8711581,5097930,4098974,3402308,2404051,3621066,3676644,6650077,9428196,8077187,4745179,3782729,5967342,4915275,5790416,5052582,4047315,3959982,5073720,6029949,4287671,5705507,2919862,2739061,3012380,5459946,4963310,4807891)
lib_size=as.data.frame(lib_size)
lib_size$sample=c(colnames(filtered_1x[,2:39]))
row.names(lib_size)=lib_size$sample
lib_size=select(lib_size,-sample)
lib$mapped
rpk.norm$samples$lib.size = lib_size[rownames(rpk.norm$samples),]
rpk.norm$samples$lib.size
rpk.norm$samples

rpk.norm <- calcNormFactors(rpk.norm)

mapped=rpk.norm$samples
mapped$sample=row.names(mapped)
comp=left_join(input,mapped,by=c("group","sample"))
comp=comp%>%
  mutate(ratio=norm.factors.y/norm.factors.x)
mapped_getmm=cpm(rpk.norm)
mapped_getmm=as.data.frame(mapped_getmm)
mapped_getmm$gene=row.names(mapped_getmm)
mapped_getmm=mapped_getmm%>%
  gather(-gene,key="sample",value="geTMM")%>%
  filter(gene=="STM_0716_E_M_E069_A_bin.67_k121_1854211_13"|
           gene=="STM_0716_E_M_E058_A_bin.18_k121_1959993_17"|
           gene=="STM_0716_E_M_E026_E030_E034_A_bin.12_k121_7316787_117"|
           gene=="STM_0716_E_M_E058_A_bin.18_k121_584664_4"|
           gene=="20120600_E2D_38_c_000000010023_10"|
           gene=="STM_0716_E_M_E026_E030_E034_A_bin.12_k121_7316787_28")
mapped_getmm=mapped_getmm%>%mutate(type="mapped")

full_join(mapped_getmm,input_getmm)%>%
  full_join(.,z,by=c("gene"="name","sample","geTMM","type"))%>%
  ggplot()+
  geom_point(aes(x=sample,y=geTMM,color=type),alpha=0.5)+
  geom_line(aes(x=sample,y=geTMM,group=type),alpha=0.5)+
  theme_classic()+
  facet_grid(gene~type,scales="free_x")+
  coord_flip()

z=count5%>%filter(name=="STM_0716_E_M_E069_A_bin.67_k121_1854211_13"|name=="STM_0716_E_M_E058_A_bin.18_k121_1959993_17"|name=="STM_0716_E_M_E026_E030_E034_A_bin.12_k121_7316787_117"|name=="STM_0716_E_M_E058_A_bin.18_k121_584664_4"|name=="20120600_E2D_38_c_000000010023_10"|name=="STM_0716_E_M_E026_E030_E034_A_bin.12_k121_7316787_28")%>%
  gather(-name,-length,key="sample",value="geTMM")%>%mutate(type="count")%>%select(-length)
  
#

getmms <- cpm(rpk.norm)

getmms_df = as.data.frame(getmms)
range(getmms_df) #0.000 8082.522


getmms_df$gene=row.names(getmms)
colnames(getmms_df)
getmms_df = getmms_df[,c(39,1:38)]
colnames(getmms_df)
#write_xlsx(getmms_df, 'getmms_ge5_1X.xlsx') # too big
write.table(getmms_df, 'getmms_ge5_1X_REV_MAPPEDLIB.txt', sep="\t")



comp$input=c(100000000,98819300,95400980,78691532,100000000,88540314,98819300,100000000,83168470,100000000,99235256,91608098,93570944,80222604,81306340,100000000,100000000,100000000,84557472,100000000,88283948,100000000,100000000,100000000,100000000,96316002,63972406,61512094,100000000,100000000,100000000,100000000,100000000,100000000,100000000,78973250,100000000,100000000)

tidy_getmm_ge22_cat.un%>%
  filter(gene=="20100900_E1D_10_c_000000000189_50")%>%
  ggplot()+
  geom_point(aes(x=sample,y=geTMM))+
  theme_classic()
