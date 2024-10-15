library(writexl)
library(readxl)
library(tidyr)
library(dplyr)
library(edgeR)
library(ggplot2)



# read in htseq table
# this table has rows with all zeros removed
# find this table on zenodo: 10.5281/zenodo.13375971
counts = read.delim("htseq_2302MAGs_100M_97_REVSTRANDED_no0s.txt",sep="\t",header = FALSE)
colnames(counts)=c("gene","STM_0716_E_M_E002","STM_0716_E_M_E003","STM_0716_E_M_E004","STM_0716_E_M_E025","STM_0716_E_M_E027","STM_0716_E_M_E033","STM_0716_E_M_E034","STM_0716_E_M_E035","STM_0716_E_M_E050","STM_0716_E_M_E051","STM_0716_E_M_E052","STM_0716_E_M_E058","STM_0716_E_M_E059","STM_0716_E_M_E060","STM_0716_E_M_E062","STM_0716_E_M_E063","STM_0716_E_M_E064","STM_0716_E_M_E070","STM_0716_E_M_E071","STM_0716_E_M_E072","STM_0716_E_M_E121","STM_0716_E_M_E122","STM_0716_E_M_E123","STM_0716_E_M_E129","STM_0716_E_M_E130","STM_0716_E_M_E131")

# get lib size
lib=as.data.frame(colSums(counts[,2:27]))
lib$no_feature=t(counts[1760697,2:27])
lib$ambig=t(counts[1760698,2:27])
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
counts=counts[1:1760696,]
row.names(counts)=counts[,1]
counts=counts[,-1]

# read in gene lengths
len=read.delim("gene_lengths.txt",sep="\t",header = FALSE)
colnames(len)=c("gene","length")

# moving on to filtering
counts$sum=rowSums(counts)
range(counts$sum) # range 1 to 974827

filtered_counts=counts %>%
  filter(sum>0)
range(filtered_counts$sum) # range 1 to 974827, because i already removed rows of all zero

#  making read counts less than 5 = 0
counts5 = filtered_counts
counts5[counts5==1|counts5==2|counts5==3|counts5==4] <-0

pa_counts=ifelse(counts5>0,1,0)
pa_counts=as.data.frame(pa_counts)
pa_counts=pa_counts[,-27] # get rid of old sum column
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
filtered_1x = filtered_1x[,c(28,1:26)]
colnames(filtered_1x)

# convert to geTMM
count=filtered_1x
count=left_join(count,len,by=c("name"="gene"))
count5=count%>%distinct()
rpk = ((count5[,2:27]*10^3)/count5$length)
row.names(rpk)=count5$name
group <- c(rep("A",ncol(rpk)))

rpk.norm <- DGEList(counts=rpk,group=group)
lib
lib_size=c(100000000,98819300,95400980,
           78691532,100000000,
           83168470,100000000,99235256,
           91608098,93570944,80222604,
           100000000,84557472,100000000,
           88283948,100000000,100000000,
           63972406,61512094,100000000,
           100000000,100000000,100000000,
           78973250,100000000,100000000)
lib_size=as.data.frame(lib_size)
lib_size$sample=c(colnames(filtered_1x[,2:27]))
row.names(lib_size)=lib_size$sample
lib_size=select(lib_size,-sample)

rpk.norm$samples$lib.size = lib_size[rownames(rpk.norm$samples),]
rpk.norm$samples$lib.size
rpk.norm$samples
rpk.norm <- calcNormFactors(rpk.norm)
getmms <- cpm(rpk.norm)

getmms_df = as.data.frame(getmms)
range(getmms_df) #0.000 9683.861


getmms_df$gene=row.names(getmms)
colnames(getmms_df)
getmms_df = getmms_df[,c(27,1:26)]
colnames(getmms_df)
write.table(getmms_df, 'getmms_ge5_1X_REV_26samples.txt', sep="\t")

