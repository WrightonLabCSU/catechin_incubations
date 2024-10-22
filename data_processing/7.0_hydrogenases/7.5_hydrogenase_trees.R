library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(readxl)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(cowplot)
library(phytools)

#####
# reading in FeFe tree
#####

# read in tree
tree = read.tree("7.3_combined_ref_pulled_fefe.mafftauto.gappyout.faa.treefile")
tree
# midpoint root?
tree=midpoint.root(tree)

# read in hydrogenase info
annotation = read.delim("7.6_fefe_hydrogenase_placements.txt",sep="\t",header = TRUE)

# undecorated tree
p=ggtree(tree, layout="circular", branch.length = "none")


# annotate tree
p+
  geom_fruit(data=annotation,geom=geom_tile,
             mapping=aes(y=label,x=0,fill=HYDDB),show.legend=TRUE,width=5,offset=0.1) +
scale_fill_manual(values=c("#fb6a4a","#fc9272","#fcbba1","#fee0d2","#fff5f0","#cb181d","#6a51a3","#54278f","#3f007d"))+
  new_scale_fill() + 
  geom_fruit(data=annotation,geom=geom_tile,
             mapping=aes(y=label,x=0,fill=reference),show.legend=TRUE,width=5,offset=0.1)+
  scale_fill_manual(values=c("white","black"))



#####
# reading in NiFe tree
#####

# read in tree
tree = read.tree("7.4_combined_ref_pulled_nife.mafftauto.gappyout.faa.treefile")
tree
# midpoint root?
tree=midpoint.root(tree)

# read in hydrogenase info
annotation = read.delim("7.7_nife_hydrogenase_placements.txt",sep="\t",header = TRUE)

# undecorated tree
p=ggtree(tree, layout="circular", branch.length = "none")


# annotate tree
p+
  geom_fruit(data=annotation,geom=geom_tile,
             mapping=aes(y=label,x=0,fill=HYDDB),show.legend=TRUE,width=5,offset=0.1) +
scale_fill_manual(values=c("#00441b","#006d2c","#238b45","#41ae76","#66c2a4","#99d8c9","#ccece6","#e5f5f9","#f7fcfd","#67001f","#980043","#ce1256","#e7298a","#df65b0","#662506","#993404","#cc4c02","#ec7014","#08306b","#08519c","#2171b5","#4292c6","#6baed6","#9ecae1","#c6dbef","#deebf7","#f7fbff","black"))+
  new_scale_fill() + 
  geom_fruit(data=annotation,geom=geom_tile,
             mapping=aes(y=label,x=0,fill=reference),show.legend=TRUE,width=5,offset=0.1)+
  scale_fill_manual(values=c("white","black"))

