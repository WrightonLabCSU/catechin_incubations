library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(Maaslin2)


paths=read_excel("Supplementary_Data_3.xlsx",sheet="carbon_cycle_pathways")

# read in metadata
metadata=read_excel("Supplementary_Data_3.xlsx",sheet="Metatranscriptomes")
metadata=metadata%>%
  select(Sample,Treatment,Time)


##################
## MaAslin2 = day 7
##################
# features as columns and samples as rows.
wide_fx_sums=paths%>%
  filter(time==7)%>%
  select(-time,-treat)%>%
  spread(key=path,value=sum,fill=0)
wide_fx_sums=as.data.frame(wide_fx_sums)
row.names(wide_fx_sums)=wide_fx_sums[,1]
wide_fx_sums=wide_fx_sums[,-1]

#Metadata file
#Formatted with features as columns and samples as rows.
meta=metadata%>%
  filter(Time==7)
meta=as.data.frame(meta)
row.names(meta)=meta$sample

Maaslin2(input_data = wide_fx_sums,
         input_metadata = meta,
         output="paths_day07_maaslin2",
         fixed_effects = c("treat"),
         plot_heatmap = FALSE,
         plot_scatter = FALSE)
##################
## MaAslin2 = day 14
##################
# features as columns and samples as rows.
wide_fx_sums=paths%>%
  filter(time==14)%>%
  select(-time,-treat)%>%
  spread(key=path,value=sum,fill=0)
wide_fx_sums=as.data.frame(wide_fx_sums)
row.names(wide_fx_sums)=wide_fx_sums[,1]
wide_fx_sums=wide_fx_sums[,-1]

#Metadata file
#Formatted with features as columns and samples as rows.
meta=metadata%>%
  filter(Time==14)
meta=as.data.frame(meta)
row.names(meta)=meta$sample

Maaslin2(input_data = wide_fx_sums,
         input_metadata = meta,
         output="paths_day14_maaslin2",
         fixed_effects = c("treat"),
         plot_heatmap = FALSE,
         plot_scatter = FALSE)

##################
## MaAslin2 = day 21
##################
# features as columns and samples as rows.
wide_fx_sums=paths%>%
  filter(time==21)%>%
  select(-time,-treat)%>%
  spread(key=path,value=sum,fill=0)
wide_fx_sums=as.data.frame(wide_fx_sums)
row.names(wide_fx_sums)=wide_fx_sums[,1]
wide_fx_sums=wide_fx_sums[,-1]

#Metadata file
#Formatted with features as columns and samples as rows.
meta=metadata%>%
  filter(Time==21)
meta=as.data.frame(meta)
row.names(meta)=meta$sample

Maaslin2(input_data = wide_fx_sums,
         input_metadata = meta,
         output="paths_day21_maaslin2",
         fixed_effects = c("treat"),
         plot_heatmap = FALSE,
         plot_scatter = FALSE)

##################
## MaAslin2 = day 35
##################
# features as columns and samples as rows.
wide_fx_sums=paths%>%
  filter(time==35)%>%
  select(-time,-treat)%>%
  spread(key=path,value=sum,fill=0)
wide_fx_sums=as.data.frame(wide_fx_sums)
row.names(wide_fx_sums)=wide_fx_sums[,1]
wide_fx_sums=wide_fx_sums[,-1]

#Metadata file
#Formatted with features as columns and samples as rows.
meta=metadata%>%
  filter(Time==35)
meta=as.data.frame(meta)
row.names(meta)=meta$sample

Maaslin2(input_data = wide_fx_sums,
         input_metadata = meta,
         output="paths_day35_maaslin2",
         fixed_effects = c("treat"),
         plot_heatmap = FALSE,
         plot_scatter = FALSE)
