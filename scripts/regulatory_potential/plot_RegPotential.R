## Plotting the "Regulatory Potential" of TE vs. TE Age 
## Abin Abrham
## Jan 30 2018


# install.packages("dplyr")
# install.packages("ggplot2")

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(reshape)
library(RColorBrewer)
rm(list = ls())

#lining up linages 
# boeroeutheira goes in between eutehria and euarchontoglires: https://en.wikipedia.org/wiki/Boreoeutheria
# Hominidae is between  Hominoidea and homininae
# Rodentia goes after Euarchontoglires and before eurarchonta https://en.wikipedia.org/wiki/Euarchonta
# Muridae goes right after Rodentia

# ---- LOAD DATA ------ 
setwd("~/Google Drive/transfer/scripts")
REG_POTEN_FILE="~/Google Drive/transfer/data/RegPotential_wLin.tsv"
regPotential_wLin = read_delim(REG_POTEN_FILE, "\t", escape_double = FALSE, trim_ws = TRUE)

# ---- TIDY UP DATA ------ 
df = as.data.table(regPotential_wLin)

#consider this: mean_sf$Lineage = factor(mean_sf$Lineage, levels = levels(factor(mean_sf$Lineage))[c(13, 11,17,5,1,15,12,4,3,14,6,16,2,9,7,8,10)])
levels(df$Lineage) = levels(df$Lineage)[c( 13, 11,17,5,1,15,12,4,3,14,6,16,2,9,7,8,10)]
df = transform(df, test=do.call(rbind, strsplit(Family, '/', fixed=TRUE)), stringsAsFactors=F)
names(df)[7:8] = c("mainFamily", "subFamily")
df[,TEcoord:=NULL]


# ---- ANALYZE DATA ------ 
ggplot(df, aes(x=RegPotential_score)) + geom_density() + ggtitle("Raw RegPotential Score Distribution")
boxplot(df$RegPotential_score, outline=FALSE, main="Raw Reg Potential Scores w/o Outliers")
boxplot(df$RegPotential_score, outline=TRUE, main="Raw Reg Potential Scores w/ Outliers")


boxplot(RegPotential_score ~ Lineage, data=df, las=2, cex.lab=0.4)
boxplot(RegPotential_score ~ Lineage, outline=FALSE, data=df, las=2, cex.axis=0.8, main ="Reg Potential w/o Outliers", ylab="regulatory potential")

# plot boxplots without outliers 
qplot(df$RegPotential_score, df$Lineage, geom="boxplot")
ggplot(df, aes(x=Lineage, y=RegPotential_score)) + geom_boxplot() + scale_y_continuous(limits = c(0, 15))
ggplot(df, aes(x=Lineage, y=RegPotential_score)) + geom_violin() + scale_y_continuous(limits = c(0, 15))

#Summarize by TE main family and Lineage 
mean_df = df[,.(avg=mean(RegPotential_score)),by=.(Lineage,mainFamily)]
mean_df$Lineage = factor(mean_df$Lineage, levels = levels(factor(mean_df$Lineage))[c(13, 11,17,5,1,15,12,4,3,14,6,16,2,9,7,8,10)])
mean_df$mainFamily = factor(mean_df$mainFamily, levels = levels(factor(mean_df$mainFamily))[c(1,2,4,5,3,8,9,7,6,10)])


# plot
p = ggplot(mean_df, aes(Lineage, avg))
p = p + geom_jitter(aes(color=mainFamily, shape=mainFamily), width=0.2)
p = p + scale_color_manual(values=c("#F3350D", "#F3350D", "#00C509",
                                 "#00C509", "#0D53F3", "#56B4E9",
                                 "#56B4E9", "#E8820F", "#076F63",
                                 "#999999")) 
p = p + scale_shape_manual(values=c(16,1,17,2,18,22,0,11,8,3))
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="top")
p = p + labs(y="Mean Regulatory Potential")
p

# now lets break it up by sub families
mean_sf = df[,.(avg=mean(RegPotential_score)), by=.(Lineage,subFamily,mainFamily)]
mean_sf$Lineage = factor(mean_sf$Lineage, levels = levels(factor(mean_sf$Lineage))[c(13, 11,17,5,1,15,12,4,3,14,6,16,2,9,7,8,10)])
mean_sf$mainFamily = factor(mean_sf$mainFamily, levels = levels(factor(mean_sf$mainFamily))[c(1,2,4,5,8,9,3,7,6,10)])

#plot
sf = ggplot(mean_sf, aes(Lineage, avg))
sf = sf + geom_jitter(aes(color=mainFamily), width=0.2, size =0.8)
sf = sf + theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="top")+ guides(colour = guide_legend(nrow = 1))
sf = sf + scale_color_manual(values=brewer.pal(10,"Paired")) 
sf = sf + labs(y="Mean Regulatory Potential", title = "Mean of Regulatory Potential by Sub Families", subtitle= "Colored by Main Family",
               caption= "each point represent one sub-family")
sf


