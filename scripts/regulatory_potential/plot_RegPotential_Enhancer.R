## Plotting the "Regulatory Potential" of TE comapred to Enhancer Activity
## Abin Abrham
## Jan 31 2018


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

# lining up linages 
# boeroeutheira goes in between eutehria and euarchontoglires: https://en.wikipedia.org/wiki/Boreoeutheria
# Hominidae is between  Hominoidea and homininae
# Rodentia goes after Euarchontoglires and before eurarchonta https://en.wikipedia.org/wiki/Euarchonta
# Muridae goes right after Rodentia

# ---- LOAD DATA ------ 
setwd("~/Google Drive/transfer/scripts")
REG_POTEN_FILE="~/Google Drive/transfer/data/intersect_regPot_enhancer.tsv"
regPotential_wLin = read_delim(REG_POTEN_FILE, "\t", col_names = FALSE)

# ---- TIDY UP DATA ------ 
df = as.data.table(regPotential_wLin)
names(df) = c("TE_chr","TE_start", "TE_end", "RegPotential_score", "TE", "TEFamily", "Lineage", "TF_count", 
              "enhancer_chr","enhancer_start","enhancer_end", "enhancer_seq",
              "overlap","blank")

df[,c("enhancer_chr","enhancer_start","enhancer_end","overlap","blank"):=NULL]
df$Lineage = factor(df$Lineage, levels = levels(factor(df$Lineage))[c(13, 11,17,5,1,15,12,4,3,14,6,16,2,9,7,8,10)])
df = transform(df, test=do.call(rbind, strsplit(TEFamily, '/', fixed=TRUE)), stringsAsFactors=F)
names(df)[10:11] = c("mainFamily", "subFamily")

df$enhancer_seq[df$enhancer_seq=="."] = FALSE
df$enhancer_seq[df$enhancer_seq != FALSE] = TRUE
df$enhancer_seq = factor(df$enhancer_seq)


# ---- ANALYZE DATA ------ 
# Quick plots exploring best representation of distribution 
# ggplot(df, aes(x=RegPotential_score)) + geom_density() + ggtitle("Raw RegPotential Score Distribution")
# boxplot(df$RegPotential_score ~ df$enhancer_seq, outline=FALSE, main="Raw Reg Potential Scores w/o Outliers", ylab="Raw Regulatory Potential Score", xlab="Enhancer Overlap")
# boxplot(df$RegPotential_score ~ df$enhancer_seq, outline=TRUE, main="Raw Reg Potential Scores w/ Outliers", ylab="Raw Regulatory Potential Score", xlab="Enhancer Overlap")
# ggplot(df, aes(Lineage, RegPotential_score, fill=enhancer_seq)) + geom_violin() + scale_y_log10()
# ggplot(df, aes(x=Lineage, y=RegPotential_score)) + geom_boxplot() + scale_y_continuous(limits = c(0, 15))


# RAW DATA plot 
rd = ggplot(df, aes(x=Lineage, y=RegPotential_score, color=mainFamily)) + geom_jitter(alpha=0.1, size=0.35) + scale_y_continuous(limits = c(0, 15))
rd = rd + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", legend.key=element_rect(fill='black'))+ guides(colour = guide_legend(nrow = 1, override.aes = list(size=4, alpha=1)))
rd = rd + labs(y="Raw Regulatory Potential Score", title = "Raw Regulatory Potential of TEs", subtitle= "TE grouped by Date of Exapation", caption= "each point represent one specific TE")
rd
ggsave("Jan31_2018_RawRegPot_colored_by_TEFamily.png", limitsize=FALSE)

# violin plots of enhancer vs non enhancer reg potential distribution 
edf = df[enhancer_seq==TRUE]
erd = ggplot(edf, aes(x=Lineage, y=RegPotential_score, color=mainFamily)) + geom_jitter(alpha=0.1, size=0.35) + scale_y_continuous(limits = c(0, 15))
erd = erd + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", legend.key=element_rect(fill='black'))+ guides(colour = guide_legend(nrow = 1, override.aes = list(size=4, alpha=1))) 
erd = erd + labs(y="Raw Regulatory Potential Score", title = "Raw Regulatory Potential of TEs that overlap Enhancers", subtitle= "TE grouped by Date of Exapation", caption= "each point represent one specific TE")
erd
ggsave("Jan31_2018_RawRegPot_TEwEnh_colored_by_TEFamily.png", limitsize=FALSE)




# Mean RegPotential Score grouped by TE main family across each lineage 
mean_df = df[,.(avg=mean(RegPotential_score)),by=.(Lineage,mainFamily)]
mean_df$Lineage = factor(mean_df$Lineage, levels = levels(factor(mean_df$Lineage))[c(13, 11,17,5,1,15,12,4,3,14,6,16,2,9,7,8,10)])
mean_df$mainFamily = factor(mean_df$mainFamily, levels = levels(factor(mean_df$mainFamily))[c(1,2,4,5,3,8,9,7,6,10)])

# plot
p = ggplot(mean_df, aes(Lineage, avg))
p = p + geom_jitter(aes(color=mainFamily, shape=mainFamily), size=2, width=0.2)
p = p + scale_color_manual(values=c("#F3350D", "#F3350D", "#00C509",
                                    "#00C509", "#0D53F3", "#56B4E9",
                                    "#56B4E9", "#E8820F", "#076F63",
                                    "#999999")) 
p = p + scale_shape_manual(values=c(16,1,17,2,18,22,0,11,8,3))
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top")
p = p + labs(y="Mean Regulatory Potential",  title = "Mean of Regulatory Potential by TE Main Families",
             caption= "each point represent one main-family")
p

# Mean RegPotential Score grouped by TE sub-family across each lineage 
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


