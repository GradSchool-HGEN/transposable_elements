### Analyze PCA for TFs on LTRs
### Feb 16 2018
### Abin Abraham 


library(data.table)
library(tidyr)
library(dplyr)
library(FactoMineR) 
library(factoextra)
library(corrplot) 

ROOT_PATH = "/Users/abin-personal/Google Drive/GD_transfer/data/"
SAVE_PATH = "/Users/abin-personal/Google Drive/GD_transfer/output/"
RAW_FILE = "intersect_HERV_TFBS_facet_enhancers.tsv"


TFdf = read.table(file.path(ROOT_PATH,RAW_FILE), "\t", header=FALSE)
TFdf = as.data.table(TFdf)
names(TFdf) = c("TE_chr", "TE_start", "TE_end", "TE", "TE_family", "TF_chr", "TF_start", "TF_end", "TF", "TE_info", "enhc_chr", "enhc_start", "enhc_end", "enhc_overlap")
TFdf = unite(TFdf, 'TEid', c("TE_chr","TE_start","TE_end"), sep = c(":","-"), remove = FALSE)
TFdf$enhc_overlap <- ifelse( TFdf$enhc_chr == ".", FALSE, TRUE)
head(TFdf) 

small = TFdf # a subset to make testing easier

df=small[, .(.N), by = .(TEid, TF)]
df = spread(df, TF, N )
df = as.data.frame(df)
df[is.na(df)] = 0 
head(df)

#Variable Map
df.pca = PCA(df[,2:dim(df)[2]], scale.unit = TRUE, ncp = 3, graph = TRUE)
fviz_pca_var(df.pca, labelsize =2)



var = get_pca_var(df.pca) # get info on variables into the PCA 
##   Name
## 1 "var$coord"
## 2 "var$cor"
## 3 "var$cos2"
## 4 "var$contrib" "contributions of the variables"


# Contributions of variables to PC1
fviz_contrib(df.pca, choice = "var", axes = 1, top = 50)
# Contributions of variables to PC2
fviz_contrib(df.pca, choice = "var", axes = 2, top = 50)



fviz_eig(df.pca, addlabels = TRUE, ylim = c(0, 120))
head(df.pca$coord, 4)

head(var$cos2, 4)
corrplot(var$contrib, is.corr=FALSE)

# Contributions of variables to PC1
fviz_contrib(df.pca, choice = "var", axes = 1, top = 50)
# Contributions of variables to PC2
fviz_contrib(df.pca, choice = "var", axes = 2, top = 50)

# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(123)
res.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(df.pca, col.var = grp,
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster", labelsize=3)

res.desc <- dimdesc(df.pca, axes = c(1,2), proba = 0.05)


# Description of dimension 1
res.desc$Dim.1
# Description of dimension 2
res.desc$Dim.2

