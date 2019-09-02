library(tibble)
library(ggplot2)
library(class)
library(cluster)
library(factoextra)
library("NbClust")

data <- read.csv("USArrests.csv")

df <- scale(USArrests)     # Standaryzacja zbioru

#Wyznaczamy optymalną ilość klastrów; metoda ward.D2
res.nbclust <- NbClust(data = df, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 15, 
                       
                       method = "ward.D2", index = "all", alphaBeale = 0.1)


fviz_nbclust(res.nbclust, ggtheme = theme_minimal())

#Z analizy wynika, że optymalną ilością klastrów będzie  k = 2


######
#Metoda grupowania fuzzy clustering zwana też c-średnich
######

res.fanny <- fanny(df, 2)  # fuzzy clustering dla k = 2

head(res.fanny$clustering) # Grupy obserwacji


fviz_cluster(res.fanny, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(),
             legend = "right")

########
#Metoda grupowania K-średnich
#######
set.seed(123)
kus_arrests<-kmeans(df, centers =2, nstart = 50)


plot(x=df[,1], y=df[,2], col=kus_arrests$cluster)
points(kus_arrests$centers, pch=3, cex=2)
clusplot(df, kus_arrests$cluster, color = T, labels = 2, main = 'Cluster Plot')
