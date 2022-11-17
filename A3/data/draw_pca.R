library(ggplot2)
data <- read.table('./step5-pca.eigenvec',header=T)
data <- data[,3:ncol(data)]
names(data) <- paste0('PC',1:ncol(data)) # change the name of columns into PC1, PC2, ..., PCn
ggplot(data=data)+
  geom_point(aes(x=PC1, y=PC2))+
  labs(title="PCA plot", x ="PC1", y = "PC2")+
  xlim(c(-0.1,0.1))+
  ylim(c(-0.5,0.5))

values <- read.table('./step5-pca.eigenval')
# calculate the explained proportion
values$Percentage <- values$V1/sum(values$V1)
values$PC <- paste0(1:nrow(values))
ggplot(data=values, aes(x=PC, y=Percentage))+
  geom_col()+
  ylim(c(0,0.5))+
  labs(title="Explained variance of each eigenvalue")
