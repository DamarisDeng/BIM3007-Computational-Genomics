library(ggplot2)
data <- read.table('./step8-dif-miss.missing', header=T)
data$logP <- log10(data$P)
ggplot(data=data,aes(x=logP))+
  geom_freqpoly()+
  xlim(c(-5,0))+
  labs(x='log10P differential missingness', 
       y='Count', 
       title='Summary of differential missingness')
