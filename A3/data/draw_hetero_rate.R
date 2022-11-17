library(ggplot2)
data <- read.table('step3-pruned.het', header=T)
data$n <- 1:nrow(data)
ggplot(data=data, aes(x=n,y=F))+
  geom_hline(aes(yintercept=0.05), color='red',linetype='dashed')+
  geom_hline(aes(yintercept=-0.05), color='red',linetype='dashed')+
  geom_point(size=0.3, alpha=0.5)+
  ggtitle("Summary of the heterozygosity rate")+
  xlab('')+
  ylab('Heterozygosity Rate')+
  ylim(c(-0.2, 0.2))+
  xlim(c(0, 1000))
library(dplyr)
filtered_data <- data %>% filter(F < -0.05 | F > 0.05)

library(ggplot2)
data <- read.table('step3.het', header=T)
data$n <- 1:nrow(data)
ggplot(data=data, aes(x=n,y=F))+
  geom_hline(aes(yintercept=0.05), color='red',linetype='dashed')+
  geom_hline(aes(yintercept=-0.05), color='red',linetype='dashed')+
  geom_point(size=0.3, alpha=0.5)+
  ggtitle("Summary of the heterozygosity rate")+
  xlab('')+
  ylab('Heterozygosity Rate')+
  ylim(c(-0.2, 0.2))+
  xlim(c(0, 1000))
# theme(axis.title.y = element_blank())
