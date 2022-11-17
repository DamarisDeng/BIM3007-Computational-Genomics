# plot PI_hat
library(ggplot2)

table <- read.table('./step4-IBD.genome',header=T)
ggplot(data=table,aes(PI_HAT))+
  geom_histogram(binwidth=0.005)+
  xlim(c(0,0.1))+
  labs(title='Summary of PI hat')

