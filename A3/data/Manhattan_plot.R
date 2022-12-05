library("qqman")

args <- commandArgs(trailingOnly = TRUE)

# check if both input and output file names are provided 
if (length(args)<2) 
{
  stop("Please provided input and output file names", call.=FALSE)
}

#input and output file names
input <- args[1]
output <- args[2]

gwasResults <- read.table(input, head=TRUE)
jpeg(output, width=8, height=4, units="in", res=200)

manhattan(gwasResults, main = "Manhattan Plot (Before QC)" , col = c("blue4", "orange3"), annotateTop=TRUE, logp=TRUE, annotatePval=1E-20, ylim=c(0,25))
dev.off()