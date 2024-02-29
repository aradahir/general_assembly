#!/usr/bin/env Rscript

library(argparse)
library(ggplot2)

parser <- ArgumentParser(description= 'add depth plots')

parser$add_argument('--input', '-i', help= 'input depth table file name')
parser$add_argument('--output', '-o', help= 'output depth plot file name')
parser$add_argument('--name', '-n', help= 'sample name')
xargs<- parser$parse_args()

dat <- read.table(xargs$input)
consensus_name = xargs$output
name = xargs$name

pdf(file = consensus_name, width = 4, height = 4)
plot(x= dat$V2, y = dat$V3, xlab = 'position(bp)', ylab = 'depth coverage', type = "S", 
     main = paste(name),ylim = c(0, max(dat$V3)+100))
dev.off()
