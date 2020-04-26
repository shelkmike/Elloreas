#!/usr/bin/Rscript
args <- commandArgs(TRUE)
infile_name <-args[1]
outfile_name <-args[2]

jpeg(outfile_name,width=1000,height=1000, quality = 90)

dots = read.table(infile_name,header=T)

par(mar = c(5, 5, 5, 5)) # Setting the margin on all sides
plot(dots,type="l",xlab="position (bp)",ylab="position (bp)",main="Dot plot for the final contig produced by Elloreas",cex.main=2, cex.lab=1.5)