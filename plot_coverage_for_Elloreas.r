#!/usr/bin/Rscript
args <- commandArgs(TRUE)
infile_name <-args[1]
outfile_name <-args[2]
text_for_the_title <-args[3]

pdf(outfile_name,width=10,height=10) #in "pdf" width and height are in inches

a<-read.table(infile_name)
x<-a[,2]
y<-a[,3]

plot(x,y,xlim=c(0,max(x)),ylim=c(0,max(y)), type='l',main = text_for_the_title,xlab="position (bp)",ylab="coverage", xaxt='n', xaxs = "i", yaxs = "i", col="red")
#adding ticks to the X axis, one tick per 1000 bp.
x_marks_large=seq(0, max(x), by=1000)
axis(1, at=x_marks_large, labels=x_marks_large, col.axis="black", las=2, cex.axis=1.0, tck=-.02)