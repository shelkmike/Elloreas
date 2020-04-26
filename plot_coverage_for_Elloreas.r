#!/usr/bin/Rscript
args <- commandArgs(TRUE)
infile_name <-args[1]
outfile_name <-args[2]

jpeg(outfile_name,width=1000,height=480, quality = 90)

a<-read.table(infile_name)
x<-a[,2]
y<-a[,3]
plot(x,y,xlim=c(0,max(x)),ylim=c(0,max(y)), type='l',main = "Coverage plot",xlab="position (bp)",ylab="coverage", xaxt='n', xaxs = "i", yaxs = "i", col="red")
#добавляем риски по оси икс через 1000 нуклеотидов
x_marks_large=seq(0, max(x), by=1000)
axis(1, at=x_marks_large, labels=x_marks_large, col.axis="black", las=0, cex.axis=1.0, tck=-.02)