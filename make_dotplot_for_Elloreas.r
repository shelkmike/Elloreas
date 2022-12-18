#!/usr/bin/Rscript
args <- commandArgs(TRUE)
infile_name <-args[1]
outfile_name <-args[2]
text_for_the_title <-args[3]

pdf(outfile_name,width=10,height=10) #in "pdf" width and height are in inches

a<-read.table(infile_name,header=T)
x<-a[,1]
y<-a[,2]

#the function to calculate maximum when "N/A" values are present among the data
my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

#calculating maximum, while ignoring "N/A" values
maximum_value_of_x<-my.max(x)
maximum_value_of_y<-my.max(y)

par(mar = c(5, 5, 5, 5)) # Setting the margin on all sides
plot(x,y,xlim=c(0,maximum_value_of_x),ylim=c(0,maximum_value_of_y),type="l",main=text_for_the_title,xlab="position (bp)",ylab="position (bp)", xaxs = "i", yaxs = "i", cex.main=2, cex.lab=1.5, las=2)


