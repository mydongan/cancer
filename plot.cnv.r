#rm(list=ls())
#setwd('D:/desk/R_tmp')

#dd <- read.delim('rawdata.txt',header=T,sep="\t")
#info <- read.delim("plussize",header=F,sep="\t")
allargs<- commandArgs(trailingOnly = T)
input<-allargs[1]
output <-allargs[2]

dd<-read.delim(input,header=T,sep='\t')
info<-read.delim("/data/database/reference/Homo_sapiens/ncbi/GRCh37.p13/other/database/Sig_gene/plussize",header=F,sep="\t")


library(plotrix);
library(colorRamps);

dat<-c();
for(i in 1:nrow(dd) ) {
  for(j in 5:length(dd[i,]) ) {
    samp <- colnames(dd)[j]
    chr <- dd[i,1]
    if(chr )
      str <- dd[i,2]+info[which(info[,1]==chr),3]
    en <- dd[i,3]+info[which(info[,1]==chr),3]
    value <- dd[i,j]
    data <- data.frame(sample=samp,start=str,end=en,CN=value)
    dat <- rbind(dat,data)
  }
}


my_palette <- t(col2rgb(colorRampPalette(c("blue", "white", "red"))(n = 5)))
colors<-color.scale(as.numeric(dat$CN),my_palette[,1],my_palette[,2],my_palette[,3]);
color.legend<-color.gradient(my_palette[,1],my_palette[,2],my_palette[,3],nslices=30);
lev <- seq(-1,1,length.out = 30)
getcolor <- function(CNV)
{
  print(CNV)
  if(CNV==-1){return(color.legend[1])}else if(CNV == 1){return(color.legend[30])}
  for(i in 1:length(lev)-1)
  {
    if( CNV>=lev[i] && CNV<=lev[i+1] )
    {
      return(color.legend[1+i])
    }
  }
}
getcolor(0.03)
colors <- c()
for(cn in dat$CN){
  colors=c(colors,getcolor(cn))
}

dat$color <- colors;


len <- (ncol(dd)-4)

sample=colnames(dd)[5:ncol(dd)]
max = info[nrow(info),2] + info[nrow(info),3]
per=(max/60)
chr=info[,1]


  pos=c(as.numeric((info[,2]/2+info[,3]))/per)
  pos[c(16:22)] <- pos[c(16:22)]-0.1
  pos[22] <- 60
  pos[21] <- 58.7
  pos[20] <- 57.3
png(paste0(output,'.png'),width=1600,height=(len*20 + 120))
par(mar=c(4,20,4,4))
plot(1:len,1:len,type="n",ylim=c(0,len+1),xlim=c(0,60),bty="n",yaxt="n",xaxt="n",xlab="Chromosome",ylab="",xaxs="i");
axis(2,at=(c(1:len)-0.5),labels=sample,las=2,col = "NA", col.ticks = "NA");
axis(1,at=(pos),labels=chr,las=1,col = "NA", col.ticks = "NA");
axis(1,at=(unique(as.numeric(c(max,info[,3])))/per),labels=NA,las=1,col = "black", col.ticks = "black");
for (i in 1:nrow(dat)) {
  subset<-dat[i,];
  x.left<-subset$start/per;
  x.right<-subset$end/per;
  y.top<- which(subset$sample == sample);
  y.bottom <- y.top-1
  col<-as.character(subset$color);
  rect(x.left,y.bottom,x.right,y.top,col=col,border="NA");
}
lines(c(0,60), c(0,0), type = 'l')
lines(c(0,60), c(len,len), type = 'l')
lines(c(0,0), c(0,len), type = 'l')
lines(c(60,60), c(0,len), type = 'l')
color.legend(49,len+1,59,len+1.5,c(-1,1),color.legend,gradient="x");
dev.off()
pos=c(as.numeric((info[,2]/2+info[,3]))/per)
pos[c(16:22)] <- pos[c(16:22)]-0.1
pos[22] <- 60
pos[21] <- 58.7
pos[20] <- 57.3
pdf(paste0(output,'.pdf'),width=17.8,height=17.8*0.6)
par(mar=c(4,20,4,4))
plot(1:len,1:len,type="n",ylim=c(0,len+1),xlim=c(0,60),bty="n",yaxt="n",xaxt="n",xlab="Chromosome",ylab="",xaxs="i");
axis(2,at=(c(1:len)-0.5),labels=sample,las=2,col = "NA", col.ticks = "NA");
axis(1,at=(pos),labels=chr,las=1,col = "NA", col.ticks = "NA");
axis(1,at=(unique(as.numeric(c(max,info[,3])))/per),labels=NA,las=1,col = "black", col.ticks = "black");
for (i in 1:nrow(dat)) {
  subset<-dat[i,];
  x.left<-subset$start/per;
  x.right<-subset$end/per;
  y.top<- which(subset$sample == sample);
  y.bottom <- y.top-1
  col<-as.character(subset$color);
  rect(x.left,y.bottom,x.right,y.top,col=col,border="NA");
}
lines(c(0,60), c(0,0), type = 'l')
lines(c(0,60), c(len,len), type = 'l')
lines(c(0,0), c(0,len), type = 'l')
lines(c(60,60), c(0,len), type = 'l')
color.legend(49,len+1,59,len+1.5,c(-1,1),color.legend,gradient="x");
dev.off()
