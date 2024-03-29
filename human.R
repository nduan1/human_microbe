#####to get relative abundance of each nuclectide location and use gam model to create a smooth lines for visualization
#human.R will run by submitting job human.sh
#input file 
# e.g. SRR13215371_chr21_ab.tsv
#chr21	1	0
#chr21	2	0
#chr21	3	0
#chr21	4	0
#....
library(dplyr)
library(ggplot2)
library(gsubfn)
library(mgcv)
library(gratia)
library(fBasics)
library(ggplot2)
library(tidymv)
library(tidyverse)
library(data.table)
library(getopt)

### get input parameters ##
spec = matrix(c('ifile','i',1,"character",'ofile','o',1,"character",'mfile','m',1,"character"),byrow=TRUE,ncol=4)
opt=getopt(spec)
inputfile <- opt$ifile
outputfile <- opt$ofile
image <- opt$mfile

#########preprocess function#########
slidingwindowplot <- function(windowsize, inputseq)
{
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(starts)
  chunkbps <- numeric(n)
  chunkstats<- numeric(n)
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
    chunkmean <- mean(chunk)
    chunkstdv<-sd(chunk)
    chunkbps[i] <- chunkmean
    chunkstats[i]<-chunkstdv
  }
  return (list(starts,chunkbps,chunkstats))
}
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
binSize<-300000
column.types <- c("character", "numeric", "numeric")

#########data input and process#########
#input dataset
dat<-as.data.frame(fread(inputfile, header=FALSE, sep="\t"))
myvector_all<-as.vector(as.matrix(dat[3]))
windowAll<-slidingwindowplot(binSize,myvector_all)
df<-data.frame(windowAll[[1]],windowAll[[2]],windowAll[[3]])
colname<-c("x","mean","sd")
colnames(df)<-colname
df[,-1] <- remove_outliers(df$mean)
df[, -1] <- lapply(df[, -1], function(x){(x/sum(x, na.rm=TRUE))*100} )
df_trans <- df[,-3]
df_trans$mean[df_trans$mean==0] <- NA
########model generate##########
gam_mod_pos <- gam(((mean))~s(x,bs="cr", k=100),
                   data = na.omit(df_trans), method = "REML", family = Gamma(link = "log"))
#########model fit check########
modelfit<-appraise(gam_mod_pos, method = "simulate")
##########output dataset##########
df_new <-  mutate(as.data.frame(na.omit(df_trans)),
                  resid = resid(gam_mod_pos),
                  predt = predict(gam_mod_pos))

write.table(df_new, file = outputfile,quote = F, col.names = F,row.names = F, sep='\t')

pdf(image)
ggplot(as.data.frame(df_new),aes(x = x, y = log(mean))) + geom_point() + scale_color_manual(values = c("blue"))+ xlab("position") + ylab("log(Coverage percentage)") + theme_bw() + geom_line(aes(x=x, y = (predt)),data = as.data.frame(df_new),col="red")
dev.off()
# the output df_new format is 
#e.g. SRR13215377_chr21_ab.output.tsv
#5100001 1.11268025489968        0.741023005722934       -0.552977021113403
#5400001 0.281181482163505       -0.652138280856498      -0.537340749809931
#5700001 0.34084241152694        -0.507810726366282      -0.521655342654666
#6000001 0.361442226883144       -0.471590380699301      -0.505895293300911
#6300001 0.452842596876521       -0.287669754277603      -0.490061559283455
#....

# already marked the significant position at each chrmosome use red lines based on the output from human_stats.R
