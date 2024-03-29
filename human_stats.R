###use the relative abundance to run the Mann-Whitney U test to find the significant location and generate figures
library(fs)
library(tidyverse)
library(colorspace) 
library(RColorBrewer)
library(directlabels)
library(ggplot2)
library(getopt)
library(stats)
### get input parameters ##
spec = matrix(c('ifile','i',1,"character",'ofile','o',1,"character",'mfile','m',1,"character"),byrow=TRUE,ncol=4)
opt=getopt(spec)
inputfile <- opt$ifile
outputfile <- opt$ofile
image <- opt$mfile

#1 read in files
meta_data <- read.csv("/data/EmiolaLab/duann2/human_proj/human_R_proj/dataset/metadata.csv",header = T,sep = ",")
#metadata.csv
#         Run    BioSample    Sex  Tumor
# SRR13215457 SAMN16954469   male Cancer
# SRR13215456 SAMN16954470   male Cancer
# SRR13215438 SAMN16954471   male Cancer
# SRR13215427 SAMN16954472   male Cancer
# SRR13215416 SAMN16954473 female Cancer
# SRR13215405 SAMN16954474   male Cancer
#...
file_paths <- dir_ls(paste("/data/EmiolaLab/duann2/human_proj/human_R_proj/results/",inputfile,sep = ""))
#here inputfile=a folder chr n contains 107 samples (SRR13215377_chr1_ab.output.tsv, SRR13215376_chr1_ab.output.tsv, SRR13215375_chr1_ab.output.tsv...)
# SRR13215377_chr1_ab.output.tsv
#1	0.117292428001571	0.165227430928148	-2.30388505863838
#900001	0.0969952465618428	-0.0540381459465307	-2.27856405841745
#1200001	0.094335312141923	-0.0898698290455549	-2.26966336153287
#1500001	0.104381918859161	0.000669332504155599	-2.26036806738177
#1800001	0.101916614679013	-0.0328400436697558	-2.25057952664753
#....

file_contents <- list()
for (i in seq_along(file_paths)) {
  file_contents[[i]] <- read.csv(
    file = file_paths[[i]],header=FALSE, sep="\t", 
  )
  colnames(file_contents[[i]]) <- c("position","mean","resid","predt")
  file_contents[[i]] <- mutate(file_contents[[i]],
                               sample=rep(substr(path_file(file_paths[i]),1,11),dim(file_contents[[i]])[1]))
}

file_contents <- set_names(file_contents, file_paths)#give each list a name

outdat <- data.frame(matrix("",nrow=0,ncol=3))
for (i in 1:length(file_paths)){
  outdat <- rbind(outdat, file_contents[[i]])
}
outdat_all <- left_join(outdat,meta_data,by=c("sample"="Run"))
####setup the color####
pal_cancer <- sequential_hcl(n=43,h1=0,h2=-100, c1=80,c2= 40,l2=75, alpha = 0.3)
pal_normal <- sequential_hcl(n=52,h1=240,h2=130, c1=30,c2=33, l2=95, alpha = 1)
cols <- c(pal_normal,pal_cancer)
outdat_all.sub <- subset(outdat_all,sample!="SRR13215428"&sample!="SRR13215429"&
                           sample!="SRR13215430"&sample!="SRR13215431"&
                           sample!="SRR13215432"&sample!="SRR13215436"&
                           sample!="SRR13215439"&sample!="SRR13215440"&
                           sample!="SRR13215441"&sample!="SRR13215442"&
                           sample!="SRR13215448"&sample!="SRR13215449")

pdf(image,width = 20, height = 5)
ggplot(as.data.frame(outdat_all.sub),aes(x = position, y = log(mean), color=(sample))) + scale_color_manual(values = cols)+ ggtitle(as.character(inputfile))+xlab("Genomeic Location (10Kbp)") + ylab("log(Coverage percentage)") + theme_bw() + theme(legend.position = "none")+geom_line(aes(x=position, y = (predt),color=(sample)),data = as.data.frame(outdat_all.sub),size=0.3)
dev.off()

###Mann-Whitney U test######
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

posi <- unique(outdat_all.sub$position)
n <- length(posi)
npvals <- numeric(n)
for (i in seq_along(posi)){
  sub.data.normal <- subset(outdat_all.sub,position==as.character(posi[i])&Tumor=="Normal")
  sub.data.cancer <- subset(outdat_all.sub,position==as.character(posi[i])&Tumor=="Cancer")
  sub.data.normal[,2] <- remove_outliers(sub.data.normal$mean)
  sub.data.cancer[,2] <- remove_outliers(sub.data.cancer$mean)
  man.data <- na.omit(rbind(sub.data.normal,sub.data.cancer))
  man.data$Tumor <- as.factor(man.data$Tumor)
  len<- length(levels(man.data$Tumor))
  if (len==2){
  	man <- wilcox.test(mean~Tumor,data=man.data)
  	npval <- man$p.value
  	npvals[i] <- npval
  }
  else {
  	npvals[i] <- NA
  }
}

p.adj <- p.adjust(npvals)
sig.data <- as.data.frame(cbind(posi,p.adj))
sig.data.sub <- subset(na.omit(sig.data),p.adj<0.05)
colnames(sig.data.sub) <- c("position","p-value")

write.table(sig.data.sub, file = outputfile,quote = F, col.names = F,row.names = F, sep='\t')
##output would be significant position and adjusted p value
#2100001	0.00184458897844473
#3600001	0.0242266878102262
#5100001	0.0334101219595682
#5400001	0.0419355963002726
#8400001	0.0113710026093534
#15900001	0.033519833491915
#19500001	0.0214854240077041
#...





