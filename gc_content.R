## Get GC content ready for every chromosome

library(getopt)
library(seqinr)

spec = matrix(c('ifile','i',1,"character",'ofile','o',1,"character"),byrow=TRUE,ncol=4)
opt=getopt(spec)
inputfile <- opt$ifile
outputfile <- opt$ofilehuam


chr_fa <- read.fasta(inputfile,forceDNAtolower = T)#fasta file

slidingwindowgc <- function(windowsize, inputseq)
{
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  gc_dat <- as.data.frame(matrix(ncol = 2,nrow = length(starts)))
  colnames(gc_dat) <- c("pos","gc_content")
  gc_dat[,1] <- as.character(starts)
  n <- length(starts)
  for (i in 1:as.numeric(n)) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)] 
    chunkgc <- GC(chunk)
    gc_dat[i,2] <- chunkgc
  }
  return(gc_dat)
}
binSize<-20000

chr_windowAll<-slidingwindowgc(binSize,chr8_fa[[1]])
write.table(chr_windowAll, file = outputfile,quote = F, col.names = F,row.names = F, sep='\t')

#the output would be pos and gc_content
