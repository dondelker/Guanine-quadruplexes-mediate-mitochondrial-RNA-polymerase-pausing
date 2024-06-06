# Plots PROseq data on the mitochondrial genome.
# Load libraries and select genome
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(reshape2)
genome = BSgenome.Hsapiens.UCSC.hg19

source("bigWig.R")

# Define genome region 
loc = GRanges("chrM:450-17019:-")

# Add PROseq data
Ctrl = as.data.frame(bw.coverage.stranded.chrM("A549", loc, flip.minus=T)[[1]])
Salt = as.data.frame(bw.coverage.stranded.chrM("A549_Salt", loc, flip.minus=T)[[1]])
G4s = as.data.frame(bw.coverage.stranded.chrM("MT_G4s", loc, flip.minus=F)[[1]])

# Combine datasets into a single dataframe and create column for coordinates
df3 = cbind(Ctrl,Salt, G4s)
names(df3) = c("Ctrl", "Salt","G4s")
df3$G4s <- replace(df3$G4s, df3$G4s == 0, NA)
df3$G4s <- df3$G4s + 40
df3$coordinates = seq(nrow(df3))

# Plot PROseq reads with G4 values 
ggplot(df3, aes(x=coordinates)) + 
      geom_col(aes(y=Ctrl), color="blue", size=1) +
      geom_point(aes(y=G4s), color="red", size=0.4) +
      theme(panel.background=element_blank(), panel.grid.major=element_blank(), 
      panel.grid.minor=element_blank(), axis.line=element_line(size=0.8, 
      color="black", linetype=1)) + xlab("chrM:1-16569 (- strand)") + ylab("RPM")
