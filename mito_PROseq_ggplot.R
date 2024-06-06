# Plots PROseq data on the mitochondrial genome.
# Load libraries and select genome
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(reshape2)
genome = BSgenome.Hsapiens.UCSC.hg19

source("bigWig.R")

# Define genome region 
loc = GRanges("chrM:10451-10951:-")
loc = GRanges(loc, seqinfo=seqinfo(genome))

# Add PROseq data
CTRLC = as.data.frame(bw.coverage.stranded.chrM("CtrlC_FB", loc, flip.minus=T)[[1]])
CTRLD = as.data.frame(bw.coverage.stranded.chrM("CTRLD_FB", loc, flip.minus=T)[[1]])
CTRLG = as.data.frame(bw.coverage.stranded.chrM("CtrlG_FB", loc, flip.minus=T)[[1]])
CTRLH = as.data.frame(bw.coverage.stranded.chrM("CTRLH_FB", loc, flip.minus=T)[[1]])
CTRLL = as.data.frame(bw.coverage.stranded.chrM("CtrlL_FB", loc, flip.minus=T)[[1]])

# Combine datasets into a single dataframe and create column for coordinates
df6 = cbind(CTRLC, CTRLD, CTRLG, CTRLH, CTRLL)
names(df6) = c("CTRLC", "CTRLD", "CTRLG", "CTRLH", "CTRLL")
df6$coordinates = seq(nrow(df6))
final_df <- melt(df6, measure.vars = c("CTRLC", "CTRLD", "CTRLG", "CTRLH", "CTRLL"))

# Plot PROseq reads 
ggplot(final_df, aes(x=coordinates, y=value)) + 
      geom_col(color="red", size=1) +
      facet_grid(variable ~., scales = "free") +
      theme(panel.background=element_blank(), panel.grid.major=element_blank(), 
      panel.grid.minor=element_blank(), axis.line=element_line(size=0.8, 
      color="black", linetype=1)) + xlab("chrM:10451-10951 (- strand)") + ylab("RPM")
