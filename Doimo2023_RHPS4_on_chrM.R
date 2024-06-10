# Plots RHPS4 data from Doimo 2023.

library(rtracklayer)

library(dplyr)
library(ggplot2)
library(reshape2)
library(tibble)
library(tidyr)

library(cowplot)

library(BSgenome.Hsapiens.UCSC.hg38)

library(bigWigUtils)


# where files are
data.dir = "/home/jtburd/data/seq/mitochondria/Doimo2023/G4_chrM_files/"
tsv.dir = paste0(data.dir, "chrM_readDepth/")
# possibly not using bigWig files
bw.dir = paste0(data.dir, "chrM_bigWigs/")

# region to plot
region = GRanges("chrM:1-16569:+")

genome = BSgenome.Hsapiens.UCSC.hg38

# read in all of the TSV files
rd = list()
for (f in list.files(tsv.dir)) {
    sample.name = gsub("_\\d+.tsv$", "", f)
    sample.name = gsub("INDUCED", "Induced", sample.name)
    x = read.table(paste0(tsv.dir, f))
    rd[[sample.name]] = x[,3]
}
rd = do.call(cbind, rd)

# read in all the bigWig files, and convert to a matrix, "rd.bw"
# (short for "read depth, from bigWig")
if (!exists("rd.bw")) {
    bw = list()
    for (f in list.files(bw.dir)) {
        sample.name = gsub("_chrM.bw$", "", f)
        bw[[sample.name]] = as.numeric(bw.read.depth(
            paste0(bw.dir, f), region)[[1]])
    }
    rd.bw = do.call(cbind, bw)
}

# rd, with each column standardized (to have mean=0, s.d.=1)
# (and converted to a data frame)
rd.1 = data.frame(scale(rd.bw))
# get various averages
rd.1$Induced_FLAG = (rd.1$Induced_FLAG_1 + rd.1$Induced_FLAG_2) / 2
rd.1$Induced_Input = (rd.1$Induced_Input_1 + rd.1$Induced_Input_2) / 2
rd.1$RHPS4_FLAG = (rd.1$RHPS4_FLAG_1 + rd.1$RHPS4_FLAG_2) / 2
rd.1$RHPS4_Input = (rd.1$RHPS4_Input_1 + rd.1$RHPS4_Input_2) / 2
# ...and ratios
rd.1$RHPS4_diff = rd.1$RHPS4_FLAG - rd.1$RHPS4_Input
rd.1$Induced_diff = rd.1$Induced_FLAG - rd.1$Induced_Input

# PROseq mean
proseq.mean = bw.read.depth.stranded("/home/jtburd/pub/hg19/stats/PROseq/Ctrl5_mean.plus.bw",
    GRanges("chrM:1-16569"))[[1]]

# x-coordinates, for plotting
rd.1$position = c(1:16569)

png("readDepth_vs_bigWig.png", width=6, height=6, units="in", res=75)
plot(rd.bw[,"RHPS4_FLAG_1"], rd[,"RHPS4_01_FLAG"], col=rgb(0,0,0,0.3))

dev.off()

# get averages of read depth
rd.average = list()
for (sample1.name in grep("_01_", colnames(rd), value=TRUE)) {
    sample2.name = gsub("_01_", "_02_", sample1.name)
    avg.name = gsub("_01_", "_", sample1.name)
    rd.average[[avg.name]] = (rd[,sample1.name] + rd[,sample2.name]) / 2
}
rd.average = do.call(cbind, rd.average)

# convert to long form
rd.long = rd.average %>% as_tibble %>%
    rowid_to_column("position") %>%
    gather(key="sample", value="read.depth", -1)
# remove one outlier
rd.long = rd.long[ rd.long$position != 3107 , ]

png("Doimo2023_averaged_signals.png",
    width=11, height=5, units="in", res=75)
g = (ggplot(rd.long, aes(x=position, y=read.depth))
    + geom_line(aes(color=sample)) + theme_classic())
plot(g)
dev.off()

# trying to normalize each of these to sum to 1
rd.proportion = t( t(rd.average) / apply(rd.average, 2, sum) )

# convert to long form
rd.long.2 = rd.proportion %>% as_tibble %>%
    rowid_to_column("position") %>%
    gather(key="sample", value="read.depth", -1)
# remove one outlier
rd.long.2 = rd.long.2[ rd.long.2$position != 3107 , ]

png("Doimo2023_summing_to_1.png",
    width=11, height=5, units="in", res=75)
g = (ggplot(rd.long.2, aes(x=position, y=read.depth)) 
    + geom_line(aes(color=sample)) + theme_classic())
plot(g)
dev.off()

# ratios of that normalized data
rd.ratio = data.frame(
    RHPS4 = rd.proportion[,"RHPS4_FLAG"] / rd.proportion[,"RHPS4_Input"],
    Induced = rd.proportion[,"Induced_FLAG"] / rd.proportion[,"Induced_Input"])
# ..., that in long form
rd.ratio.long = rd.ratio %>% as_tibble %>%
    rowid_to_column("position") %>%
    gather(key="sample", value="ratio", -1)
# remove one outlier
rd.ratio.long = rd.ratio.long[ rd.ratio.long$position != 3107 , ]
# order the sample name
rd.ratio.long$sample = factor(rd.ratio.long$sample, levels=c("RHPS4", "Induced"))

# attempt at plot like Figure 5C
pdf("Doimo2023_RHPS4_on_chrM.pdf", width=9, height=4)
g = (ggplot(rd.ratio.long, aes(x=position, y=ratio)) 
    + geom_line(aes(color=sample))
    + scale_color_manual(values=c("green", "darkgrey"))
    + geom_hline(yintercept=1, color="black", alpha=0.4) + theme_classic())
plot(g)
dev.off()

# ... similarly, zoomed in to region near a G4
pdf("Doimo2023_RHPS4_on_chrM_zoomed.pdf", width=7, height=6)
if (TRUE) {
xlim = c(6188, 6428)
#     was "xlim(6188, 6388)"
proseq.data = data.frame(position=c(1:16569), PROseq=proseq.mean)
g1 = (ggplot(proseq.data, aes(x=position, y=PROseq))
    + geom_area(color="#0000ff", fill="#0000ff", linewidth=0.25)
    + xlim(xlim[1], xlim[2])
    + theme_classic())
g2 = (ggplot(rd.ratio.long, aes(x=position, y=ratio)) 
    + geom_line(aes(color=sample), size=1.25)
    + scale_color_manual(values=c("green", "darkgrey"))
    + xlim(xlim[1], xlim[2]) + ylim(0, 2)
    + xlab("Position") + ylab("FLAG / Input ratio")
    + geom_hline(yintercept=1, color="black", alpha=0.4)
    + theme_classic()
    + theme(legend.position="bottom"))
g = plot_grid(g1, g2, align="v", axis="b", ncol=1)
plot(g)
}
# take 2, using base R graphics




dev.off()

# write out bigWig files (with other chromosomes set to 0)
# XXX this is a hacky way to get read depth
y = import.bw(paste0(bw.dir, "Induced_Input_1_chrM.bw"), as="RleList")
y = y * 0

# utility to set chrM to some number
rle.with.chrM = function(chrM.read.depth) {
    y1 = y
    y1[["chrM"]] = as.numeric(chrM.read.depth)
    y1
}

dir.create("bw", showWarnings=FALSE)
export.bw(rle.with.chrM(rd.ratio[,"RHPS4"]),
    "bw/RHPS4_ratio_hg38.bw")
export.bw(rle.with.chrM(rd.ratio[,"Induced"]),
    "bw/Induced_ratio_hg38.bw")

