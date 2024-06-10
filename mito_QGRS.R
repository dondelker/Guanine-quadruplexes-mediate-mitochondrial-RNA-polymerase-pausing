# Runs QGRS on sequence (including mitochondrial).

library(rtracklayer)

library(BSgenome.Hsapiens.UCSC.hg19)
genome.hg19 = BSgenome.Hsapiens.UCSC.hg19

seqlevelsStyle(genome.hg19) = "Ensembl"

# where the binary is
QGRS.bin = "../git/qgrs-cpp/bin/qgrs"

# Pattern for 2 or 3 Gs, with a spacer of 1-7 nt
G4.pattern = {
  G = "g{2,3}"
  S = "[ACGT]{1,7}"
  paste0(G,S,G,S,G,S,G)
}

# Runs QGRS on a sequence.
#   genome1: genome object to use for getting sequence
#   g: a GRanges giving the region for which to get sequence
#   Returns: a GRanges object, giving where the matches are
#     (For now, including all matches, which may overlap.)
run.QGRS = function(genome1) function(g) {
  stopifnot(length(g)==1)
  seq1 = getSeq(genome1, g)
  # run QGRS
  seq.temp = tempfile()
  write(as.character(seq1), file=seq.temp)
  r = read.csv(
    pipe(paste(QGRS.bin, "-g g -v -csv -i", seq.temp)),
    as.is=T, strip.white=T)[,1:8]
  file.remove(seq.temp)
  # if result is empty, return an empty GRanges
  if (nrow(r)==0) {
    return(GRanges())
  }
  # create a GRanges of where the matches correspond to
  # (note that T1..T4 are 0-based)
  g4.pos = NULL
  if (as.character(strand(g))=="+")
    g4.pos = GRanges(seqnames(g),
      IRanges(start=start(g)+r$T1, width=nchar(r$SEQ)),
      "+", seqinfo=seqinfo(genome1))
  if (as.character(strand(g))=="-")
    g4.pos = GRanges(seqnames(g),
      IRanges(end=end(g)-r$T1, width=nchar(r$SEQ)),
      "-", seqinfo=seqinfo(genome1))
  # include info about matches
  match.info = r[,c("ID", "TS", "GS", "SEQ")]
  colnames(match.info) = c("ID", "num.tetrads", "G.score", "sequence")
  mcols(g4.pos) = match.info
  # also add on whether this has a spacer
  g4.pos$has.1.to.7.spacer = 1 * grepl(G4.pattern, g4.pos$sequence)
  g4.pos
}

# Wrapper for run.QGRS() which also filters out
# duplicate sites (due to wrapping around on MT).
#   genome1: genome object to use for getting sequence
#   g: a GRanges of sites
#   Returns: result of calling run.QGRS, with duplicate
#     locations on MT removed (based on coordinates)
run.QGRS.mito = function(genome1) function(g) {
  r = run.QGRS(genome1)(g)
  i = start(r) < 1
  r[ i ] = shift(r[i], seqlengths(genome1)[["MT"]])
  r
}

# run on both strands of MT
mt.region = GRanges(
  c("MT:1-16571:+", "MT:1-16571:-", "MT:-200-200:+", "MT:-200-200:-"))
r = lapply(mt.region, run.QGRS.mito(genome.hg19))
# concatenate these
r = r[ sapply(r, length) > 0 ]
r = unlist(GRangesList(r))
# remove duplicates
r = r[ !duplicated(r) ]
# save to a file
write.csv(as.data.frame(r), "seq/proseq2/mitochondria/G4/mito_QGRS.csv", na="")

