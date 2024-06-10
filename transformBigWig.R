# Utilities for transforming a track in various ways.

library(parallel)

source("R/file/bigWig.R")

# Applies a function to a circular sequence (e.g. chrM).
#   f: the function
#   margin: the amount to add on to each end of the numbers
#     (assumes that the input is longer than this)
#   x: the input, as an RleList
#   Returns: f(x), assuming (for the sake of windowed
#     calculations) that x "wraps around"
apply.circular = function(f, margin=1000) function(x) {
  n = length(x)
  # pad x
  x1 = c(x[ (n+1-margin) : n ], x, x[ 1 : margin ])
  # compute f on the padded vector
  y1 = f(x1)
  # return f of just the central chunk
  y1[ (margin+1) : (n+margin) ]
}

# Applies a function to (stranded) bigWig coverage files.
#   f: function to apply (to each chromosome of an RleList)
#   input.dir: where to read files from
#   output.dir: where to write transformed files to
#   flip.minus: if TRUE, flip the minus strand, both when
#     reading and writing files
#   Side effects: writes files with the same name (skipping
#     files which are already present)
write.transformed.big.wig.stranded =
    function(f, input.dir, output.dir, flip.minus=TRUE) {
  system(paste0("mkdir -p ", output.dir))
  bw.file = get.stranded.bw.experiments(input.dir)
  for(bw in bw.file) {
    input.bw.base = paste0(input.dir, bw)
    output.bw.base = paste0(output.dir, bw)
    if (file.exists(paste0(output.bw.base, ".plus.bw")) &&
        file.exists(paste0(output.bw.base, ".minus.bw")))
      next
    cat(output.bw.base, "\n")
    x = import.bw.stranded(input.bw.base)
    if (flip.minus)
      x = rlelist.flip.minus(x)
    y = mclapply(x, f, mc.cores=62, mc.preschedule=F)
    y = as(y, "RleList")
    if (flip.minus)
      y = rlelist.flip.minus(y)
    export.bw.stranded(y, output.bw.base)
  }
}

# Writes a bigWig file, which is a function of some other files.
# (Computes for all chromosomes in parallel, based on
# "mc.cores" option).
#   output.file: name of output file
#   input.files: names of input files
#   f: function to apply
#   Side effects: reads in input files, computes f of them, and
#     writes them out as a bigWig file. (If output.file already
#     exists, does nothing).
write.transformed.bw = function(output.file, input.files, f) {
  if (file.exists(output.file))
    return()
  cat("writing ", output.file, "\n")
  # read in data
  x = mclapply(input.files,
    function(bw) import.bw(bw, as="RleList"))
  # get list of chromosomes present in all of the input files
  chr = names(x[[1]])
  for (i in 1:length(x))
    chr = intersect(chr, names(x[[i]]))
  chr = sort(chr)
  # compute f(x) (restricting to common chromosomes)
  # XXX this is sort of inelegant
  if (length(x) == 1)
    y = mcmapply(f, x[[1]][chr])
  if (length(x) == 2)
    y = mcmapply(f, x[[1]][chr], x[[2]][chr])
  if (length(x) == 3)
    y = mcmapply(f, x[[1]][chr], x[[2]][chr], x[[3]][chr])
  if (length(x) == 4)
    y = mcmapply(f,
      x[[1]][chr], x[[2]][chr], x[[3]][chr], x[[4]][chr])
  # write out y
  y = as(y, "RleList")
  export.bw(y, output.file)
}

# Like write.transformed.bw, but writes a stranded file.
#   output.base: base name of output bigWig files
#   input.files: base names of input files.
#     (??? currently these are assumed to be stranded;
#     but these could be a mix of stranded and unstranded,
#     determined by file name)
#   f: function to apply
#   Side effects: writes out a stranded bigWig file
write.transformed.bw.stranded =
    function(output.base, input.files, f) {
  write.transformed.bw(paste0(output.base, ".plus.bw"),
    paste0(input.files, ".plus.bw"), f)
  write.transformed.bw(paste0(output.base, ".minus.bw"),
    paste0(input.files, ".minus.bw"), f)
}

# Computes a z-score.
#   window.size: size of window for statistics (this uses
#     S4Vectors::runmean()'s definition, which is the total width
#     on both sides of a location; thus "window.size=201" means
#     a site, plus 100 bp on each side)
#   min.coverage: only include sites with coverage above this cutoff
#   min.z: only include z-scores above this number
#   Returns: function to compute Z
z.score = function(window.size=201, min.coverage=1, min.z=3) function(x) {
  # get windowed mean and variance
  m = S4Vectors::runmean(x, window.size, endrule="constant")
  v = S4Vectors::runmean((x-m)^2, window.size, endrule="constant")
  # compute z-score
  z = (x - m) / sqrt(v)
  # convert NaNs to 0
  runValue(z)[ is.nan(runValue(z)) ] = 0
  # mask out cases with coverage or Z too low
  z[ (x < min.coverage) | (z < min.z) ] = 0
  z
}

# Computes coverage relative to local background.
#   window.size: size of window for statistics
#   min.coverage: only include sites with coverage above this cutoff
#   Returns: function to compute relative coverage
rel.coverage = function(window.width, min.coverage=1) function(x) {
  # get windowed mean
  m = S4Vectors::runmean(x, window.size, endrule="constant")
  # compute ratio
  r = x / m
  # convert NaNs to 0
  runValue(r)[ is.nan(runValue(r)) ] = 0
  # mask out cases with coverage too low
  r[ x < min.coverage ] = 0
  r
}

