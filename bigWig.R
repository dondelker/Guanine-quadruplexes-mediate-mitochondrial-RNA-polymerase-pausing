
#' Gets read depth as a list of numeric vectors.
#'
#' Consider using bw.read.depth.1, which should give the
#' same results, and is much faster.
#' @param bigWigName name of a bigWig file (including the
#' ".bw" suffix)
#' @param region the regions, as GRanges objects
#' @return a list of read depth
#' @keywords bigWig
#' @examples
#' \dontrun{bw.read.depth("file.bw", GRanges("1:1000-2000:+"))}
#' @import rtracklayer
#' @export
bw.read.depth = function(bigWigName, region) {
  # ??? also allow passing in a BigWigFile object?)
  bigWig = BigWigFile(bigWigName)
  # get read depth as an RleList object (ignoring warning that
  # "'which' contains seqlevels not known to BigWig file")
  rle = suppressWarnings(
    import(bigWig, which=region, as="RleList"))
  get.rle(rle, region)
}

#' Gets read depth as a list of numeric vectors, for a pair of
#' stranded bigWig files.
#'
#' Consider using bw.read.depth.1, which should give the
#' same results, and is much faster.
#' @param bigWigName name of the plus-strand bigWig file
#' (e.g. including the ".plus.bw" suffix; other suffixes
#' may be supported later)
#' @param region the regions, as GRanges objects
#' @param abs.value if TRUE, then this returns the absolute value
#'   of the numbers (this is to avoid the problem that often the
#'   the - strand is stored negated, for viewing in e.g. IGV).
#' @return a list of read depth
#' @keywords bigWig
#' @examples
#' @import rtracklayer
#' @export
bw.read.depth.stranded = function(bigWigName, region, abs.value=T) {
  # FIXME
  # - support other popular suffixes?
  # - allow specifying the minus file name as well?
  #   (which hopefully will prevent having to add too many
  #   more suffixes)
  bigWigPlusFile = bigWigName
  bigWigMinusFile = sub("\\.plus\\.bw$", ".minus.bw", bigWigPlusFile)
  r = as.list(rep(NA, length(region)))
  plus = (as.character(strand(region)) != "-")
  r[ plus ] = bw.read.depth(bigWigPlusFile, region[plus])
  r[ !plus ] = bw.read.depth(bigWigMinusFile, region[!plus])
  if (abs.value)
    r = lapply(r, abs)
  names(r) = names(region)
  r
}


#' Gets read depth as a list of numeric vectors.
#' @param bigWigName name of a bigWig file (including the
#' ".bw" suffix)
#' @param region the regions, as GRanges objects
#' @return a list of read depth
#' @keywords bigWig
#' @examples
#' \dontrun{bw.read.depth("file.bw", GRanges("1:1000-2000:+"))}
#' @import rtracklayer
#' @export
bw.read.depth.1 = function(bigWigName, region) {
  # ??? also allow passing in a BigWigFile object?)
  bigWig = BigWigFile(bigWigName)
  # get read depth as a NumericList object (ignoring warning if
  # "'which' contains seqlevels not known to BigWig file")
  y = suppressWarnings(
    import(bw, which=region, as="NumericList"))
  y = revElements(y, strand(region)=="-")
  y
}

#' Gets read depth as a list of numeric vectors, for a pair of
#' stranded bigWig files.
#' @param bigWigName name of the plus-strand bigWig file
#' (e.g. including the ".plus.bw" suffix; other suffixes
#' may be supported later)
#' @param region the regions, as GRanges objects
#' @param abs.value if TRUE, then this returns the absolute value
#'   of the numbers (this is to avoid the problem that often the
#'   the - strand is stored negated, for viewing in e.g. IGV).
#' @return a list of read depth
#' @keywords bigWig
#' @examples
#' @import rtracklayer
#' @export
bw.read.depth.stranded.1 = function(bigWigName, region, abs.value=T) {
  # FIXME
  # - support other popular suffixes?
  # - allow specifying the minus file name as well?
  #   (which hopefully will prevent having to add too many
  #   more suffixes)
  bigWigPlusFile = bigWigName
  bigWigMinusFile = sub("\\.plus\\.bw$", ".minus.bw", bigWigPlusFile)
  r = as.list(rep(NA, length(region)))
  plus = (as.character(strand(region)) != "-")
  r[ plus ] = bw.read.depth.1(bigWigPlusFile, region[plus])
  r[ !plus ] = bw.read.depth.1(bigWigMinusFile, region[!plus])
  if (abs.value)
    r = lapply(r, abs)
  names(r) = names(region)
  r
}


#' Reads in read depth for a stranded pair of bigWig files,
#' as an RleList with ".plus" or ".minus" at the ends of
#' chromosome names (or sequence names).
#' This is intended to simplify genome-wide math (e.g.,
#' taking the ratio of numbers across the genome).
#' @param bigWigPlusFile name of the plus-strand bigWig file
#' (e.g. including the ".plus.bw" suffix; other suffixes
#' may be supported later)
#' @param bigWigMinusFile name of the minus-strand bigWig file
#'   (if not present, this will be guessed)
#' @param abs.value if True, this takes tha absolute value
#'   of all of the numbers
#' @return "flat" RleList, with ".plus" and ".minus"
#'   suffixed onto chromosome names. (The list will be
#'   sorted by these modified chromosome names).
#' @keywords bigWig
#' @examples
#' @import rtracklayer
#' @export
#' @seealso import.bw to import an unstranded bigWig file
import.bw.stranded = function(bigWigPlusFile, bigWigMinusFile=NULL, abs.value=TRUE) {
  # FIXME add other file suffixes (e.g. ".for.bw" and ".rev.bw"
  # from deepTools)
  if (is.null(bigWigMinusFile))
    bigWigMinusFile = sub("\\.plus\\.bw$", ".minus.bw", bigWigPlusFile)
  # read in coverage
  read.depth.plus = import(bigWigPlusFile, as="RleList")
  read.depth.minus = import(bigWigMinusFile, as="RleList")
  # alter names of these
  names(read.depth.plus) = paste0(names(read.depth.plus), ".plus")
  names(read.depth.minus) = paste0(names(read.depth.minus), ".minus")
  r = c(read.depth.plus, read.depth.minus)
  if (abs.value)
    r = abs(r)
  r[ sort(names(r)) ]
}

#' Writes out coverage from an RleList corresponding to a
#' stranded dataset, split into two bigWig files with names
#' ending in ".plus.bw" and ".minus.bw". (This is the opposite
#' of import.bw.stranded).
#' @param x an RleList object, with strands stored in chromosomes
#' with ".plus" and ".minus" as suffixes.
#' @param base.name base name for output files (which will be
#'   suffixed with ".plus.bw" and ".minus.bw", respectively).
#' @keywords bigWig
#' @examples
#' @import rtracklayer
#' @export
export.bw.stranded = function(x, base.name) {
  x.plus = x[ grepl("\\.plus$", names(x)) ]
  names(x.plus) = gsub("\\.plus", "", names(x.plus))
  export.bw(x.plus, paste0(base.name, ".plus.bw"))

  x.minus = x[ grepl("\\.minus$", names(x)) ]
  names(x.minus) = gsub("\\.minus", "", names(x.minus))
  export.bw(x.minus, paste0(base.name, ".minus.bw"))
}

#' Negates entries of an RleList which end in ".minus"
#' (and are presumably on the "-" strand).
#' This is to deal with pairs of bigWig files for stranded
#' data, in which the "-" strand is negated for viewing in IGV.
#' (This is intended to work with import.bw and export.bw).
#' @param x an RleList, whose names are suffixed with
#'   ".plus" or ".minus"
#' @keywords bigWig
#' @examples
#' @import rtracklayer
#' @export
rlelist.flip.minus = function(x) {
  i = grepl("\\.minus$", names(x))
  x[i] = x[i] * -1
  x
}

#' Sorts names of an RleList to match the ordering of, e.g.,
#'   a GRanges object.
#' @param x an RleList (with names suffixed with ".plus" and ".minus")
#' @param chrom.order: the order in which to put the chromosomes (as a
#'   vector of strings)
#' @return the RleList, sorted (first by chrom.order, then with ".plus"
#'   before ".minus")
#' @export
sort.rlelist = function(x, chrom.order) {
  # create the ordering
  a = paste0(rep(chrom.order, each=2), c(".plus",".minus"))
  # only include things in the RleList
  a1 = a[ a %in% names(x) ]
  x[ a1 ]
}

