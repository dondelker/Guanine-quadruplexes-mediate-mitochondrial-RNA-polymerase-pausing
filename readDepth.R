# FIXME add GRanges <-> Rle functions, and
# "add strand names" function

#' Gets read depth of an RleList object at selected regions.
#' This is almost like subscripting an RleList with
#' a GRanges object, except that:
#' - the Rles for the minus strand are reversed, and
#' - the names of the result are set from the names
#'   of the region.
#' @param rle an RleList object
#' @param region at which to get coverage of (as a GRanges object)
#' @return a (plain) list of the coverage at "region"
#'   (with the same names, if any)
#' @keywords RleList
#' @import rtracklayer
#' @export
get.rle = function(rle, region) {
  r = as.list(rle[region])
  i = (as.character(strand(region)) == "-")
  r[i] = lapply(r[i], rev)
  names(r) = names(region)
  r
}

#' Converts read depth to a matrix.
#' @param read.depth read depth, as a list of Rle objects (as
#' returned by bw.read.depth or bw.read.depth.stranded)
#' @param num.windows the number of windows
#' @return a matrix, with one row per entry in read.depth
#'   (and the same names), and num.windows columns
#' @keywords RleList
#' @import rtracklayer
#' @export
windowed.matrix = function(read.depth, num.windows) {
  rle.to.matrix(lapply(read.depth, windowize(num.windows)))
}

#' Converts a list of equal-length Rle objects to a matrix.
#' @param x list of Rle objects
#' @return a matrix with one row per element of x
#' @import rtracklayer
rle.to.matrix = function(x) {
  sizes = sapply(x, length)
  stopifnot(all(sizes == sizes[1]))
  r = do.call(rbind, lapply(x, as.numeric))
  rownames(r) = names(x)
  r
}

#' Given a vector, averages it over windows.
#' @param n the number of windows
#' @param x the numbers to windowize
#' @return a vector of length n, giving a windowed
#' average of x (or all NAs, if x is empty).
windowize = function(num.windows) function(x) {
  # convert from Rle to numeric. (This usually seems faster,
  # except for really long regions).
  x = as.numeric(x)
  # if x is empty, return all NAs
  if (length(x)==0)
    return(rep(NA, num.windows))
  # if more windows are requested than there are
  # numbers, expand the numbers in x
  if (num.windows > length(x)) {
    i = trunc( c(0:(num.windows-1)) * (length(x) / num.windows) ) + 1
    return(as.vector(x[i]))
  }
  # compute cumulative sum
  s = c(0, cumsum(x))
  # XXX these indices are a bit opaque
  a = trunc( c(0:(num.windows-1)) * (length(x) / num.windows) ) + 1
  b = trunc( c(1:num.windows) * (length(x) / num.windows) ) + 1
  # return sum-in-a-region, divided by the region size
  as.vector(s[b] - s[a]) / (b - a)
}

