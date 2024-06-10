

# Test reading numbers from a bigWig file
test_that("can read some numbers from a bigWig file", {
  expect_equal(2 * 2, 4)
})

test_that("read an unstranded bigWig file (from UCSC; requires Internet)", {
# an arbitrary unstranded track from UCSC
  # a region, on both strands
  g = GRanges(c("chr7:5566779-5570232:+", "chr7:5566779-5570232:-"))
  y = bw.read.depth(
    "http://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeSunyRipSeqGm12878Elavl1SigRep1.bigWig",
    g)
  # these should be the same, except reversed
  expect_equal( y[[1]], rev(y[[2]]) )
})

