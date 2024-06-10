# Computes BG4 difference.

library(parallel)
library(rtracklayer)

options(mc.cores=10)

# directory containing the bigWig files
data.dir = "../pub/hg19/bigWig_RPM/G4/Mao2018/"
# names of the files to read in
bw.file = list.files(data.dir, pattern=".clean.bw")

# read in all the bigWig files
if (!exists("bw")) {
  # first, read in the files
  bw = mclapply(bw.file, function(f)
    import.bw(paste0(data.dir, f), as="RleList"))
  # simplify file names
  sample.name = gsub(".*K562_", "", bw.file)
  sample.name = gsub("(\\.rmdup)?\\.clean\\.bw", "", sample.name)
  names(bw) = sample.name
}

# average BG4 IP, and input ("IP")
bg4 = (bw[["asynch_a_701_504"]]
  + bw[["asynch_c_703_504"]]
  + bw[["P9_Async_a_701_517"]]
  + bw[["P9_Async_b_701_502"]]
  + bw[["P9_Async_c_701_503"]]) / 5
input = (bw[["asynch_IP_704_504"]]
  + bw[["P9_Async_IP_702_504"]]) / 2
# compute difference
bg4.minus.input = bg4 - input
# write out difference (in same directory)
export.bw(bg4.minus.input, paste0(data.dir, "BG4_minus_input.bw"))

