source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
biocLite("GenomicRanges")
biocLite("Rsamtools")
biocLite("GenomicFeatures")
biocLite("chromstaR")
biocLite("truncnorm")
biocLite("Repitools")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(Rsamtools)
library(GenomicFeatures)
library(biomaRt)
library(GenomicRanges)
library(Repitools)
library(chromstaR)
library(BSgenome.Hsapiens.UCSC.hg19)

bamfile <- '/Users/JYang/Desktop/Stanford/WeissmanTilingScreenGuideSequence.sort.bam'

#read in entire BAM file
bam <- scanBam(bamfile)
#names of the BAM fields
names(bam[[1]])
#distribution of BAM flags
table(bam[[1]]$flag)
#function for collapsing the list of lists into a single list
#as per the Rsamtools vignette
.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#store names of BAM fields
bam_field <- names(bam[[1]])

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

## Read the file into a GRanges object
reads <- readBamFileAsGRanges(bamfile, pairedEndReads=FALSE, remove.duplicate.reads=TRUE)
#print(head(reads))

#bgr <- BAM2GRanges(bamfile ,
#                   what = c("qname",  "flag",   "rname",  "strand", "pos",    "qwidth", "mapq",   "cigar",  "mrnm",   "mpos",  
#                            "isize",  "seq",    "qual"),
#                   flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE),
#                   verbose = TRUE)

start(reads) <- start(reads) - 5
end(reads) <- end(reads) + 5

print(head(reads))

seqs = getSeq(BSgenome.Hsapiens.UCSC.hg19, reads, as.character=TRUE) 
head(seqs)

write.csv(seqs, file = "/Users/JYang/Desktop/Stanford/seqs.csv")
