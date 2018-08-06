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

####################################
#FILES:
#
#Input:
#WeissmanTilingScreenGuideSequence.sort.bam
#
#Output:
#modified_sequences_unfiltered.csv
#modified_sequences.csv (filters out unmapped sequences)

####################################
#get genomic sequence from Tiling Screen

bamfile <- '/Users/JYang/Desktop/Stanford/WeissmanTilingScreenGuideSequence.sort.bam'

#read in entire BAM file
bam <- scanBam(bamfile)
#names of the BAM fields
bam_field <- names(bam[[1]])
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

#go through each BAM field and unlist
full_list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

#store as data frame
bam_df <- do.call("DataFrame", full_list)
names(bam_df) <- bam_field

## Read the file into a GRanges object
bgr <- BAM2GRanges(bamfile, 
                   flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE), 
                   verbose = TRUE)
#keep original GRanges object, modify copy
bgr.mod <- bgr

#if (-)ve strand, +8 to start, and +N to end make sequence length = 38
#if (+)ve strand, +8 to end, and +N to start to make sequence length = 38
mod_seq <- function(x){
  for(i in 1:length(x)){
    if(as.logical(strand(x)[i] == "-")){
      start(x)[i] = start(x)[i] - 8
      if(width(x)[i] != 38){
        padding_value = 38 - width(x)[i]
        end(x)[i] = end(x)[i] + padding_value
      }
    }
    if(as.logical(strand(x)[i] == "+")){
      end(x)[i] = end(x)[i] + 8
      if(width(x)[i] != 38){
        padding_value = 38 - width(x)[i]
        start(x)[i] = start(x)[i] - padding_value
      }
    }
  }
  return(x)
}

mod_bgr <- mod_seq(bgr.mod)
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, mod_bgr, as.character=TRUE)
write.csv(seqs, file = "/Users/JYang/Desktop/Stanford/modified_sequences_unfiltered.csv")
seqs.mod <- seqs

#install.packages("gsubfn")
library(gsubfn) 
#check which do not have GG in 32, 33 pos
sum(strapply(seqs, "(..).....$", simplify = TRUE) != "GG")

#indices to remove (line below used in preprocessing_for_models.R)
#flag == 4 is unmapped
remove_idx <- which(bam_df$flag == 4)
remove_idx2 <- which(strapply(seqs, "(..).....$", simplify = TRUE) != "GG")

condition <- lapply(seqs.mod, function(x) strapply(x, "(..).....$", simplify = TRUE) == "GG")
seqs_final <- seqs.mod[unlist(condition)]
write.csv(seqs_final, file = "modified_sequences.csv")

#create sequence logo
#install.packages("ggseqlogo")
require(ggplot2)
require(ggseqlogo)
ggseqlogo(seqs_final, method = 'prob')
