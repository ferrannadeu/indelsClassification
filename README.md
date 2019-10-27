# indelsClassification

#### R function for the classification of small insertions and deletions for mutational signature analysis

The indelsClassification R function performs the classification of small insertions and deletions (indels) in 84 classes following the criteria described in [The Repertoire of Mutational Signatures in Human Cancer](https://www.biorxiv.org/content/early/2018/05/15/322859) (Alexandrov *et al.*, bioRxiv 2018). An extra indel class has been added to account for "complex insertions/deletions" (see the *PCAWG7_indel_classification_2017_12_08_modified_by_FNadeu.xlsx* document for details).

![alt text](https://github.com/ferrannadeu/indelsClassification/blob/master/indelsClassification.jpeg "indelsClassification output")


### Requirements

The following R libraries are required:

* GenomicRanges
* BSgenome
* BSgenome.Hsapiens.UCSC.hg19
* Rsamtools


### Input

The input matrix/data frame must have at least these four columns (header names can be anything):

1. Chromosome: both annotations ("1" or "chr1") are allowed
2. Position of the indel
3. Reference sequence of the indel
4. Alternate sequence of the indel

Note that indels must be normalized, and that SNVs and/or MNVs are not allowed.


### Outputs

The indelsClassification function returns a list of two objects:

* Matrix with the 5' and 3' context of the indels and their classification (note that the nucleotide at 'position' is included both in the 5’ and 3’ context sequences)
* Graphical representation of the classification


### Run instructions

For a detailed explanation see the *indelsClassification.html* R Markdown document.

Briefly, in your R script:

```r
# load libraries
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Rsamtools)

# source function
source("indelsClassification.R")

# read input file with chr, position, ref, alt for some indels
INDELS <- read.table("indelsClassification_toyExample.tsv", sep = "\t", header = T, stringsAsFactors = F) 

# set display method
#pdf(file=paste0('/path/to/image.pdf'), width = 15, height = 4)
#png(filename=paste0('/path/to/image.png'), width = 1500, height = 400, type=c("cairo-png"))
pdf(NULL) # screen output

# run the function
out <- indelsClassification(mat = INDELS)
while (!is.null(dev.list()))  dev.off()

# check outputs
head(out[[1]])
out[[2]]

# grid::grid.newpage() # run this in your R script to open a new page for the next plot
```


### Note

The script seems to work well and we constantly check its output when running it on new indels. Potential bugs will be corrected and improvements may be added. Feel free to use it and share it!


### Bugs, comments and improvements

Bugs, comments and improvements can be send to *nadeu@clinic.cat*. They will be very much appreciated!


### Acknowledgments

I thank Francesco Maura for the time he invested in reviewing the output of this function as well as for his comments and suggestions.
