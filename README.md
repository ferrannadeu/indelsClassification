# indelsClassification

#### R function for the classification of small insertions and deletions for mutational signature analysis

The indelsClassification R function performs the classification of small insertions and deletions (indels) in 84 classes following the criteria described in [The Repertoire of Mutational Signatures in Human Cancer](https://www.biorxiv.org/content/early/2018/05/15/322859) (Alexandrov *et al.*, bioRxiv 2018). An extra indel class has been added to account for "complex insertions/deletions" (see the *PCAWG7_indel_classification_2017_12_08_modified_by_FNadeu.xlsx* document for details).


### Libraries

The following R libraries are required:

* GenomicRanges
* BSgenome
* BSgenome.Hsapiens.UCSC.hg19
* Rsamtools

### Input

The input matrix/data frame must have at least these four columns (header names can be anything):

1. Chromosome: both annotations "1" and "chr1" are allowed
2. Position of the indel
3. Reference sequence of the indel
4. Alternate sequence of the indel

Note that indels must be normalized, and that SNVs and/or MNPs are not allowed.

### Outputs

The indelsClassification function returns a list of two objects:

* Matrix with the 5' and 3' context of the indels and their classification (note that the nucleotide at 'position' is included both in the 5’ and 3’ context sequences)
* Graphical representation of the classification

### Run instructions

For a detailed explanation see the *indelsClassification.html* Rmarkdown document.

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

# run the function
out <- indelsClassification(mat = INDELS)

# check outputs
head(out[[1]])
out[[2]]
```

### Note

The script seems to work well, but we constantly check it when running it on new samples/indels. Potential bugs will be corrected and improvements may be added. Feel free to use it and share it!

### Bugs, comments or improvements

Bugs, comments and improvements might be send to *nadeu@clinic.ub.es*. They will be very much appreciated!

### Aknowledgments

I appreciate to Francesco Maura the time invested in reviewing the output of this function as well as for appropriate comments and suggestions.









