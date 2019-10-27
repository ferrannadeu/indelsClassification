indelsClassification <- function(mat, classColumn=NULL, subClassColumn=NULL, repeatsColumns=NULL, title="Classification of small indels", subtitle=FALSE, yScale="percentage"){
  
  if(is.null(classColumn) | is.null(subClassColumn) | is.null(repeatsColumns)){
    if(nrow(mat) > 0) {
      mat <- mat[,1:4]
      colnames(mat) <- c("CHR", "POSITION", "REF", "ALT")
      if(!startsWith(as.character(mat$CHR[1]), "chr")){ mat$CHR <- paste0("chr", mat$CHR) }
      
      if(any(nchar(mat$REF) == 1 & nchar(mat$ALT) == 1)){ stop("SNVs in input are not allowed!") }
      if(any(nchar(mat$REF) == nchar(mat$ALT))){ stop("MNVs in input are not allowed!") }
      
      mat$context5p <- NA
      mat$context3p <- NA
      mat$indelClass <- NA
      mat$indelSubClass <- NA
      mat$indelRepeats <- NA
      
      # context at 5' of POSITION (POSITION incldued)
      mx <- apply(mat, 1, function(x){ max(c(nchar(x[3]), nchar(x[4])))-1})
      gr1 <- GRanges(mat$CHR, IRanges(start=as.numeric(mat$POSITION)-mx, end=as.numeric(mat$POSITION)), strand = "+")
      context <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gr1)
      context <- as.data.frame(context)$x
      mat$context5p <- context
      
      # context at 3' of POSITION (POSITION incldued)
      mx <- apply(mat, 1, function(x){ max(c(nchar(x[3]), nchar(x[4])))})
      gr1 <- GRanges(mat$CHR, IRanges(start=as.numeric(mat$POSITION), end=as.numeric(mat$POSITION)+mx*6), strand = "+")
      context <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gr1)
      context <- as.data.frame(context)$x
      mat$context3p <- context
      
      for(i in 1:nrow(mat)){
        chrom <- mat$CHR[i]
        position <- mat$POSITION[i]
        ref <- mat$REF[i]
        alt <- mat$ALT[i]
        
        # complex indels
        if(nchar(ref) > 1 && nchar(alt) > 1){
          mat$indelClass[i] <- "complex indel"
          mat$indelSubClass[i] <- "."
          mat$indelRepeats[i] <- "."
        }
        
        # insertions
        else if(nchar(ref) < nchar(alt)){ 
          ins <- substr(alt, 2, nchar(alt))
          repeats <- NA
          context3p <- substring(mat$context3p[i], 2)
          for(c in 5:0){ 
            if(!isEmpty(grep(pattern = paste0("^(", ins,"){", c, "}"), x = context3p, perl = T))){
              if(grep(pattern = paste0("^(", ins,"){", c, "}"), x = context3p, perl = T)[1]==1){ 
                repeats <- c; break 
              }
            }
          }
          if(nchar(alt)-nchar(ref)==1){
            cl <- "1bp insertion"
            sc <- ins
            if(ins == "G"){ sc <- "C:G" }
            else if(ins == "C"){ sc <- "C:G" }
            else if(ins == "T"){ sc <- "T:A" }
            else if(ins == "A"){ sc <- "T:A" }
            if(repeats == 5){ repeats <- "5+" }
          }
          else{
            cl <- ">1bp insertion at repeats"
            sc <- nchar(alt)-nchar(ref)
            if(sc >= 5) { sc <- "5+" }
            if(repeats == 5){ repeats <- "5+" }
          }
          mat$indelClass[i] <- cl
          mat$indelSubClass[i] <- sc
          mat$indelRepeats[i] <- repeats
        }
        
        # deletions
        else{
          del <- substr(ref, 2, nchar(ref))
          repeats <- NA
          context3p <- substring(mat$context3p[i], 2)
          for(c in 6:1){ 
            if(!isEmpty(grep(pattern = paste0("^(", del,"){", c, "}"), x = context3p, perl = T))){
              if(grep(pattern = paste0("^(", del,"){", c, "}"), x = context3p, perl = T)[1]==1){ 
                repeats <- c; break 
              }
            }
          }
          
          if(nchar(ref)-nchar(alt) == 1){ # 1bp deletion
            cl <- "1bp deletion"
            sc <- del
            if(del == "G"){ sc <- "C:G" }
            else if(del == "C"){ sc <- "C:G" }
            else if(del == "A"){ sc <- "T:A" }
            else if(del == "T"){ sc <- "T:A" }
            if(repeats == 6){ repeats <- "6+" }
          }
          
          else{ # >1bp deletion at repeats or deletions with microhomology
            
            if(repeats > 1){ # deletion at repeats (2 or more repeats)
              cl <- ">1bp deletion at repeats"
              sc <- nchar(ref)-nchar(alt)
              if(sc >= 5) { sc <- "5+" }
              if(repeats == 6){ repeats <- "6+" }
            }
            
            else{ # if repeats == 1, deletion with repeat size of 1 or deletions with microhomology?
              # check homology at 3':
              homology3p = NA
              context3p <- substring(mat$context3p[i], 2)
              for(c in (nchar(del)-1):1){
                delNew <- paste0(del, substring(del, 1, c))
                if(!isEmpty(grep(pattern = paste0("^(", delNew,"){", 1, "}"), x = context3p, perl = T))){
                  if(grep(pattern = paste0("^(", delNew,"){", 1, "}"), x = context3p, perl = T)[1]==1){ 
                    homology3p <- c; break 
                  }
                }
              }
              if(is.na(homology3p)){ homology3p <- 0 }
              
              # check homology at 5':
              homology5p <- NA
              context5p <- mat$context5p[i]
              for(c in 2:nchar(del)){
                delNew <- substring(del, c, nchar(del))
                if(!isEmpty(grep(pattern = paste0("(", delNew,"){", 1, "}$"), x = context5p, perl = T))){
                  if(grep(pattern = paste0("(", delNew,"){", 1, "}$"), x = context5p, perl = T)[1]==1){ 
                    homology5p <- nchar(delNew); break 
                  }
                }
              }
              if(is.na(homology5p)){ homology5p <- 0 }
              
              # get maximum homology
              if(homology5p > homology3p){ homology <- homology5p }
              else { homology <- homology3p }
              
              
              if(homology > 0){ # if homology, deletion with microhomology
                cl <- "deletion with microhomology"
                sc <- nchar(ref)-nchar(alt)
                if(sc >= 5) { sc <- "5+" }
                repeats <- homology
                if(repeats >=5){ repeats <- "5+" }
              }
              else{ # else, >1bp deletion at repeats, 1 repeat unit
                cl <- ">1bp deletion at repeats"
                sc <- nchar(ref)-nchar(alt)
                if(sc >= 5) { sc <- "5+" }
                if(repeats == 6){ repeats <- "6+" }
              }
            }
          }
          mat$indelClass[i] <- cl
          mat$indelSubClass[i] <- sc
          mat$indelRepeats[i] <- repeats
        }
      }
    }
  }
  
  # if class, subclass and repeats columns are already defined
  else{
    if(nrow(mat) > 0){
      mat <- as.data.frame(mat)
      mat1 <- mat[,1:4]
      mat1 <- cbind(mat1, mat[,c(classColumn, subClassColumn, repeatsColumns)])
      mat <- mat1
      colnames(mat) <- c("CHR", "POSITION", "REF", "ALT", "indelClass", "indelSubClass", "indelRepeats")
    }
  }
  
  if(nrow(mat) > 0){
    classes <- data.frame(matrix(c(rep("1bp deletion", 12),
                                   rep("1bp insertion", 12),
                                   rep(">1bp deletion at repeats", 24),
                                   rep(">1bp insertion at repeats", 24),
                                   rep("deletion with microhomology", 11),
                                   "complex indel",
                                   
                                   rep("C:G", 6), rep("T:A", 6),
                                   rep("C:G", 6), rep("T:A", 6),
                                   rep( c(2:4, "5+"), each=6),
                                   rep( c(2:4, "5+"), each=6),
                                   c(2, rep(3, 2), rep(4,3), rep("5+",5)),
                                   ".",
                                   
                                   rep( c(1:5, "6+"), 2),
                                   rep( c(0:4, "5+"), 2),
                                   rep( c(1:5, "6+"), 4),
                                   rep( c(0:4, "5+"), 4),
                                   c(1, 1:2, 1:3, 1:4, "5+"),
                                   "."
    ), ncol = 3, byrow = F))
    colnames(classes) <- c("class", "subclass", "repeats")
    classes$counts <- 0
    
    for(i in 1:nrow(classes)){
      classes[i, "counts"] <- nrow(mat[mat$indelClass == classes$class[i] & mat$indelSubClass == classes$subclass[i] & mat$indelRepeats == classes$repeats[i],])
    }
    classes$type <- paste(classes$class, classes$subclass, sep="_")
    classes$class <- factor(classes$class, levels = unique(classes$class))
    classes$subclass <- factor(classes$subclass, levels = unique(classes$subclass))
    classes$type <- factor(classes$type, levels = unique(classes$type))
    
    tit <- title
    if(subtitle){ subtit <-  paste0("Indels input = ", nrow(mat), " | Indels annotated = ", sum(classes$counts)) }
    else{ subtit <- NULL }
    if(yScale== "count"){ yLab <- "Count" }
    else{ 
      yLab <- "Percentage (%)"
      classes$counts <- classes$counts/sum(classes$counts)*100
    }
    
    cls <- c(rep("#FDBE6F", 6), rep("#FF8001", 6), rep("#B0DD8B", 6), rep("#36A12E", 6), 
             rep("#FDCAB5", 6), rep("#FC8A6A", 6), rep("#F14432", 6), rep("#BC141A", 6),
             rep("#D0E1F2", 6), rep("#94C4DF", 6), rep("#4A98C9", 6), rep("#1764AB", 6),
             rep("#E2E2EF", 1), rep("#B6B6D8", 2), rep("#8683BD", 3), rep("#61409B", 5),
             rep("#262626", 1))
    
    # pdf(NULL)
    dev.control(displaylist="enable")
    
    par(mar=c(4.1,1.1,4.1,0.1))
    mp <- barplot(classes$counts, 
                  las=1, axes=F,
                  col = cls, space=0.75)
    title(main=tit, line = 3, adj=0, cex.main=0.8)
    title(sub=subtit, cex.sub=0.7, line=2.35, adj=0)
    axis(2, line = -2, las=2, cex.axis=0.6)
    mtext(yLab, side = 2, line = 0.25, cex=0.66)
    axis(1, tick = F, at = as.vector(mp), labels = c(as.character(classes$repeats)[1:83], ""), cex.axis=0.55, line = -0.5)
    rect(xleft = as.vector(mp)[c(seq(1,67,6), 73, 74, 76, 79, 84)]-0.5, 
         ybottom = rep(0-max(classes$counts)*0.05,17), 
         xright = as.vector(mp)[c(seq(6,72,6), 73, 75, 78, 83, 84)]+0.5,
         ytop =  rep(0-max(classes$counts)*0.02,17), 
         col = unique(cls), xpd=T)
    rect(xleft = as.vector(mp)[c(seq(1,67,6), 73, 74, 76, 79, 84)]-0.5, 
         ybottom = rep(max(classes$counts)*1.035,17), 
         xright = as.vector(mp)[c(seq(6,72,6), 73, 75, 78, 83, 84)]+0.5, 
         ytop =  rep(max(classes$counts)*1.085,17), 
         col = unique(cls), xpd=T)
    text(x = c((as.vector(mp)[1:83]+(as.vector(mp)[2:84]- as.vector(mp)[1:83])/2)[seq(3,73,6)],
               as.vector(mp)[73], as.vector(mp)[74]+0.875, as.vector(mp)[c(77, 81)]),
         y=max(classes$counts)*1.06, 
         c("C", "T", "C", "T", rep(c(2:4, "5+"),3)), xpd=T, cex = 0.6)
    text(x=c((as.vector(mp)[1:83]+(as.vector(mp)[2:84]- as.vector(mp)[1:83])/2)[c(6,18,36,60)],
             as.vector(mp)[c(78, 84)]),
         y=max(classes$counts)*1.15, 
         c("1bp deletion", "1bp insertion", ">1bp deletion at repeats\n(deletion length)", 
           ">1bp insertion at repeats\n(insertion length)", 
           "Deletion with microhomology\n(deletion length)", "Complex\nindels"),
         xpd=T, cex = 0.64)
    text(x=c((as.vector(mp)[1:83]+(as.vector(mp)[2:84]- as.vector(mp)[1:83])/2)[c(6,18,36,60)],
             as.vector(mp)[78]),
         y=0-max(classes$counts)*0.16, 
         c("Homopolymer length", "Homopolymer length", "Number of repeat units", 
           "Number of repeat units", "Microhomology length"),
         xpd=T, cex = 0.62)
    
    p <- recordPlot()
    invisible(dev.off())
  }
  
  else{ p <- NULL }
  return(list(mat, p))
}
