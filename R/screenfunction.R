#' Screen a sequence using Lastz
#'
#' You provide a reference and a querry and the function compute the percentage of the reference which is covered by the querry
#'
#' @param reference the reference sequence that you want to screen
#' @param querry the querry sequence
#'
#' @return numeric value of the percentage of the reference sequence which is covered by the querry
#' @import GenomicRanges IRanges Biostrings
#'
#' @export
screenLastz <- function(reference,querry)
{
  try(unlink("temp", recursive = TRUE))
  dir.create('temp',showWarnings = F)
  myarg <- paste0(reference,'[multiple] ',querry,' --ambiguous=iupac --notransition --step=100 --nogapped ‑‑format=rdotplot > temp/result.maf')
  system2(command= 'lastz',args=myarg)

  last <- try(read.table("temp/result.maf"), silent = T)
  if (class(last) == "data.frame")
  {
    start.stop <- as.numeric(na.omit(suppressWarnings(as.numeric(as.character(last$V1)))))
    start <- start.stop[seq(1,length(start.stop),by=2)]
    stop <- start.stop[seq(2,length(start.stop),by=2)]
    GR <- GRanges(seqnames = 'seq',ranges = IRanges(start = start,end = stop))
    GR.disjoin <- disjoin(GR)

    hitlength <- sum(width(GR.disjoin))
    seqlength <- sum(width(readDNAStringSet(reference)))

    percentage <- 100*round(hitlength/seqlength,3)
  }
  else (percentage <- 0)
  unlink('temp',recursive = T)
  return(percentage)
}



#' Screen a sequence using Blast
#'
#' You provide a reference and a querry and the function compute the percentage of the reference which is covered by the querry
#'
#' @param reference the reference sequence that you want to screen. Fasta file in one or severa sequences
#' @param querry the querry sequence. Fasta file in one or severa sequences
#'
#' @return numeric value of the percentage of the reference sequence which is covered by the querry
#' @import GenomicRanges IRanges Biostrings
#'
#' @export
screenBlast <-function (reference, querry,min.pc.ident)
{
  try(unlink("temp", recursive = TRUE))
  dir.create("temp")
  dir.create("temp/dbblast")
  myarg <- paste0("-in ", reference, " -out temp/dbblast/db -dbtype nucl")
  system2(command = "makeblastdb", args = myarg, stdout = F)
  myarg <- paste0("-query ", querry, " -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt \"7 qacc bitscore qlen length pident qstart qend sacc sstart send \"")
  system2(command = "blastn", args = myarg)
  blast <- try(read.table("temp/blast.txt", comment.char = "#"), silent = T)
  if (class(blast) == "data.frame")
  {
    colnames(blast) <- c("querry.access", "bitscore", "querry.length", "alignment.lenght", "pc.ident.",
                         "querry.start", "querry.end", "subject.access", "subject.start", "subject.end")
    blast <- blast[blast$pc.ident.>min.pc.ident,]
    start <- blast$subject.start
    end <- blast$subject.end
    new.start <- start
    new.end <- end
    new.start[start>end] <- end[start>end]
    new.end[start>end] <- start[start>end]
    data.frame(start,end,new.start,new.end)
    GR <- GRanges(seqnames = blast$subject.access,ranges = IRanges(start =new.start ,end = new.end))
    GR.disjoin <- disjoin(GR)
    hitlength <- sum(width(GR.disjoin))
    seqlength <- sum(width(readDNAStringSet(reference)))
    percentage <- 100*round(hitlength/seqlength,3)
  }
  else{percentage <- 0}
  unlink('temp',recursive = T)
  return(percentage)
}

