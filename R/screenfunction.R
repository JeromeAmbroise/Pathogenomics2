#' Screen a sequence
#'
#' You provide a reference and a querry and the function compute the percentage of the reference which is covered by the querry
#'
#' @param reference the reference sequence that you want to screen
#' @param querry the querry sequence
#'
#' @return numeric value of the percentage of the reference sequence which is covered by the querry
#'
#' @export
screenLastz <- function(reference,querry)
{
  try(unlink("temp", recursive = TRUE))
  dir.create('temp',showWarnings = F)
  myarg <- paste0(reference,' ',querry,' --ambiguous=iupac --notransition --step=100 --nogapped ‑‑format=rdotplot > temp/result.maf')
  system2(command='/home/jerome/Documents/0-installation-files/lastz/src/lastz',args=myarg)

  last <- try(read.table("temp/result.maf"), silent = T)
  if (class(last) == "data.frame")
  {
    start.stop <- as.numeric(na.omit(suppressWarnings(as.numeric(as.character(last$V1)))))
    start <- start.stop[seq(1,length(start.stop),by=2)]
    stop <- start.stop[seq(2,length(start.stop),by=2)]
    GR <- GRanges(seqnames = 'seq',ranges = IRanges(start = start,end = stop))
    GR.disjoin <- disjoin(GR)

    hitlength <- sum(width(GR.disjoin))
    seqlength <- width(readDNAStringSet(reference))
    percentage <- 100*round(hitlength/seqlength,3)
  }
  else (percentage <- 0)
  unlink('temp',recursive = T)
  return(percentage)
}
