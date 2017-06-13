#' Searching for motifs using regular expressions (REGEX)
#'
#' This function uses searches a list of \code{SeqFastadna} objects for sequences with RxLR or CRN motifs.
#' @param sequence A list of \code{SeqFastadna} objects from \code{seqinr} \code{\link{read.fasta}}.
#' @param motif A character string indicating the motif to be searched (RxLR or CRN). Default is RxLR
#' @keywords regex effector
#' @export
#' @examples
#' fasta.file <- system.file("extdata", "test_infestans.fasta", package = "effectR")
#' ORF <- seqinr::read.fasta(fasta.file)
#' regex.search(ORF)
#'
regex.search <- function(sequence=sequence, motif="RxLR"){
 if (unique(unlist(lapply(sequence, class))) != "SeqFastadna") {
    stop("The object is not a list of sequences read by seqinr.")
  }
  seq <- lapply(sequence, function (x) paste(unlist(x),collapse = ""))
    regex <- list()
    for (i in 1:length(seq)){
      if (motif == "RxLR"){
        regex[[i]] <- unlist(gregexpr(seq[[i]], pattern="^\\w{12,60}r\\wlr\\w{6,10}eer", perl = T,ignore.case = T))
      } else if (motif == "CRN"){
        regex[[i]] <- unlist(gregexpr(seq[[i]], pattern="^\\w{1,90}LFLAK\\w+", perl = T,ignore.case = T))
      }
    #percentage <- percentage + 1/length(seq)*100
    }
    regex <- as.data.frame(do.call(rbind, regex))
    regex$seq <- names(seq)
    regex <- regex[!regex$V1 < 0, ]
    regex <- sequence[seqinr::getName(sequence) %in% regex$seq]
    if (length(regex) == 0){
      stop(paste0("No ",motif, " sequences found."))
    }
    return(regex)
}
