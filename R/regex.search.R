#' Searching for motifs using regular expressions (REGEX)
#'
#' This function uses searches a list of \code{SeqFastadna} objects for sequences with RxLR or CRN motifs.
#' @param sequence A list of \code{SeqFastadna} objects from \code{seqinr} \code{\link{read.fasta}}. The \code{SeqFastadna} object must be comprised by amino acid sequences, not DNA sequences
#' @param motif A character string indicating the motif to be searched. Motifs for two cytoplasmic effectors are added to the function: \code{RxLR} or \code{CRN} effectors. Each of these motifs are associated with a by-default REGEX (\code{reg.pat}). A third option, \code{custom}, allows for the search of custom motifs. The \code{custom} option requires the specification of the motif REGEX pattern in the \code{reg.pat} option, in a \code{\link{regex}} format. Default is RxLR
#' @param reg.pat A character string indicating the REGEX pattern for the \code{custom} motif. The specification of the REGEX pattern in must be in \code{\link{regex}} format. Required for \code{custom} option of \code{motif}
#' @keywords regex effector
#' @export
#' @examples
#' fasta.file <- system.file("extdata", "test_infestans.fasta", package = "effectR")
#' ORF <- seqinr::read.fasta(fasta.file)
#' rxlr.cand <- regex.search(ORF)
#' custom.cand <- regex.search(ORF, motif = "custom", reg.pat ="^\\w{12,60}r\\wlr\\w{6,10}eer")
#'
regex.search <- function(sequence, motif = "RxLR", reg.pat = NULL){
 if (unique(unlist(lapply(sequence, class))) != "SeqFastadna") {
    stop("The object is not a list of sequences read by seqinr.")
  }
  seq <- lapply(sequence, function (x) paste(unlist(x),collapse = ""))
    regex <- list()
    if (motif %in% c("RxLR","CRN") & !is.null(reg.pat)){
      message(paste0("Custom REGEX patterns are not supported with the 'CRN' or 'RxLR' motif options.\n The package will use the default REGEX patterns used to search for ", motif, " motifs."))
      Sys.sleep(2)
    }
    for (i in 1:length(seq)){
      if (motif == "RxLR"){
        reg.pat <- "^\\w{10,40}\\w{1,96}R\\wLR\\w{1,40}eer"
      } else if (motif == "CRN"){
        reg.pat <- "^\\w{1,90}LFLAK\\w+"
      } else if (motif == "custom"){
        if (is.null(reg.pat)){
          stop("No custom REGEX pattern found.\n The 'custom' option requires a mandatory REGEX pattern")
        } else {
        reg.pat <- reg.pat
        }
      }
      regex[[i]] <- unlist(gregexpr(seq[[i]], pattern = reg.pat, perl = T ,ignore.case = T))
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
