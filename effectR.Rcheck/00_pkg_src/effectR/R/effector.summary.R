#' Returns non-redundant sequences from hmm.search or regex.search and generates a motif table
#'
#' This function summarize the results from \code{\link{regex.search}} or \code{\link{hmm.search}}.
#' @param hmm.result A list of \code{SeqFastadna} objects obtained from \code{\link{regex.search}} or \code{\link{hmm.search}}
#' @param motif A character string indicating the motif of interest. Motifs for two cytoplasmic effectors are added to the function: \code{RxLR} or \code{CRN} effectors. Each of these motifs are associated with a by-default REGEX (\code{reg.pat}). A third option, \code{custom}, allows for the search of custom motifs. The \code{custom} option requires the specification of the motif REGEX pattern in the \code{reg.pat} option, in a \code{\link{regex}} format. Default is RxLR
#' @param reg.pat A character string indicating the REGEX pattern for the \code{custom} motif. The specification of the REGEX pattern in must be in \code{\link{regex}} format. Required for \code{custom} option of \code{motif}
#' @keywords regex effector
#' @export
#' @return  A list of two objects: Summary motif table and non-redundant sequences (only with results of \code{\link{hmm.search}})
#' @examples
#' \dontrun{
#'
#' fasta.file <- system.file("extdata", "test_infestans.fasta", package = "effectR")
#' ORF <- seqinr::read.fasta(fasta.file)
#' REGEX <- regex.search(ORF, motif='RxLR')
#' candidate.rxlr <- hmm.search(original.seq = fasta.file, regex.seq=REGEX, num.threads = 2)
#' effector.summary(candidate.rxlr)
#' # Custom motifs
#' reg.pat <- "^\\w{50,60}[w,v]"
#' REGEX <- regex.search(sequence = ORF, motif = "custom", reg.pat = reg.pat)
#' candidate.custom <- hmm.search(original.seq = fasta.file, regex.seq = REGEX)
#' effector.summary(candidate.custom, motif = "custom", reg.pat = reg.pat)
#'
#'}
effector.summary <- function (hmm.result, motif="RxLR", reg.pat=NULL){
  motif.out <- list()
  if (length(hmm.result) == 4){
    hmm.result <- hmm.result[names(hmm.result) %in% c("REGEX","HMM","HMM_Table")]
  }
  if (length(hmm.result) == 3) {
    summ.dataf <- c(hmm.result[[1]],hmm.result[[2]])
    consensus.seq <- summ.dataf[!duplicated(seqinr::getName(summ.dataf))]
    motif.out$consensus.sequences <- consensus.seq
  } else {
    consensus.seq <- hmm.result
  }
  sequences <- lapply(consensus.seq, function (x) paste(unlist(x),collapse = ""))
  if (motif == "RxLR"){
    rxlr.num <- as.numeric(gsub(pattern = " ", replacement = "",unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="r\\wlr", perl = T,ignore.case = T))), function (x) length(x)))))
    rxlr.motif <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="r\\wlr", perl = T,ignore.case = T))), function (x) paste(x,collapse = ",")))
    eer.num <- as.numeric(gsub(pattern = " ", replacement = "", unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="eer", perl = T,ignore.case = T))), function (x) length(x)))))
    eer.motif <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="eer", perl = T,ignore.case = T))), function (x) paste(x,collapse = ",")))
    motifs <- data.frame(seqinr::getName(consensus.seq),rxlr.num,rxlr.motif,eer.num,eer.motif, stringsAsFactors = F)
  } else if (motif == "CRN") {
    lflak.num <- as.numeric(gsub(pattern = " ", replacement = "", unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="lflak", perl = T,ignore.case = T))), function (x) length(x)))))
    lflak.motif <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="lflak", perl = T,ignore.case = T))), function (x) paste(x,collapse = ",")))
    hvlv.num <- as.numeric(gsub(pattern = " ", replacement = "",unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="hvlv", perl = T,ignore.case = T))), function (x) length(x)))))
    hvlv.motif <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="hvlv", perl = T,ignore.case = T))), function (x) paste(x,collapse = ",")))
    seq.length <-
    motifs <- data.frame(seqinr::getName(consensus.seq),lflak.num,lflak.motif,hvlv.num,hvlv.motif, stringsAsFactors = F)
  } else if (motif == "custom"){
    if (is.null(reg.pat)){
      stop("No custom REGEX pattern found.\n The 'custom' option requires a mandatory REGEX pattern")
    } else {
      custom.num <- as.numeric(gsub(pattern = " ", replacement = "", unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern=reg.pat, perl = T,ignore.case = T))), function (x) length(x)))))
      custom.motif <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern=reg.pat, perl = T,ignore.case = T))), function (x) paste(x,collapse = ",")))
      motifs <- data.frame(seqinr::getName(consensus.seq),custom.num,custom.motif, stringsAsFactors = F)
    }
  }
  motifs[motifs == -1] <- NA
  if (motif %in% c("RxLR","CRN")){
    motifs[,2][is.na(motifs[,3])] <- 0
    motifs[,4][is.na(motifs[,5])] <- 0
    motifs <- motifs[order(motifs[,2],motifs[,4],decreasing = T),]
    motifs$summary <- "No MOTIFS"
    motifs$summary[motifs[,2] > 0 & motifs[,4] > 0] <- "Complete"
    motifs$summary[motifs[,2] != 0 & motifs[,4] == 0] <- if (motif == "RxLR"){c("Only RxLR motif")} else if (motif == "CRN"){("Only LFLAK motif")}
    motifs$summary[motifs[,2] == 0 & motifs[,4] != 0] <- if (motif == "RxLR"){c("Only EER motif")} else if (motif == "CRN"){("Only HVLV motif")}
    if (motif == "RxLR"){
      colnames(motifs) <- c("Sequence ID","RxLR number","RxLR position","EER number","EER position","MOTIF")
    } else if (motif == "CRN") {
      colnames(motifs) <- c("Sequence ID","LFLAK number","LFLAK position","HVLV number","HVLV position","MOTIF")
    }
  } else if (motif == "custom"){
    motifs[,2][is.na(motifs[,3])] <- 0
    motifs <- motifs[order(motifs[,2],decreasing = T),]
    motifs$summary <- "No MOTIFS"
    motifs$summary[motifs[,2] > 0] <- "Custom motif"
    colnames(motifs) <- c("Sequence ID","Motif number","Motif position","MOTIF")
  }
  motifs <- data.frame(motifs, length = unlist(lapply(motifs[,1], function (x) length(unlist(consensus.seq[seqinr::getName(consensus.seq) %in% x])))))
  rownames(motifs) <- NULL
  motif.out$motif.table <- motifs
  return(motif.out)
}


