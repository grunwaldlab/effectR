#' Returns non-redundant sequences from hmm.search or regex.search and generates a motif table
#'
#' This function summarize the results from \code{\link{regex.search}} or \code{\link{hmm.search}}.
#' @param hmm.result A list of \code{SeqFastadna} objects obtained from \code{\link{regex.search}} or \code{\link{hmm.search}}
#' @param motif A character string indicating the motif to be searched (RxLR or CRN). Default is RxLR
#' @keywords regex effector
#' @export
#' @return  A list of two objects: Summary motif table and non-redundant sequences (only with results of \code{\link{hmm.search}})
#' @examples
#' fasta.file <- system.file("extdata", "test_infestans.fasta", package = "effectR")
#' ORF <- seqinr::read.fasta(fasta.file)
#' REGEX <- regex.search(ORF, motif='RxLR')
#' candidate.rxlr <- hmm.search(original.seq = fasta.file, regex.seq=REGEX, mafft.path="/usr/local/bin/", hmm.path="/usr/local/bin/")
#' effector.summary(candidate.rxlr)
#'
effector.summary <- function (hmm.result=candidate.rxlr, motif="RxLR"){
  motif.out <- list()
  if (length(hmm.result) == 3) {
    summ.dataf <- c(hmm.result[[1]],hmm.result[[2]])
    consensus.seq <- summ.dataf[!duplicated(seqinr::getName(summ.dataf))]
    motif.out$consensus.sequences <- consensus.seq
  } else {
    consensus.seq <- hmm.result
  }
  sequences <- lapply(consensus.seq, function (x) paste(unlist(x),collapse = ""))
  if (motif == "RxLR"){
    rxlr.num <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="r\\wlr", perl = T,ignore.case = T))), function (x) length(x)))
    rxlr.motif <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="r\\wlr", perl = T,ignore.case = T))), function (x) paste(x,collapse = ",")))
    eer.num <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="[ED][ED][KR]", perl = T,ignore.case = T))), function (x) length(x)))
    eer.motif <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="[ED][ED][KR]", perl = T,ignore.case = T))), function (x) paste(x,collapse = ",")))
    motifs <- data.frame(seqinr::getName(consensus.seq),rxlr.num,rxlr.motif,eer.num,eer.motif, stringsAsFactors = F)
  } else if (motif == "CRN") {
    lflak.num <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="lflak", perl = T,ignore.case = T))), function (x) length(x)))
    lflak.motif <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="lflak", perl = T,ignore.case = T))), function (x) paste(x,collapse = ",")))
    hvlv.num <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="hvlv", perl = T,ignore.case = T))), function (x) length(x)))
    hvlv.motif <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="hvlv", perl = T,ignore.case = T))), function (x) paste(x,collapse = ",")))
    motifs <- data.frame(seqinr::getName(consensus.seq),lflak.num,lflak.motif,hvlv.num,hvlv.motif, stringsAsFactors = F)
  }
  motifs[motifs == -1] <- NA
  motifs <- data.frame(apply(motifs, 2, as.character), stringsAsFactors = F)
  motifs[,2][is.na(motifs[,3])] <- 0
  motifs[,4][is.na(motifs[,5])] <- 0
  motifs <- motifs[order(motifs[,2],motifs[,4],decreasing = T),]
  motifs$summary <- "No MOTIFS"
  motifs$summary[motifs[,2] > 0 & motifs[,4] > 0] <- "Complete"
  motifs$summary[motifs[,2] != 0 & motifs[,4] == 0] <- if (motif == "RxLR"){                                      c("Only RxLR motif")} else if (motif == "CRN"){("Only LFLAK motif")}
  motifs$summary[motifs[,2] == 0 & motifs[,4] != 0] <- if (motif == "RxLR"){                                      c("Only EER motif")} else if (motif == "CRN"){("Only HVLV motif")}
  if (motif == "RxLR"){
    colnames(motifs) <- c("Sequence ID","RxLR number","RxLR position","EER number","EER position","MOTIF")
  } else if (motif == "CRN") {
    colnames(motifs) <- c("Sequence ID","LFLAK number","LFLAK position","HVLV number","HVLV position","MOTIF")
  }
  motif.out$motif.table <- motifs
  return(motif.out)
}

