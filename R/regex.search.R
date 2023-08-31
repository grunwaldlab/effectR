#' Searching for motifs using regular expressions (REGEX)
#'
#' This function uses searches a list of \code{SeqFastadna} objects for sequences with RxLR or CRN motifs.
#' @param sequence A list of \code{SeqFastadna} objects from \code{seqinr} \code{\link{read.fasta}}. The \code{SeqFastadna} object must be comprised by amino acid sequences, not DNA sequences
#' @param motif A character string indicating the motif to be searched. Motifs for two cytoplasmic effectors are added to the function: \code{RxLR} or \code{CRN} effectors. Each of these motifs are associated with a by-default REGEX (\code{reg.pat}). These motifs are adapted from Haas et al. (2009). Additional REGEX are available using "Win2007" or "Whisson2007", as described in Haas et al. 2009 and presented in Table S-A2.
#'
#' A third option, \code{custom}, allows for the search of custom motifs. The \code{custom} option requires the specification of the motif REGEX pattern in the \code{reg.pat} option, in a \code{\link{regex}} format.
#'
#' Default \code{motif} is RxLR
#'
#' @param reg.pat A character string indicating the REGEX pattern for the \code{custom} motif. The specification of the REGEX pattern in must be in \code{\link{regex}} format. Required for \code{custom} option of \code{motif}
#' @keywords regex effector
#' @references Haas, B.J., Kamoun, S., Zody, M.C., Jiang, R.H., Handsaker, R.E., Cano, L.M., Grabherr, M.,
#'  Kodira, C.D., Raffaele, S., Torto-Alalibo, T. and Bozkurt, T.O., 2009.
#'  Genome sequence and analysis of the Irish potato famine pathogen Phytophthora infestans.
#'  Nature, 461(7262), p.393.
#'
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
        reg.pat <- "^\\w{10,40}\\w{1,96}R\\wLR\\w{1,40}EER"
      } else if (motif == "CRN"){
        reg.pat <- "^\\w{1,90}LFLAK\\w+"
      } else if (motif == "win2007" | motif == "Win2007") {
        # Haas et al. 2009: "a signal peptide within residues 1-30
        # followed by an RXLR motif within residues 30-60"
        reg.pat <- "^\\w{1,30}\\w{1,26}R\\wLR"
      } else if (motif == "whisson2007" | motif == "Whisson2007"){
        # Haas et al. 2009: "allow for a signal peptide between residues 10-40,
        # followed by the RXLR motif within the next 100 residues,
        # followed by the EER motif, allowing D and K replacements to E and R"
        reg.pat <- "^\\w{10,40}\\w{1,96}R\\wLR\\w{1,40}[ED][ED][KR]"
      } else if (motif == "whisson2007_rxlr" | motif == "Whisson2007_rxlr" |
                 motif == "whisson2007_RXLR" | motif == "Whisson2007_RXLR") {
        # Same as Whisson but ONLY RXLR, for flexibility
        reg.pat <- "^\\w{10,40}\\w{1,96}R\\wLR"
      } else if (motif == "whisson2007_eer" | motif == "Whisson2007_eer" |
                 motif == "whisson2007_EER" | motif == "Whisson2007_EER") {
        # Same as Whisson but only the EER motif, to add flexibility
        reg.pat <- "^\\w{10,40}\\w{1,100}\\w{1,40}[ED][ED][KR]"
      } else if (motif == "ai2020" | motif == "Ai2020") {
        # Ai et al. 2020 modified Whisson to find Pythium RXLRs since
        # "... some functionally verified RXLRs contain degenerate dEER motifs (Supplementary Table S2)"
        # > "The original regex model (SP.{1,96}R.LR.{1,40}[ED][ED][KR])
        #    was modified (SP.{1,40}R.LR.{1,40}([ED].[ED][KR]|[ED][ED].{0,3}[KR])"
        reg.pat <- "^\\w{10,40}\\w{1,40}R\\wLR\\w{1,40}(([ED]\\w[ED][KR])|([ED][ED]\\w{0,3}[KR]))"
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
