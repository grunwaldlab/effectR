#' Plots the relative frequencies of each position for \code{hmmsearch} table.
#'
#' This function plots the results from \code{\link{hmm.search}} as a barplot with amino acids in the x axis and the relative frequency of each amino acid in the y axis
#' @param hmm.table The HMM profile table resulting from \code{\link{hmm.search}}
#' @keywords regex effector plot
#' @export
#' @examples
#' \dontrun{
#'
#' fasta.file <- system.file("extdata", "test_infestans.fasta", package = "effectR")
#' ORF <- seqinr::read.fasta(fasta.file)
#' REGEX <- regex.search(ORF, motif='RxLR')
#' candidate.rxlr <- hmm.search(original.seq = fasta.file, regex.seq=REGEX, num.threads = 2)
#' hmm.logo(candidate.rxlr$HMM_Table)
#'}

hmm.logo <- function (hmm.table) {
  # Plot
  hmm <- hmm.table
  colnames(hmm) <- hmm[1,]
  hmm <- hmm[-(1:5),]
  hmm[,1] <- suppressWarnings(as.numeric(hmm[,1]))
  hmm <- hmm[hmm[,1]%%1==0, ]
  hmm <- suppressWarnings(hmm[!is.na(as.numeric(hmm[,2])), ])
  rownames(hmm) <- hmm[,1]
  hmm <- hmm[,-1]
  hmm <- data.frame(sapply(hmm, function(x) as.numeric(as.character(x))), stringsAsFactors = F)
  hmm.sums <- apply(hmm,2,function (x) max(x)/x)
  hmm.sums <- apply(hmm.sums,2,function (x) x/sum(x))
  hmm.sums <- apply(hmm.sums, 2, function (x) x - mean(as.numeric(as.character(x))))
  hmm.sums[hmm.sums < 0] <- 0
  # Melt HMM
  hmm.m <- reshape2::melt(t(hmm.sums))
  colnames(hmm.m) <- c("element","position","bits")
  hmm.m$bits <- as.numeric(as.character(hmm.m$bits))
  p <- ggplot2::ggplot(hmm.m, ggplot2::aes_string(x="position", y="bits", fill="element")) + ggplot2::geom_bar(stat = "identity",position = "stack",width=1, alpha=0.5) + ggplot2::geom_text(ggplot2::aes_string(label="element", size="bits"), position='stack') +  viridis::scale_fill_viridis(discrete=TRUE) + ggplot2::theme_bw() + ggplot2::ylab("Relative Frequency (bits)") + ggplot2::guides(fill=FALSE)
  return(p)
}
