#' Searching for motifs using HMM searches
#'
#' This function uses MAFFT and HMMER to search for sequences with RxLR or CRN motifs using hidden markov models.
#' @param original.seq The path for the original six-frame translation FASTA file
#' @param regex.seq A list of \code{SeqFastadna} objects resulting from \code{\link{regex.search}}. The HMM profile will be constructed using these sequences
#' @param mafft.path Local path of folder containing the MAFFT binary executable file or the executable file itself. If not specified, then MAFFT must be in the program search path.
#' @param hmm.path Local path of  folder containing the HMMER binaries.  If not specified, then HMMER executables must be in the program search path.
#' @param num.threads Number of threads to be used by MAFFT
#' @param seed The seed to used with HMMER commands. Set this to get the same output each time
#' @keywords regex effector
#' @export
#' @return A list of three elements: REGEX candidate effectors, HMM candidate effectors, and HMM results table.
#' @examples
#'
#'\dontrun{
#'
#' fasta.file <- system.file("extdata", "test_infestans.fasta", package = "effectR")
#' ORF <- seqinr::read.fasta(fasta.file)
#' REGEX <- regex.search(ORF, motif="RxLR")
#' candidate.rxlr <- hmm.search(original.seq = fasta.file, regex.seq = REGEX)
#' }
#' @details
#' \code{hmm.search} uses the results from \code{\link{regex.search}} to search for motifs of interest using Hidden Markov Models after aligning the sequences with MAFFT.
#' After the multiple sequence alignment is complete, the function constructs a HMM profile using the alignment data. The HMM profile is in the original list of \code{SeqFastadna} objects to obtain the best HMM results with sequences with RxLR or CRN motifs.
#' @note
#' If MAFFT/HMMER are not the program search path, the user has to specify the path for the MAFFT and the HMMER executable binaries and specify them in the \code{mafft.path} and \code{hmm.path}
hmm.search <-  function(original.seq, regex.seq, mafft.path = NULL, num.threads = 2, hmm.path = NULL, seed = sample(1:10000, 1)){
  set.seed(seed)
  sequences <- regex.seq
  if (unique(unlist(lapply(sequences, class))) != "SeqFastadna") {
    stop("The object is not a list of sequences read by seqinr.")
  }
  if (length(regex.seq) < 4){
    stop("Not enough sequences for HMM step. At least 4 sequences are required.")
  }
  # All variable names
  # time.stamp <- gsub(format(Sys.time(), "%a_%b_%d_%X"), pattern = ":", replacement = "_")
  file.name <- c("REGEX.fasta")
  mafft.out.name <- c("MAFFT.fasta")
  hmmbuild.out <- c("hmmbuild.hmm")
  hmmsearch.out <- c("hmmsearch.txt")
  original.dir <- getwd()
  original.seq <- paste0(original.seq)

  # TMP dir
  tmp.dir <- tempdir()
  setwd(tmp.dir)
  on.exit({
    unlink(file.path(tmp.dir, "*"), recursive = TRUE)
    setwd(original.dir)
  })

  # MAFFT alignment
  cat("Starting MAFFT alignment.\n")
  cat("---\n")
  Sys.sleep(1)
  cat("Executing MAFFT\nPlease be patient\n")
  seqinr::write.fasta(sequences = seqinr::getSequence(sequences), names=seqinr::getName(sequences), file.out = file.name)
  mafft.command <- c(get_mafft_path(mafft.path),
                     "--legacygappenalty",
                     "--genafpair",
                     "--maxiterate", "1000",
                     "--thread", num.threads,
                     "--quiet", file.name)
  system2(mafft.command, stdout = mafft.out.name)
  cat("MAFFT alignment finished!")
  cat("\n")

  # HMM
  cat("Starting HMM\n")
  cat("---\n")
  cat("Creating HMM profile\n\n")
  if(file.exists(mafft.out.name) == F){
    stop("No MAFFT alignment found")
  }

  ## HMM build
  unlink(file.path(tmp.dir, hmmbuild.out))
  if (Sys.info()[['sysname']] %in% "Windows"){
  fasta_to_stockholm(fasta.file = mafft.out.name)
  stock.name <- gsub(mafft.out.name, pattern = ".fasta", replacement = ".stockholm"
    )
  hmmbuild_command <- c(get_hmmer_path("hmmbuild.exe", hmm.path),
                        "--amino",
                        "--seed", seed,
                        hmmbuild.out,
                        stock.name)
  } else {
    hmmbuild_command <- c(get_hmmer_path("hmmbuild", hmm.path),
                          "--amino",
                          "--seed", seed,
                          hmmbuild.out,
                          mafft.out.name)
  }
  system2(hmmbuild_command, stdout = F)

  ## HMM Press
  if (Sys.info()[['sysname']] %in% "Windows"){
    hmmpress_command <- c(get_hmmer_path("hmmpress.exe", hmm.path),
                        hmmbuild.out)
  } else {
    hmmpress_command <- c(get_hmmer_path("hmmpress", hmm.path),
                          hmmbuild.out)
  }
  system2(hmmpress_command)
  cat("HMM profile created.\n")

  ## HMM search
  cat("\nStarting HMM searches\n")
  if (Sys.info()[['sysname']] %in% "Windows"){
    hmmsearch_command <- c(get_hmmer_path("hmmsearch.exe", hmm.path),
                           "-T", "0",
                           "--seed", seed,
                           "--tblout", hmmsearch.out,
                           hmmbuild.out,
                           original.seq)
    } else {
    hmmsearch_command <- c(get_hmmer_path("hmmsearch", hmm.path),
                         "-T", "0",
                         "--seed", seed,
                         "--tblout", hmmsearch.out,
                         hmmbuild.out,
                         original.seq)
  }
  system2(hmmsearch_command)
  cat("\n")
  cat("hmmsearch finished!\n")

  ## Reading in hmm results
  hmm.hits <- utils::read.delim(hmmsearch.out, comment.char = "#", header = F, sep = "", blank.lines.skip = T, stringsAsFactors = F)[,1]
  if (Sys.info()[['sysname']] %in% "Windows"){
  hmm.table <- utils::read.table(hmmbuild.out, blank.lines.skip = T, skip = 14, sep = "", fill = T, stringsAsFactors = F)
  } else {
  hmm.table <- utils::read.table(hmmbuild.out, blank.lines.skip = T, skip = 16, sep = "", fill = T, stringsAsFactors = F)
  }
  total.seq <- seqinr::read.fasta(original.seq)
  hmm.seq <- total.seq[seqinr::getName(total.seq) %in% hmm.hits]
  total.seq <- list(regex.seq, hmm.seq, hmm.table)
  names(total.seq) <- c("REGEX","HMM","HMM_Table")
  cat(paste0("\nHMM search done. \n---\n\nTotal of sequences found in REGEX: ", length(total.seq[[1]]), "\n"))
  cat(paste0("Total of sequences found in HMM: ",  length(total.seq[[2]]), "\n"))
  cat(paste0("Total of redundant hits: ", sum(duplicated(unlist(lapply(total.seq[c(1,2)], function (x) seqinr::getName(x))))),"\n"))
  cat(paste0("Number of effector candidates: ",length(unique(unlist(lapply(total.seq[c(1,2)], function (x) seqinr::getName(x)))))))

  return(total.seq)
}
