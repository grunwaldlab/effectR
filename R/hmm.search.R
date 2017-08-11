#' Searching for motifs using HMM searches
#'
#' This function uses MAFFT and HMMER to search for sequences with RxLR or CRN motifs using hidden markov models.
#' @param original.seq The path for the original six-frame transaltion FASTA file
#' @param regex.seq A list of \code{SeqFastadna} objects resulting from \code{\link{regex.search}}. The HMM profile will be constructed using these sequences
#' @param mafft.path Local path of the MAFFT binary executable file
#' @param hmm.path Local path of the HMMER binaries
#' @param num.threads Number of threads to be used by MAFFT
#' @keywords regex effector
#' @export
#' @return A list of three elements: REGEX candidate effectors, HMM candidate effectors, and HMM results table.
#' @examples
#'
#'\dontrun{
#' fasta.file <- system.file("extdata", "test_infestans.fasta", package = "effectR")
#' ORF <- seqinr::read.fasta(fasta.file)
#' REGEX <- regex.search(ORF, motif="RxLR")
#' candidate.rxlr <- hmm.search(original.seq = fasta.file, regex.seq=REGEX,
#'                   mafft.path="/usr/local/bin/", hmm.path="/usr/local/bin/", num.threads = 2)
#' }
#' @details
#' \code{hmm.search} uses the results from \code{\link{regex.search}} to search for motifs of interest using hidden markov models after aligning the sequences with MAFFT.
#' After the multiple sequence alignment is complete, the function constructs a HMM profile using the alignment data. The HMM profile is in the original list of \code{SeqFastadna} objects to obtain the best HMM results with sequences with RxLR or CRN motifs.
#' @note
#' The user has to specify the path for the MAFFT and the HMMER executable binaries and specify them in the \code{mafft.path} and \code{hmm.path}

hmm.search <-  function(original.seq = "file.fasta", regex.seq = sequences, mafft.path="/usr/local/bin/", num.threads = 2, hmm.path="/usr/local/bin/"){
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

  # MAFFT alignment
  cat("Starting MAFFT alignment.\n")
  cat("---\n")
  cat(paste0("Detecting if MAFFT is installed in the specified mafft.path: "), mafft.path,"\n")
  Sys.sleep(1)
  if (file.exists(file.path(mafft.path, "mafft"))){
    cat("Executing MAFFT\nPlease be patient\n")
    seqinr::write.fasta(sequences = seqinr::getSequence(sequences), names=seqinr::getName(sequences), file.out = file.name)
    mafft.command <- paste0(mafft.path,"mafft --legacygappenalty --genafpair --maxiterate 1000 --thread ", num.threads ," --quiet ",file.name," > ",mafft.out.name)
    system(mafft.command)
    cat("MAFFT alignment finished!")
    cat("\n")}
  else {
    stop("MAFFT not found in the specified path.\n Please check your MAFFT path (go to sh and ask `which mafft`) or MAFTT installation")}
  # HMM
  cat("Starting HMM\n")
    cat("---\n")
  cat("Creating HMM profile\n\n")
    if(file.exists(mafft.out.name) == F){
      stop("No MAFFT alignment found")
    }

  ## HMM build
  required_programs <- c("hmmbuild", "hmmpress", "hmmsearch")
  if (any(! file.exists(file.path(hmm.path, required_programs)))) {
    stop(paste0("hmmbuild not found in ", hmm.path, "\nCheck your HMMER installation path\n"), .call = FALSE)
  }

  unlink(file.path(tmp.dir, hmmbuild.out))
  hmmbuild <- paste(file.path(hmm.path, "hmmbuild"), "--amino", hmmbuild.out, mafft.out.name)
  system(hmmbuild, ignore.stdout = T, ignore.stderr = F)
  system(paste0(hmm.path,"hmmpress ", hmmbuild.out), ignore.stdout = F, ignore.stderr = F)
  system(paste0("perl -pi -e 's/ {2,}/\t/g' ",hmmbuild.out))
  cat("HMM profile created.\n")

  ## HMM search
  cat("\nStarting HMM searches\n")
  system(paste0(hmm.path, "hmmsearch -T 0 --seed 1 --tblout ", hmmsearch.out," ",hmmbuild.out," ",original.seq), ignore.stdout = F, ignore.stderr = F)
  system(paste0("perl -pi -e 's/ {2,}/\t/g' ",hmmsearch.out))
  cat("\n")
  cat("hmmsearch finished!\n")

  ## Reading in hmm results
  hmm.hits <- utils::read.delim(hmmsearch.out,comment.char = "#", header = F, sep = "\t", stringsAsFactors = F)[,1]
  hmm.table <- utils::read.table(hmmbuild.out, blank.lines.skip = T, skip = 16, sep = "", fill = T, stringsAsFactors = F)
  total.seq <- seqinr::read.fasta(original.seq)
  hmm.seq <- total.seq[seqinr::getName(total.seq) %in% hmm.hits]
  total.seq <- list(regex.seq, hmm.seq, hmm.table)
  names(total.seq) <- c("REGEX","HMM","HMM_Table")
  cat(paste0("\nHMM search done. \n---\n\nTotal of sequences found in REGEX: ", length(total.seq[[1]]), "\n"))
  cat(paste0("Total of sequences found in HMM: ",  length(total.seq[[2]]), "\n"))
  cat(paste0("Total of redundant hits: ", sum(duplicated(unlist(lapply(total.seq[c(1,2)], function (x) seqinr::getName(x))))),"\n"))
  cat(paste0("Number of effector candidates: ",length(unique(unlist(lapply(total.seq[c(1,2)], function (x) seqinr::getName(x)))))))
  unlink(file.path(tmp.dir,"*"), recursive = T)
  setwd(original.dir)
  return(total.seq)
}
