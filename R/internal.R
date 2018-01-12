#' Check if MAFFT is installed and return path
#'
#' Check if MAFFT is installed and return path to executable
#'
#' @param path Where to look for MAFFT. By default it will look on the search
#'   path.
#' @param error If \code{TRUE}, throw an error if MAFFT is not installed.
#' @param verbose If \code{TRUE}, print progress reports.
#'
#' @return TRUE/FALSE
#'
#' @keywords internal
get_mafft_path <- function(mafft.path = NULL, error = TRUE,
                           verbose = FALSE) {
  # Set default path
  if (is.null(mafft.path)) {
    path <- unname(Sys.which("mafft"))
  } else if (endsWith(mafft.path, "mafft")) {
    path <- path.expand(mafft.path)
  } else if (Sys.info()[['sysname']] %in% "Windows") {
    path <- file.path(path.expand(mafft.path),"mafft.bat")
  }  else {
    # add on executable to path if not already present
    path <- file.path(path.expand(mafft.path), "mafft")
  }


  # Check if MAFFT is installed
  is_installed <- file.exists(path) == T
  if (! is_installed && error) {
    if (is.null(mafft.path)) {
      stop(paste0("MAFFT not found in your computer's search path.",
                  "'\n Please check that MAFFT is installed and in the search path or specify the path to the MAFFT installation using the `mafft.path` option."), call. = FALSE)
    } else {
      stop(paste0("MAFFT not found in the specified path: '", path,
                  "'\n Please check your MAFFT installation."), call. = FALSE)
    }
  }
  return(path)
}

#' Check if HMMER is installed and return path
#'
#' Check if HMMER is installed and return path to an executable
#'
#' @param command Which command to get the path for. For example "hmmsearch".
#' @param path Where to look for HMMER. By default it will look on the search
#'   path.
#' @param error If \code{TRUE}, throw an error if HMMER is not installed.
#' @param verbose If \code{TRUE}, print progress reports.
#'
#' @return TRUE/FALSE
#'
#' @keywords internal
get_hmmer_path <- function(command, hmmer.path = NULL, error = TRUE,
                           verbose = FALSE) {

  # Set default path
  if (is.null(hmmer.path)) {
    path <- unname(Sys.which(command))
  } else {
    path <- file.path(path.expand(hmmer.path), command)
  }

  # Print progress
  if (verbose) {
    message(paste0("Checking if HMMER is installed in the specified path: '", path,"'"))
  }

  # Check if mafft is installed
  is_installed <- system2(path, "-h", stderr = NULL, stdout = NULL) == 0
  if (! is_installed && error) {
    if (is.null(hmmer.path)) {
      stop(paste0("HMMER not found in your computer's search path.",
                  "'\n Please check that HMMER is installed and in the search path or specify the path to the HMMER installation using the `hmm.path` option."), call. = FALSE)
    } else {
      stop(paste0("HMMER not found in the specified path: '", path,
                  "'\n Please check your HMMER installation."), call. = FALSE)
    }
  }

  return(path)
}


#' Converts FASTA to STOCKHOLM
#'
#' The HMMER binary version for windows does not recognize FASTA files to build the
#' hmm profile. This function creates a STOCKHOLM file readable by HMMER 3.0
#'
#' @param fasta.file FASTA object to be converted
#'
#' @return TRUE/FALSE
#'
#' @keywords internal
#'
fasta_to_stockholm <- function(fasta.file){
  stock.name <- gsub(fasta.file, pattern = ".fasta", replacement = ".stockholm"
  )
  seq <- seqinr::read.fasta(fasta.file)
  seq.seq <- lapply(seqinr::getSequence(seq), function (x) (paste0(x,collapse="")))
  seq.names <- seqinr::getName(seq)
  seq.final <- list()
  for (i in 1:length(seq)){
    seq.final[[i]] <- paste(seq.names[[i]],seq.seq[[i]], sep = " ")
  }
  seq.final <- unlist(seq.final)
  writeLines(c("# STOCKHOLM 1.0", seq.final,"//"), con = stock.name, sep = "\n")
}
