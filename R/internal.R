#' Check if MAFFT is installed and return path
#'
#' Check if MAFFT is installed and return path to executable
#'
#' @param path Where to look for MAFFT. By default it will look on the search
#'   path.
#' @param error If \code{TRUE}, throw an error if mafft is not installed.
#' @param verbose If \code{TRUE}, print progress reports.
#'
#' @return TRUE/FALSE
#'
#' @keywords internal
get_mafft_path <- function(mafft.path = NULL, error = TRUE,
                           verbose = FALSE) {
  # Set defualt path
  if (is.null(mafft.path)) {
    mafft.path <- "mafft"
  } else {
    mafft.path <- file.path(mafft.path, "mafft")
  }

  # Print progress
  if (verbose) {
    message(paste0("Checking if MAFFT is installed in the specified path: '", mafft.path,"'"))
  }

  # Check if mafft is installed
  is_installed <- system2(mafft.path, "--version", stderr = NULL) == 0
  if (! is_installed && error) {
    stop(paste0("MAFFT not found in the specified path: '", mafft.path,
                "'\n Please check your MAFFT installation."), call. = FALSE)
  }

  return(mafft.path)
}

#' Check if HMMER is installed and return path
#'
#' Check if HMMER is installed and return path to an executable
#'
#' @param command Which command to get the path for. For example "hmmsearch".
#' @param path Where to look for HMMER. By default it will look on the search
#'   path.
#' @param error If \code{TRUE}, throw an error if mafft is not installed.
#' @param verbose If \code{TRUE}, print progress reports.
#'
#' @return TRUE/FALSE
#'
#' @keywords internal
get_hmmer_path <- function(command, hmmer.path = NULL, error = TRUE,
                           verbose = FALSE) {
  # Set defualt path
  if (is.null(hmmer.path)) {
    hmmer.path <- command
  } else {
    hmmer.path <- file.path(hmmer.path, command)
  }

  # Print progress
  if (verbose) {
    message(paste0("Checking if HMMER is installed in the specified path: '", hmmer.path,"'"))
  }

  # Check if mafft is installed
  is_installed <- system2(hmmer.path, "-h", stderr = NULL, stdout = NULL) == 0
  if (! is_installed && error) {
    stop(paste0(command, " not found in the specified path: '", hmmer.path,
                "'\n Please check your HMMER installation."), call. = FALSE)
  }

  return(hmmer.path)
}

