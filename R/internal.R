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
    path <- unname(Sys.which("mafft"))
  } else if (endsWith(mafft.path, "mafft")) {
    path <- mafft.path
  } else {
    # add on executable to path if not already present
    path <- file.path(mafft.path, "mafft")
  }


  # Check if mafft is installed
  is_installed <- system2(path, "--version", stderr = NULL) == 0
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
    path <- unname(Sys.which(command))
  } else {
    path <- file.path(hmmer.path, command)
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

