library(effectR)
fasta.file <- system.file("extdata", "test_infestans.fasta", package = 'effectR')
ORF <- seqinr::read.fasta(fasta.file)
REGEX <- regex.search(sequence = ORF, motif = "RxLR")
REGEX.seq <- lapply(REGEX, function (x) paste(unlist(x),collapse = ""))
num.hits <- grep(REGEX.seq, pattern="^\\w{12,60}r\\wlr\\w{6,10}eer", perl = T,ignore.case = T)

get_mafft_path <- function(mafft.path = NULL, error = TRUE,
                           verbose = FALSE) {
  # Set default path
  if (is.null(mafft.path)) {
    path <- unname(Sys.which("mafft"))
  } else if (endsWith(mafft.path, "mafft")) {
    path <- mafft.path
  } else if (Sys.info()[['sysname']] %in% "Windows") {
    path <- file.path(mafft.path,"mafft.bat")
  }  else {
    # add on executable to path if not already present
    path <- file.path(mafft.path, "mafft")
  }


  # Check if mafft is installed
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

get_hmmer_path <- function(command="hmmsearch", hmmer.path = NULL, error = TRUE,
                           verbose = FALSE) {

  # Set default path
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


test.mafft <- try(get_mafft_path(), silent = T)
test.hmmer <- try(get_hmmer_path(), silent = T)

context("Generating RxLR candidates from HMM")

test_that("effectR can read FASTA alignment correctly ", {
  expect_equal(class(ORF), "list")
  expect_equal(class(ORF[[1]]), "SeqFastadna")
  expect_equal(length(ORF), 28)
})



if (class(test.mafft) != "try-error" || class(test.hmmer) != "try-error"){
  test_that("candidate.rxlr returns a list with 3 objects, 17 REGEX, 19 HMM and 19 rows in HMM table ", {
  skip_on_cran()
  candidate.rxlr <- hmm.search(original.seq = fasta.file, regex.seq = REGEX, seed = 1)
  expect_equal(length(candidate.rxlr), 3)
  expect_equal(class(candidate.rxlr), "list")
  expect_equal(names(candidate.rxlr), c("REGEX","HMM","HMM_Table"))
  expect_equal(length(candidate.rxlr$REGEX), 15)
  expect_equal(length(candidate.rxlr$HMM), 17)
})

test_that("Invalid dependency paths cause understandable errors ", {
  skip_on_cran()
  expect_error(hmm.search(original.seq = fasta.file, regex.seq = REGEX, seed = 1,
                          mafft.path = "This is not right..."),
               "MAFFT not found in the specified path")
  expect_error(hmm.search(original.seq = fasta.file, regex.seq = REGEX, seed = 1,
                          hmm.path = "/something/very/misguided"),
               "HMMER not found in the specified path")
})
}
