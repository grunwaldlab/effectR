library(effectR)
fasta.file <- system.file("extdata", "test_infestans.fasta", package = 'effectR')
ORF <- seqinr::read.fasta(fasta.file)
REGEX <- regex.search(sequence = ORF, motif = "RxLR")
REGEX.seq <- lapply(REGEX, function (x) paste(unlist(x),collapse = ""))
num.hits <- grep(REGEX.seq, pattern="^\\w{12,60}r\\wlr\\w{6,10}eer", perl = T,ignore.case = T)

context("Generating RxLR candidates from HMM")

test_that("effectR can read FASTA alignment correctly ", {
  expect_equal(class(ORF), "list")
  expect_equal(class(ORF[[1]]), "SeqFastadna")
  expect_equal(length(ORF), 27)
})

test_that("regex.search returns 17 sequences with RxLR motifs ", {
  expect_equal(length(REGEX), 17)
  expect_equal(length(num.hits), 17)
  expect_equal(num.hits, c(1:17))
})

test_that("candidate.rxlr returns a list with 3 objects, 17 REGEX, 19 HMM and 19 rows in HMM table ", {
  skip_on_cran()
  candidate.rxlr <- hmm.search(original.seq = fasta.file, regex.seq = REGEX, seed = 1)
  expect_equal(length(candidate.rxlr), 3)
  expect_equal(class(candidate.rxlr), "list")
  expect_equal(names(candidate.rxlr), c("REGEX","HMM","HMM_Table"))
  expect_equal(length(candidate.rxlr$REGEX), 17)
  expect_equal(length(candidate.rxlr$HMM), 19)
  expect_equal(length(candidate.rxlr$HMM_Table), 21)
})

test_that("Invalid dependency paths cause understandable errors ", {
  expect_error(hmm.search(original.seq = fasta.file, regex.seq = REGEX, seed = 1,
                          mafft.path = "This is not right..."),
               "MAFFT not found in the specified path")
  expect_error(hmm.search(original.seq = fasta.file, regex.seq = REGEX, seed = 1,
                          hmm.path = "/something/very/misguided"),
               "HMMER not found in the specified path")
})

test_that("MAFFT path option accepts both directory and file ", {
  skip_on_cran()
  mafft_path <- system2(c("which", "mafft"), stdout = TRUE)
  expect_equal({
    set.seed(1)
    hmm.search(original.seq = fasta.file, regex.seq = REGEX, seed = 1, mafft.path = mafft_path)
  }, {
    set.seed(1)
    hmm.search(original.seq = fasta.file, regex.seq = REGEX, seed = 1, mafft.path = dirname(mafft_path))
  })
})
