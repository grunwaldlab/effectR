library(effectR)
fasta.file <- system.file("extdata", "test_infestans.fasta", package = 'effectR')
ORF <- seqinr::read.fasta(fasta.file)
REGEX <- regex.search(sequence = ORF, motif = "RxLR")
REGEX.seq <- lapply(REGEX, function (x) paste(unlist(x),collapse = ""))

context("Generating RxLR candidates from REGEX")

test_that("effectR can read FASTA alignment correctly ", {
  expect_equal(class(ORF), "list")
  expect_equal(class(ORF[[1]]), "SeqFastadna")
  expect_equal(length(ORF), 28)
})
