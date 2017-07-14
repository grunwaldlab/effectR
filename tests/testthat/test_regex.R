library(effectR)
fasta.file <- system.file("extdata", "test_infestans.fasta", package = 'effectR')
ORF <- seqinr::read.fasta(fasta.file)
REGEX <- regex.search(sequence = ORF, motif = "RxLR")
REGEX.seq <- lapply(REGEX, function (x) paste(unlist(x),collapse = ""))
num.hits <- grep(REGEX.seq, pattern="^\\w{12,60}r\\wlr\\w{6,10}eer", perl = T,ignore.case = T)


context("Generating RxLR candidates from REGEX")

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
