library(effectR)
fasta.file <- system.file("extdata", "test_infestans.fasta", package = 'effectR')
ORF <- seqinr::read.fasta(fasta.file)
REGEX <- regex.search(sequence = ORF, motif = "RxLR")
REGEX.seq <- lapply(REGEX, function (x) paste(unlist(x),collapse = ""))
num.hits <- grep(REGEX.seq, pattern="^\\w{12,60}r\\wlr\\w{6,10}eer", perl = T,ignore.case = T)

mafft.path <- paste0(dirname(system("which mafft", intern = T)), "/")
hmm.path <- paste0(dirname(system("which hmmsearch", intern = T)), "/")
candidate.rxlr <- hmm.search(original.seq = fasta.file, regex.seq = REGEX, mafft.path = mafft.path, hmm.path = hmm.path)

summary.list <- effector.summary(candidate.rxlr)
summary.list.regex <- effector.summary(REGEX)

hmm.plot <- hmm.logo(hmm.table = candidate.rxlr$HMM_Table)


context("Testing HMM logo summary")

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
  expect_equal(length(candidate.rxlr), 3)
  expect_equal(class(candidate.rxlr), "list")
  expect_equal(names(candidate.rxlr), c("REGEX","HMM","HMM_Table"))
  expect_equal(length(candidate.rxlr$REGEX), 17)
  expect_equal(length(candidate.rxlr$HMM), 19)
  expect_equal(length(candidate.rxlr$HMM_Table), 21)
})

test_that("effector.summary returns a list with 2 objects, 19 candidates and 19 rows in summary table ", {
  expect_equal(length(summary.list), 2)
  expect_equal(class(summary.list), "list")
  expect_equal(names(summary.list), c("consensus.sequences","motif.table"))
  expect_equal(length(summary.list$consensus.sequences), 19)
  expect_equal(nrow(summary.list$motif.table), 19)
  expect_equal(unique(summary.list$motif.table$MOTIF), "Complete")
})

test_that("effector.summary can read REGEX output", {
  expect_equal(length(summary.list.regex), 1)
  expect_equal(class(summary.list.regex), "list")
  expect_equal(names(summary.list.regex), "motif.table")
  expect_equal(nrow(summary.list.regex$motif.table), 17)
  expect_equal(unique(summary.list.regex$motif.table$MOTIF), "Complete")
})

test_that("hmm.plot works ", {
  expect_equal(class(hmm.plot), c("gg","ggplot"))
})
