# effectR: An R package to call cytoplasmic effectors

[![CRAN version](http://www.r-pkg.org/badges/version/effectR)](https://cran.r-project.org/package=effectR)

# What is *effectR*

The `effectR` package is an R package designed to call oomycete RxLR and CRN effectors by searching for the motifs of interest using regular expression searches and hidden markov models (HMM). 

# What's new in version 1.0.2?

- New test dataset with *P. infestans* reference RxLR effectors 

# Overview

The `effectR` packages searches for the motifs of interest (RxLR-EER motif for RxLR effectors and LFLAK motif for CRN effectors) using a regular expression search (`REGEX`). 
These motifs used by the REGEX `effectR` search have been reported in the literature ([Haas et al., 2009](https://www.nature.com/nature/journal/v461/n7262/full/nature08358.html), [Stam et al., 2013](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0059517)).

The `effectR` package aligns the REGEX search results using [`MAFFT`](http://mafft.cbrc.jp/alignment/software/), and builds a HMM profile based on the multiple sequence alignment result using the `hmmbuild` program from [`HMMER`](http://hmmer.org/). The HMM profile is used to search across ORF of the genome of interest using the `hmmsearch` binary from `HMMER`. 
The search step will retain sequences with significant hits to the profile of interest. 
`effectR` also combines the redundant sequences found in both REGEX and HMM searches into a single dataset that can be easily exported. 
In addition, `effectR` reads and returns the HMM profile to the user and allows for the creation of a [MOTIF logo-like plot](https://en.wikipedia.org/wiki/Sequence_logo) using `ggplot2`.

## Requirements 

- R packages:
 - [`seqinr`](https://CRAN.R-project.org/package=seqinr/seqinr.pdf)
 - [`ggplot2`](http://ggplot2.org/)
 
- External software 
  - [`MAFFT`](mafft.cbrc.jp/alignment/software/)
  - [`HMMER`](http://hmmer.org/)
 
## Installation

### From CRAN

The latest version of `effectR` can be installed from CRAN. To install, make sure R is at least version 3.4.0. In the R console type

```
install.packages("effectR")
```

### From GitHub

To install `effectR` via GitHub, make sure that the `devtools` package is installed (use `install.packages("devtools")`). After installing devtools, in the R console type:

```
devtools::install_github(repo = "grunwaldlab/effectR", build_vignettes = TRUE)
library("effectR")
```

### External software (REQUIRED)

The `effectR` package uses `MAFFT` and `HMMER3` to perform the hidden markov model seach across the results from the REGEX step. These two packages should be installed before running any of the `effectR` functions.

### Downloading and installing MAFFT

MAFFT is a multiple sequence alignment program that uses Fourier-transform algorithms to align multiple sequences. We recommend downloading and installing MAFFT by following the instructions and steps in the [MAFFT installation](http://mafft.cbrc.jp/alignment/software/) web site. 

#### Linux/OS X Users

Make sure that you remember the directory in which `MAFFT` is installed, of if the installation is sucessful, make sure to obtain the path via bash/tsh/console:

```bash
which mafft
```

```
/usr/local/bin/mafft
```
For more information about MAFFT go to the MAFFT website: http://mafft.cbrc.jp/

#### Windows Users

MAFFT comes in two main distributions for windows:

- [A version that requires a UNIX-like interface called Cygwin](http://mafft.cbrc.jp/alignment/software/windows_cygwin.html)
- [An "all-in-one" version"](http://mafft.cbrc.jp/alignment/software/windows_without_cygwin.html)

Please, download and install the **all-in-one** version.
We recommend that you download and save MAFFT in your Desktop, as it will make yyour path easily accesible.

### Downloading and installing HMMER

HMMER is used for searching sequence databases for sequence homologs. It uses hidden Markov models (profile HMMs) to search for sequences with hits to similar patterns than the profile. We use three main HMMER tools: 

  - `hmmbuild` to create the HMM database from the sequences ontained in the REGEX step of `effectR`
  - `hmmpress` converts the HMM database into a format usable by other `HMMER` programs
  - `hmmsearch` to excecute the HMM search in our sequence queries basde on the HMM profile

The `effectR` package requires all of these tools. A correct `HMMER` installation will install all three programs.

#### Linux/OS X users 

We recommend downloading and installing HMMER by following the instructions and steps in the [HMMER installation](eddylab.org/software/hmmer3/3.1b2/Userguide.pdf) web site. Make sure that you remember the directory in which `HMMER` is installed, of if the installation is sucessful, make sure to obtain the path via bash/tsh/console:

```bash
which hmmbuild
which hmmpress
which hmmsearch
```

```
/usr/local/bin/hmmbuild
/usr/local/bin/hmmpress
/usr/local/bin/hmmsearch
```

For more information about HMMER go to the HMMER website: http://hmmer.org/

#### Windows users

To use the `effectR` package in Windows, the user **must** download the Windows binaries of [HMMER](http://hmmer.org/binaries/hmmer3.0_windows.zip). `effectR` will not work with any other version of HMMER.  
## Data input

The `effectR` package is designed to work with amino acid sequences in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) representing the six-frame translation of every open reading frame (ORF) of an oomycete genome. 
Using the six-frame translation of all ORF's in a genome is recommended in order to obtain as many effectors as possible from a proteome. 
To obtain the ORF for a genome, we recommend the use of EMBOSS' [`getorf`](http://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html).

`effectR` uses a list of sequences of the class `SeqFastadna` in order to perform the effector searches. 
The function `read.fasta` from the `seqinr` package reads the FASTA amino acid file into R, creating a list of `SeqFastadna` objects that represent each of the translated ORF's from the original FASTA file. 

## REGEX search

To perform the effector search, `effectR` searches for the motifs of interest found in RxLR and CRN motifs. 
We have created the function `regex.search` to perform the seach of the motif of interest. 
The function `regex.search` requires the list of `SeqFastadna` objects and the gene family of interest. 

## Multiple sequence alignment and HMMER search

To perform the HMM search and obtain all possible effector candidates from a proteome, `effectR` uses the `REGEX` results as a template to create a HMM profile and perform a search across the proteome of interest.
We have created the `hmm.search` function in order to perfomr this search.
The `hmm.search` function requires a local installation of `MAFFT` and `HMMER` in order to perform the searches. The *absolute paths* of the binaries must be specified in the `mafft.path` and `hmmer.path` options of the `hmm.search` function.
In addition, the `hmm.function` requires the path of the original FASTA file containing the translated ORF's in the `original.seq` parameter of the function. `hmm.search` will use this file as a query in the `hmmsearch` software from HMMER, and search for all sequences with hits against the HMM profile created with the REGEX results.

A default `hmm.search` object returns a list of 3 elements: 

1. The REGEX sequences used to build the HMM profile in a `SeqFastadna` class
2. The sequences from the original translated ORF files with hits to the HMM profile in a `SeqFastadna` class
3. The HMM profile table created by HMMER's `hmmbuild` as a data frame

*NEW FEATURES*:
- `hmm.search` can use a user-defined alignment file (i.e. A multiple sequence alignment performed in MUSCLE, ClustalW, etc.) and omit the alignment step
- `hmm.search` allows the user to save the multiple sequence alignment created by MAFFT within the function

More information on these new features is available in the package help (`?hmm.search`) or in the *effectR* vignette

## Obtaining non-redundant effectors and motif summaries

The user can extract all of the non-redundant sequences and a summary table with the information about the motifs using the `effector.summary` function. 
This function uses the results from either `hmm.seach` or `regex.search` functions to generate a table that includes the name of the candidate effector sequence, the number of motifs of interest (RxLR-EER or LFLAK-HVLV) per sequence and its location within the sequence. 
In addition, when the `effector.summary` function is used in an object that contains the results of `hmm.search`, the user will obtain a list of the non-reduntant sequences. If the user provides the results from `regex.search`, the function will return the motif summary table.

The motif table has a column called **MOTIF**.
This column summarizes the candidate ORF into one of 4 categories:

- Complete: The candidate ORF has both motifs of interest (RxLR + EER or LFLAK + HVLV)
- Only RxLR/Only LFLAK: The candidate ORF only has the translocation domain
- Only EER/HVLV: The candidate ORF only has the second motif of interest
- No motifs: The sequence has a hit with the HMM profile but does not have any motif of interest

### Non-redundant sequences

To export the non-redundant effector candidates that resulted from the `hmm.search` or `regex.search` functions, we use the `write.fasta` function of the `seqinr` package. 
We recomend the users to read the documentation of the [`seqinr`](https://cran.r-project.org/package=seqinr) package
Since the objects that result from the `hmm.search` or `regex.search` function are of the `SeqFastadna` class, we can use any of the function of the `seqinr` package that use this class as well. 

## Visualizing the HMM profile using a sequence logo-like plot

To determine if the HMM profile includes the motifs of interest, we have created the function `hmm.logo`.
The function `hmm.logo` reads the HMM profile (obtained from the `hmm.search` step) and uses `ggplot2` to create a bar-plot.
The bar-plot will illustrate the bits (aminoacid score) of each amino acid used to construct the HMM profile according to its consensus position in the HMM profile. 
To learn more about sequence logo plots visit this [wikipedia article](https://en.wikipedia.org/wiki/Sequence_logo).


## Adding custom motifs to effectR

The effectR package has the capability to use custom regular expressions to predict other families of genes of interest other than RxLR/CRN effector proteins. This example uses the PAAR motif (PAAR) identified in proteins associated with the terminal spike in T6SS of some bacterial species:

```
# Loading the effectR package
library(“effectR”)

# Using the read.fasta function of the seqinr package to import the V. cholerae FASTA proteome file
fasta.file <- " V_cholerae_ATCC_39315.AA.fasta.gz”
ORF <- seqinr::read.fasta(fasta.file)

# Step 1 prediction: Since the PAAR motif can occur anywhere in the sequence of the protein, the REGEX will be very simple and will only contain the PAAR motif.
REGEX <- regex.search(ORF, motif = "custom", reg.pat = "PAAR")
```

Step 1 resulted in a total of 19 predicted proteins with the PAAR motif. We can expand the number of candidate PAAR proteins using the HMM step:

```
# Expanding the search of RxLR effectors using HMM searches (step 2). All candidate effectors predicted by step 2 will be saved in the candidate.rxlr object

candidate.paar <- hmm.search(original.seq = fasta.file, regex.seq = REGEX)

```

Step 2 resulted in one additional candidate protein with a plausible PAAR motif. We can summarize all the information from effectR using the `effector.summary()` function. It will return a table with the candidate proteins, the number of PAAR motifs within each protein, and the position of said motif:

```
# Summarizing the results of effectR.

effector.summary(candidate.paar)

                                  Sequence ID Motif number Motif position        MOTIF
tr|Q9KN60|Q9KN60_VIBCH tr|Q9KN60|Q9KN60_VIBCH            2          35,71 Custom motif
sp|Q9KR02|RUVB_VIBCH     sp|Q9KR02|RUVB_VIBCH            1            145 Custom motif
sp|Q9KPV0|GLND_VIBCH     sp|Q9KPV0|GLND_VIBCH            1            398 Custom motif
sp|Q9KSQ2|HUTG_VIBCH     sp|Q9KSQ2|HUTG_VIBCH            1            269 Custom motif
sp|Q9KPU5|NUSB_VIBCH     sp|Q9KPU5|NUSB_VIBCH            1              7 Custom motif
tr|Q9KS85|Q9KS85_VIBCH tr|Q9KS85|Q9KS85_VIBCH            1            171 Custom motif
tr|Q9KUC8|Q9KUC8_VIBCH tr|Q9KUC8|Q9KUC8_VIBCH            1            232 Custom motif
tr|Q9KND6|Q9KND6_VIBCH tr|Q9KND6|Q9KND6_VIBCH            1            153 Custom motif
tr|Q9KN94|Q9KN94_VIBCH tr|Q9KN94|Q9KN94_VIBCH            1            236 Custom motif
tr|Q9KMP0|Q9KMP0_VIBCH tr|Q9KMP0|Q9KMP0_VIBCH            1             35 Custom motif
tr|Q9KPU1|Q9KPU1_VIBCH tr|Q9KPU1|Q9KPU1_VIBCH            1            179 Custom motif
tr|Q9KSK0|Q9KSK0_VIBCH tr|Q9KSK0|Q9KSK0_VIBCH            1            501 Custom motif
tr|Q9KLU5|Q9KLU5_VIBCH tr|Q9KLU5|Q9KLU5_VIBCH            1            481 Custom motif
tr|Q9KKR8|Q9KKR8_VIBCH tr|Q9KKR8|Q9KKR8_VIBCH            1            343 Custom motif
tr|Q9KPP4|Q9KPP4_VIBCH tr|Q9KPP4|Q9KPP4_VIBCH            1            811 Custom motif
tr|Q9KUF6|Q9KUF6_VIBCH tr|Q9KUF6|Q9KUF6_VIBCH            1            290 Custom motif
tr|Q9KVN5|Q9KVN5_VIBCH tr|Q9KVN5|Q9KVN5_VIBCH            1             74 Custom motif
tr|Q9KQJ5|Q9KQJ5_VIBCH tr|Q9KQJ5|Q9KQJ5_VIBCH            1            189 Custom motif
tr|Q9KNT2|Q9KNT2_VIBCH tr|Q9KNT2|Q9KNT2_VIBCH            1            146 Custom motif
tr|Q9KUM6|Q9KUM6_VIBCH tr|Q9KUM6|Q9KUM6_VIBCH            0           <NA>    No MOTIFS
```

The results illustrate that 19 out of the 20 candidate proteins have a predicted PAAR domain within its sequence, and only protein `tr|Q9KN60|Q9KN60_VIBCH` has more than 1 PAAR motif. 

Any user can add the motif of interest into the effectR package by adding a simple line of code within the `regex.search` function. We will illustrate this feature by adding the PAAR motif search as part of the `regex.search` function:

```
# The regex.search function
regex.search <- function(sequence, motif = "RxLR", reg.pat = NULL){
 if (unique(unlist(lapply(sequence, class))) != "SeqFastadna") {
    stop("The object is not a list of sequences read by seqinr.")
  }
  seq <- lapply(sequence, function (x) paste(unlist(x),collapse = ""))
    regex <- list()
    if (motif %in% c("RxLR","CRN",”PAAR”) & !is.null(reg.pat)){
      message(paste0("Custom REGEX patterns are not supported with the 'CRN' or 'RxLR' motif options.\n The package will use the default REGEX patterns used to search for ", motif, " motifs."))
      Sys.sleep(2)
    }
    for (i in 1:length(seq)){
      if (motif == "RxLR"){
        reg.pat <- "^\\w{10,40}\\w{1,96}R\\wLR\\w{1,40}eer"
      } else if (motif == "CRN"){
        reg.pat <- "^\\w{1,90}LFLAK\\w+"
      } else if (motif == "PAAR"){
        reg.pat <- “PAAR”
	} else if (motif == "custom"){
        if (is.null(reg.pat)){
          stop("No custom REGEX pattern found.\n The 'custom' option requires a mandatory REGEX pattern")
        } else {
        reg.pat <- reg.pat
        }
      }
      regex[[i]] <- unlist(gregexpr(seq[[i]], pattern = reg.pat, perl = T ,ignore.case = T))
    #percentage <- percentage + 1/length(seq)*100
    }
    regex <- as.data.frame(do.call(rbind, regex))
    regex$seq <- names(seq)
    regex <- regex[!regex$V1 < 0, ]
    regex <- sequence[seqinr::getName(sequence) %in% regex$seq]
    if (length(regex) == 0){
      stop(paste0("No ",motif, " sequences found."))
    }
    return(regex)
}
```

After including the new PAAR motif, the user can specify the PAAR motif as an option of the `motif` parameter:

```
# Loading the effectR package in R
library(“effectR”)

# Using the read.fasta function of the seqinr package to import the V. cholerae FASTA proteome file

fasta.file <- " V_cholerae_ATCC_39315.AA.fasta.gz”
ORF <- seqinr::read.fasta(fasta.file)

# Step 1 prediction: Predict proteins with the PAAR motif

REGEX <- regex.search(ORF, motif = "PAAR)

```

This customization will allow users to add any motif of interest to the *effectR* package by forking the github repository, adding the additional line into the regex.search function, adding the respective reference of the motif into the `@references` section of the R documentation within the `regex.search` function, and submitting a pull request to the package maintainers. The package maintainers will update the package, test the motif, and, if valid, add the motif to the effector.summary function before updating *effectR*. The currently available CRAN version of *effectR* only includes the RxLR and CRN motifs to facilitate the familiarization and engagement of the community with the package, but additional custom REGEX patterns will be added as the package is updated.

