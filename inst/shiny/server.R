# Load dependencies
require(shiny)
require(shinyjs)
require(seqinr)
require(ggplot2)
require(viridis)
require(reshape2)

# Max size limit:
options(shiny.maxRequestSize=70*4024^2)

# Code for refresh button
`%then%` <- shiny:::`%OR%`
jscode <- "shinyjs.refresh = function() { history.go(0); }"


### Additional scripts

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

# Get HMMER path
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



# Converts FASTA to STOCKHOLM
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

mafft.path <- dirname(get_mafft_path())
hmm.path <- dirname(get_hmmer_path("hmmsearch"))
mafft.path.shiny <- mafft.path
hmm.path.shiny <- hmm.path


# All variable names
original.seq <- c("AA.fasta")
file.name <- c("REGEX.fasta")
mafft.out.name <- c("MAFFT.fasta")
hmmbuild.out <- c("hmmbuild.hmm")
hmmsearch.out <- c("hmmsearch.txt")
mafft.path <- mafft.path.shiny
hmm.path <- hmm.path.shiny
original.dir <- getwd()

# TMP dir
tmp.dir <- tempdir()
setwd(tmp.dir)
##

##### SHINY APP START #####
shinyServer(function(input, output, session) {

  # File read in
  fasta <- reactive({
    file <- input$file1
    write.fasta(sequences = read.fasta(file$datapath), names = names(read.fasta(file$datapath)), file.out = "AA.fasta", open = "w", nbchar = 10000)
    read.fasta(file$datapath)
  })

  # UI output
  output$input_fasta <- renderPrint({
    file <- input$file1
    if(is.null(file))
      return("No FASTA file")
    else
    return (paste(length(fasta()), "Sequences read."))
  })

  # Refresh app button
  observeEvent(input$refresh, {
    js$refresh();
  })

########### Pipeline funtions ############

  #####
  # Step 1: REGEX
  ####

  # REGEX Script
  percentage <- 0
  regex.seq <- reactive({
    validate(
      need(is.null(file) == F, "No sequence file loaded") %then%
      need(file.exists("AA.fasta") != F, "No AA file found. Reload the web-app and try again")
    )
    seq <- lapply(fasta(), function (x) paste(unlist(x),collapse = ""))
    withProgress(message = "Searching for sequences that match the REGEX of interest", value=0, {
      regex <- list()
      for (i in 1:length(seq)){
        if (input$motif_sel == "RxLR"){
          regex[[i]] <- unlist(gregexpr(seq[[i]], pattern="^\\w{12,60}r\\wlr\\w{6,10}eer", perl = T,ignore.case = T))
        } else if (input$motif_sel == "CRN"){
          regex[[i]] <- unlist(gregexpr(seq[[i]], pattern="^\\w{1,90}LFLAK\\w+", perl = T,ignore.case = T))
        }
        percentage <- percentage + 1/length(seq)*100
        incProgress(1/length(seq), detail = paste0("Progress: ",round(percentage,2),"%"))
      }
      regex <- as.data.frame(do.call(rbind, regex))
      regex$seq <- names(seq)
    })
    write.fasta(sequences = fasta()[names(fasta()) %in% regex[!regex$V1 %in% -1, 2]], names = names(fasta()[names(fasta()) %in% regex[!regex$V1 %in% -1, 2]]), file.out = file.name, open = "w", nbchar = 10000)
    fasta()[names(fasta()) %in% regex[!regex$V1 %in% -1, 2]]
  })

  # REGEX UI interface
  observeEvent(input$regex_seach, {
    output$regex <- renderPrint({
      regex.seq()
      cat("Number of sequences retained by REGEX: ")
      cat(length(regex.seq()))
      cat("\n")
      session$sendCustomMessage("regex_ready", list(fileSize=floor(runif(1) * 10000)))

    })
  })

  output$downloadREGEX <- downloadHandler(
    filename <- function() {
      paste("REGEX", "fasta", sep=".")
    },
    content <- function(file) {
       file.copy("REGEX.fasta", file)
    }
  )


  #####
  # Step 2: HMM search
  ####

    # Script
  observeEvent(input$do, {
    # Step 2A:
    # MAFFT alignment
    output$mafft <- renderPrint({
      validate(
        need(length(regex.seq()) > 4, "Not enough sequences for HMM step. More than 4 sequences with REGEX required to perform HMM search")
      )
      withProgress(message = 'Performing MAFFT alignment', value = 100, {
                mafft.command <- c(get_mafft_path(mafft.path),
                           "--legacygappenalty",
                           "--genafpair",
                           "--maxiterate", "1000",
                           "--quiet",
                           file.name)
        system2(mafft.command, stdout = mafft.out.name )
      })
      cat("MAFFT aligment: Done!")
      cat("\n")
    })

    # Step 2B:
    # HMM build
    output$hmmer_press <- renderPrint({
      validate(
        need(file.exists(mafft.out.name) != F, "No HMM performed")
      )
      withProgress(message = 'Building HMM model...', value = 100, {
        unlink(file.path(tmp.dir, hmmbuild.out))
        if (Sys.info()[['sysname']] %in% "Windows"){
          fasta_to_stockholm(fasta.file = mafft.out.name)
          stock.name <- gsub(mafft.out.name, pattern = ".fasta", replacement = ".stockholm"
          )
          hmmbuild_command <- c(get_hmmer_path("hmmbuild.exe", hmm.path),
                                "--amino",
                                hmmbuild.out,
                                stock.name)
        } else {
          hmmbuild_command <- c(get_hmmer_path("hmmbuild", hmm.path),
                                "--amino",
                                hmmbuild.out,
                                mafft.out.name)
        }
        system2(hmmbuild_command, stdout = F)
        Sys.sleep(0.2)
        if (Sys.info()[['sysname']] %in% "Windows"){
          hmmpress_command <- c(get_hmmer_path("hmmpress.exe", hmm.path),
                                hmmbuild.out)
        } else {
          hmmpress_command <- c(get_hmmer_path("hmmpress", hmm.path),
                                hmmbuild.out)
        }
      system2(hmmpress_command)
      })
      cat("HMM model built!")
      cat("\n")
    })


    # Step 2C:
    # HMM search
    output$hmmer_search <- renderPrint({
      validate(
        need(file.exists("hmmbuild.hmm") != F, "No HMM performed")
      )
      withProgress(message = 'Performing HMM search...', value = 100, {
        if (Sys.info()[['sysname']] %in% "Windows"){
          hmmsearch_command <- c(get_hmmer_path("hmmsearch.exe", hmm.path),
                                 "-T", "0",
                                 "--tblout", hmmsearch.out,
                                 hmmbuild.out,
                                 original.seq)
        } else {
          hmmsearch_command <- c(get_hmmer_path("hmmsearch", hmm.path),
                                 "-T", "0",
                                 "--tblout", hmmsearch.out,
                                 hmmbuild.out,
                                 original.seq)
        }
        system2(hmmsearch_command)
      })
      Sys.sleep(0.2)
      cat("HMM search: Done!")
      cat("\n")
      session$sendCustomMessage("download_ready", list(fileSize=floor(runif(1) * 10000)))
     })
  })

  output$downloadHMM <- downloadHandler(
    filename <- function() {
      paste("HMM_files", "zip", sep=".")
    },
    content <- function(file) {
      fileConn<-file("Readme.txt")
      writeLines(c("Files used in HMM step:", "MAFFT.fasta: MAFFT alignment of REGEX candidates.","hmmbuild.hmm: HMM model of REGEX RxLR candidates", "hmmsearch.txt: Results of HMM search using input file"), fileConn)
      close(fileConn)
      system("zip HMM_data.zip MAFFT.fasta hmmsearch.txt hmmbuild.hmm Readme.txt")
      file.copy("HMM_data.zip", file)
    },
    contentType = "application/zip"
  )


  # OPTIONAL STEP: Web Logo
  observeEvent(input$logo, {
        validate(
          need(length(regex.seq()) > 4, "Not enough sequences to construct the HMM LOGO")
        )
        output$preImage <- renderPlot({
          withProgress(message = 'Constructing HMM logo...', value = 100, {
            # Plot
            if (Sys.info()[['sysname']] %in% "Windows"){
              hmm <- utils::read.table(hmmbuild.out, blank.lines.skip = T, skip = 14, sep = "", fill = T, stringsAsFactors = F)
            } else {
              hmm <- utils::read.table(hmmbuild.out, blank.lines.skip = T, skip = 16, sep = "", fill = T, stringsAsFactors = F)
            }
            colnames(hmm) <- hmm[1,]
            hmm <- hmm[-(1:5),]
            hmm[,1] <- as.numeric(hmm[,1])
            hmm <- hmm[hmm[,1]%%1==0, ]
            hmm <- hmm[!is.na(as.numeric(hmm[,2])), ]
            rownames(hmm) <- hmm[,1]
            hmm <- hmm[,-1]
            hmm <- data.frame(sapply(hmm, function(x) as.numeric(as.character(x))))
            hmm.sums <- apply(hmm,2,function (x) max(x)/x)
            hmm.sums <- apply(hmm.sums,2,function (x) x/sum(x))
            hmm.sums <- apply(hmm.sums, 2, function (x) x - mean(as.numeric(as.character(x))))
            hmm.sums[hmm.sums < 0] <- 0
            # Melt HMM
            hmm.m <- reshape2::melt(t(hmm.sums))
            colnames(hmm.m) <- c("element","position","bits")
            hmm.m$bits <- as.numeric(as.character(hmm.m$bits))
            p <- ggplot(hmm.m, aes(x=position, y=bits, fill=element)) + geom_bar(stat = "identity",position = "stack",width=1, alpha=0.5) + geom_text(aes(label=element, size=bits), position='stack') +  scale_fill_viridis(discrete=TRUE) + theme_bw() + ylab("Relative Frequency (bits)") + guides(fill=FALSE)
            print(p)
          })
          session$sendCustomMessage("logo_ready", list(fileSize=floor(runif(1) * 10000)))
        })
       })


  #####
  # Step 3: Additional steps
  ####

  # Step 3A: Validate and create final dataset
  final_fasta <- reactive({
    regex <- read.fasta("REGEX.fasta")
    if(file.exists("hmmsearch.txt")){
      hmm <- read.table("hmmsearch.txt", stringsAsFactors = F)
      combined.names <- unique(c(hmm$V1, names(regex)))
      fasta()[names(fasta()) %in% combined.names]
    } else {
      fasta()[names(fasta()) %in% names(regex)]
    }
  })

  # Step 3B: Combine datasets and create REGEX table
  # Table
  summary_table <- reactive (
    if(file.exists("hmmsearch.txt")){
    hmm <- read.delim("hmmsearch.txt", comment.char = "#", header = F, sep = "", blank.lines.skip = T, stringsAsFactors = F)[,1]
    num_hmm <- length(hmm)
    num_len <- length(regex.seq())
    combined.names <- unique(c(hmm, names(regex.seq())))
    num_com <- length(combined.names)
    final.df <- data.frame(num_len,num_hmm,num_com)
    colnames(final.df) <- c("REGEX","HMM","Total")
    final.df
  } else {
    num_len <- length(final_fasta())
    final.df <- data.frame(num_len)
    colnames(final.df) <- c("REGEX")
    final.df
  }
  )


  observeEvent(input$combine, {
    output$final_table <- renderTable({
      summary_table()
    })
  })
  # Step 3C: Motif Table
  motif_table <- reactive({
    sequences <- lapply(final_fasta(), function (x) paste(unlist(x),collapse = ""))
    # REGEX search of the motifs of interest:
    if (input$motif_sel == "RxLR"){
      rxlr.num <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="r\\wlr", perl = T,ignore.case = T))), function (x) length(x)))
      rxlr.motif <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="r\\wlr", perl = T,ignore.case = T))), function (x) paste(x,collapse = ",")))
      eer.num <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="[ED][ED][KR]", perl = T,ignore.case = T))), function (x) length(x)))
      eer.motif <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="[ED][ED][KR]", perl = T,ignore.case = T))), function (x) paste(x,collapse = ",")))
      motifs <- data.frame(seqinr::getName(final_fasta()),rxlr.num,rxlr.motif,eer.num,eer.motif, stringsAsFactors = F)
    } else if (input$motif_sel == "CRN") {
      lflak.num <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="lflak", perl = T,ignore.case = T))), function (x) length(x)))
      lflak.motif <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="lflak", perl = T,ignore.case = T))), function (x) paste(x,collapse = ",")))
      hvlv.num <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="hvlv", perl = T,ignore.case = T))), function (x) length(x)))
      hvlv.motif <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern="hvlv", perl = T,ignore.case = T))), function (x) paste(x,collapse = ",")))
      motifs <- data.frame(seqinr::getName(final_fasta()),lflak.num,lflak.motif,hvlv.num,hvlv.motif, stringsAsFactors = F)
    }
    motifs[motifs == -1] <- NA
    motifs <- data.frame(apply(motifs, 2, as.character), stringsAsFactors = F)
    motifs[,2][is.na(motifs[,3])] <- 0
    motifs[,4][is.na(motifs[,5])] <- 0
    motifs <- motifs[order(motifs[,2],motifs[,4],decreasing = T),]
    motifs$summary <- "No MOTIFS"
    motifs$summary[motifs[,2] > 0 & motifs[,4] > 0] <- "Complete"
    motifs$summary[motifs[,2] != 0 & motifs[,4] == 0] <- if (input$motif_sel == "RxLR"){c("Only RxLR motif")} else if (input$motif_sel == "CRN"){("Only LFLAK motif")}
    motifs$summary[motifs[,2] == 0 & motifs[,4] != 0] <- if (input$motif_sel == "RxLR"){c("Only EER motif")} else if (input$motif_sel == "CRN"){("Only HVLV motif")}
    if (input$motif_sel == "RxLR"){
      colnames(motifs) <- c("Sequence ID","RxLR number","RxLR position","EER number","EER position","MOTIF")
    } else if (input$motif_sel == "CRN") {
      colnames(motifs) <- c("Sequence ID","LFLAK number","LFLAK position","HVLV number","HVLV position","MOTIF")
    }
    motifs
  })

  # Render motif table
  observeEvent(input$motif, {
    output$motif_table <- renderTable({
      motif_table()
     })
    session$sendCustomMessage("motif_ready", list(fileSize=floor(runif(1) * 10000)))
  })

  # Download MOTIF table
  output$motifDownload <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = paste("Motif_table","txt",sep = "."),
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      # Write to a file specified by the 'file' argument
      write.table(x = motif_table(), file, sep = "\t", quote = F, row.names = F)
    })

  # MOTIF summary table
  observeEvent(input$motif_summ, {
    output$motif_summ_table <- renderTable({
      table(motif_table()$MOTIF)
      })
  })

  # Downloading final RXLR file
  output$downloadData <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = paste("Candidate_effectors","fasta",sep = "."),
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      # Write to a file specified by the 'file' argument
    write.fasta(sequences = final_fasta(), names = names(final_fasta()), nbchar = 10000, file.out = file)
    })
})

######### SHINY APP ENDS ##########
# EOF
