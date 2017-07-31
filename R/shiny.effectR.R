#' Function to run the effectR package using a shiny app
#'
#' @description This function will launch an interactive web interface to run the effectR package functions to obtain effectors. It runs using the \pkg{shiny} R package.
#'
#' @usage shiny.effectR()
#' @details To sucessfully run this function the user will need to set the
#' @import shiny

shiny.effectR <- function(mafft.path=mafft.path, hmm.path=mafft.path){
  if (file.exists(paste0(mafft.path, "/mafft")) == F ){
    stop(paste0("mafft not found in ",mafft.path,"\nCheck your MAFFT installation path\n"))
  }
  if (file.exists(paste0(hmm.path, "/hmmbuild")) == F ){
    stop(paste0("hmmbuild not found in ",mafft.path,"\nCheck your HMMER installation path\n"))
  }
  assign("mafft.path.shiny",mafft.path, envir = globalenv())
  assign("hmm.path.shiny",hmm.path, envir = globalenv())
  shiny::runApp(system.file("shiny",".",package = "effectR"))
  invisible(NULL)
}
