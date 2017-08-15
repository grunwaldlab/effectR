#' Function to run the effectR package using a shiny app
#'
#' @description This function will launch an interactive web interface to run the effectR package functions to obtain effectors. It runs using the \pkg{shiny} R package.
#'
#' @details To sucessfully run this function the user will need to set the
#' @import shiny
#' @param mafft.path Local path of the MAFFT binary executable file
#' @param hmm.path Local path of the HMMER binaries
#' @export
shiny.effectR <- function(mafft.path="/usr/local/bin", hmm.path="/usr/local/bin"){
  if (file.exists(file.path(mafft.path, "mafft")) == F ){
    stop(paste0("mafft not found in ",mafft.path,"\nCheck your MAFFT installation path\n"))
  }
  if (file.exists(file.path(hmm.path, "hmmbuild")) == F ){
    stop(paste0("hmmbuild not found in ",mafft.path,"\nCheck your HMMER installation path\n"))
  }
  assign("mafft.path.shiny",mafft.path, envir = globalenv())
  assign("hmm.path.shiny",hmm.path, envir = globalenv())
  shiny::runApp(system.file("shiny",".",package = "effectR"))
  invisible(NULL)
}
