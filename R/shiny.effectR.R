#' Function to run the effectR package using a shiny app
#'
#' @description This function will launch an interactive web interface to run the effectR package functions to obtain effectors. It runs using the \pkg{shiny} R package.
#'
#' @details To successfully run this function the user will need to set the
#' @import shiny
#' @param mafft.path Local path of folder containing the MAFFT binary executable file or the executable file itself. If not specified, then MAFFT must be in the program search path.
#' @param hmm.path Local path of  folder containing the HMMER binaries.  If not specified, then HMMER executables must be in the program search path.
#' @export
#'
shiny.effectR <- function(mafft.path = NULL, hmm.path = NULL){
  pkg.env <- new.env()
  shiny::runApp(system.file("shiny",".",package = "effectR"))
  invisible(NULL)
}
