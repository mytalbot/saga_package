#' Create SAGA samples target file.
#'
#' \code{saga_gentargets} will automatically create a "SampleInformation.txt" file (SIF) if none is provided. If you chose to use this function, a folder
#' with sample *.txt files is needed. Each textfile must contain the full Agilent microarray information.
#'
#' \cr PLEASE NOTE: When the targets file is generated it cannot know about your sample batches. Batch defaults to No. for each sample is 11. You
#' can use any other number to indicate your batches BUT you'll have to change the numbers manually (i.e, in Excel or in a text editor).
#'
#' @param smplpath path to the saga data folder with the user samples.
#'
#' @return \code{targets} default target file for user samples.
#'
#' @importFrom utils write.table
#'
#' @export
#'
#'


saga_gentargets    <- function(smplpath){
  setwd(smplpath)

  if(!file.exists("SampleInformation.txt")){
    filelist            <- list.files(pattern = ".txt")
    elements            <- length(filelist)

    cnames              <- c("SAMPLE_ID", "Filename", "Batch", "Group", "Vector", "TrueLabel")

    targets <- NULL
    for (i in 1:elements){

      f1                <- filelist[i]
      f1                <- gsub("*.txt", "",f1)
      f1                <- paste("X",f1, sep="")

      f2                <- filelist[i]

      targets           <- rbind(targets, data.frame(f1, f2, 11, NA, NA,NA) )
    }
    colnames(targets)   <- cnames

    write.table(targets, file = "SampleInformation.txt", row.names = FALSE, sep = "\t", quote=FALSE)

  }else{
    print("There is already a 'SampleInformation.txt' file in the folder!")
  }

}











