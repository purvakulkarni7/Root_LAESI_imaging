#' This function reads extracted spectra from a text file and writes a new text file with only mz and intensity values. This happens to all the files in the folder path provided
#'
#' @title Write new text file with mz and intensity values
#'
#' @param inputPath input folder path containing text files with data in a specific format
#' @param outputPath output folder path where all the text files will be writtem
#'
#' @return New text file with two columns - mz and intensity values separated by a tab
#'
#' @author Purva Kulkarni
#' @date 22 May 2017

exportMassSpectra <- function(inputPath, outputPath)
{
  fileList = list.files(path = inputPath, pattern = "*.txt", full.names = FALSE)

# Read all the files one by one and write the mass spectra to a new file in the output folder
for (i in 1:length(fileList))
{
  setwd(inputPath)
  values = read.table(fileList[i], skip = 4, sep = ",", colClasses=c(NA, NA, "NULL"))
  head(values)

  setwd(outputPath)
  fileName = toString(fileList[i])
  write.table(values, fileName, sep = " ", row.names = FALSE, col.names = FALSE)
}
}