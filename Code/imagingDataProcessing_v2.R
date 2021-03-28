library(shiny)
library(stringr)
library(data.table)
require(class)
library(tidyverse)
options(datatable.fread.datatable=FALSE)
numextract <- function(string){ 
  as.numeric(str_extract(string, "\\-*\\d+\\.*\\d*"))
} 
extractPatientNum <- function(string) {
  as.numeric(sapply(str_extract_all(string, "\\d+"),tail,1))
}
getNormalizedImagingPeaksA <- function(fileNames, file, controlPeak, TICnormalization, scaleNormalization, peakElbow)
{
  cat("Processing DESI-MSI imaging data \n")
  unflat.raw.datA <- lapply(seq_along(file),function(list_index, j){
    cat("Starting processing list element", list_index, "from", j[list_index], "\n")
    patient.num = numextract(fileNames$name[list_index])
    if(is.na(patient.num)) {
      cat("Not all files are named correctly. In other words,", j[list_index], "is not named appropriately for the application. \n")
      incProgress(message = "Not all filenames had a patient number!")
      stopApp()
      return(NULL)
    }
    setProgress(value = 0.1, message = paste("Starting processing patient ", patient.num, "from set A", "\n", sep = " "))
    raw.file <- fread(j[list_index], header = F, verbose = F, stringsAsFactors = F)
    #if(!is.numeric(raw.file[6,1]))
    #{
    #  cat("File has extra header rows. Please check application input module for sample file format!\n")
    #  stopApp()
    #  return(NULL)
    #}
    fileHeader <- raw.file[1:4,]
    raw.file <- raw.file[-c(1:5),] 
    rownames(raw.file) <- NULL
    start <- seq(1, by = 2, length = ncol(fileHeader) / 2)
    end <- (last(start)+1)/2
    #setProgress(message=paste("Patient processed", patient.num, sep = " "))
    sameFileSpectraList <- lapply(start, function(i, curr.df) {
      incProgress(amount = 1.0/length(start), message =  paste("Processing spectra ", (i+1)/2, " of ", end, " from patient ", patient.num, " (set A)", sep = ""))
      cat("Processing spectra", (i+1)/2, "of", end, "from", j[list_index], "\n")
      curr.spectra <- data.table(Mass = as.numeric(unlist(curr.df[,i])), Relative.Intensity = as.numeric(unlist(curr.df[,(i+1)])))
      curr.spectra <- curr.spectra %>% na.omit()
      if(peakElbow=="Yes")
      {
        curr.spectra <- curr.spectra[which(diff(sign(diff(as.numeric(unlist(curr.spectra[,2])))))==-2)+1,]
      } 
      if (TICnormalization == "Yes") 
      {
        curr.spectra[,2] <- curr.spectra[,2]/(sum(curr.spectra[,2]))
      }
      if (scaleNormalization == "Yes")
      {
        curr.spectra[,2] <- curr.spectra[,2]/max(curr.spectra[,2])
      }
      curr.spectra$Centroid <- 0
      curr.spectra$Scan.Index = (i+1)/2
      curr.spectra$Patient.Index = patient.num
      scan.info <- data.table(Scan = as.numeric(fileHeader[1,i:(i+1)][,2]), X = as.numeric(fileHeader[2,i:(i+1)][,2]),
                              Y = as.numeric(fileHeader[3,i:(i+1)][,2]), Disease.State = 0, Patient = patient.num)
      list(curr.spectra, scan.info)
    },curr.df = raw.file)
    sameFileSpectraList
  }, j = file)
  raw.data.setA <- unlist(unflat.raw.datA,recursive=F)
  rm(unflat.raw.datA)
  raw.data.setA = do.call(Map, c(f = rbind, raw.data.setA))
  cat("All data corresponding to the first group has been read and processed! \n")
  raw.data.setA
}
getNormalizedImagingPeaksB <- function(fileNames, file, controlPeak, TICnormalization, scaleNormalization, peakElbow)
{
  unflat.raw.datB <- lapply(seq_along(file),function(list_index, j){
    patient.num = numextract(fileNames$name[list_index])
    if(is.na(patient.num)) {
      cat("Not all files are named correctly. In other words,", j[list_index], "is not named appropriately for the application. \n")
      stopApp()
      return(NULL)
    }
    cat("Starting processing list element", list_index, "from", j[list_index], "\n")
    setProgress(value = 0.1, message = paste("Starting processing patient", patient.num, "from set B", "\n", sep = " "))
    raw.file <- fread(j[list_index], header = F, verbose = F, stringsAsFactors = F)
    #if(!is.numeric(raw.file[6,1]))
    #{
    #  cat("File has extra header rows. Please check application input module for sample file format!\n")
    #  stopApp()
    #  return(NULL)
    #}
    fileHeader <- raw.file[1:4,]
    raw.file <- raw.file[-c(1:5),] 
    rownames(raw.file) <- NULL
    start <- seq(1, by = 2, length = ncol(fileHeader) / 2)
    end <- (last(start)+1)/2
    sameFileSpectraList <- lapply(start, function(i, curr.df) {
      incProgress(amount = 1.0/length(start), message =  paste("Processing spectra ", (i+1)/2, " of ", end, " from patient ", patient.num, " (set B)", sep = ""))
      cat("Processing spectra", (i+1)/2, "of", end, "from", j[list_index], "\n")
      curr.spectra <- data.table(Mass = as.numeric(unlist(curr.df[,i])), Relative.Intensity = as.numeric(unlist(curr.df[,(i+1)])))
      curr.spectra <- curr.spectra %>% na.omit()
      if(peakElbow=="Yes")
      {
        curr.spectra <- curr.spectra[which(diff(sign(diff(as.numeric(unlist(curr.spectra[,2])))))==-2)+1,]
      } 
      #curr.spectra <- curr.spectra[which(diff(sign(diff(as.numeric(unlist(curr.spectra[,2])))))==-2)+1,]
      if (TICnormalization == "Yes") 
      {
        curr.spectra[,2] <- curr.spectra[,2]/(sum(curr.spectra[,2]))
      }
      if (scaleNormalization == "Yes")
      {
        curr.spectra[,2] <- curr.spectra[,2]/max(curr.spectra[,2])
      }
      curr.spectra$Centroid <- 0
      curr.spectra$Scan.Index = (i+1)/2
      curr.spectra$Patient.Index = patient.num
      scan.info <- data.table(Scan = as.numeric(fileHeader[1,i:(i+1)][,2]), X = as.numeric(fileHeader[2,i:(i+1)][,2]),
                              Y = as.numeric(fileHeader[3,i:(i+1)][,2]), Disease.State = 1, Patient = patient.num)
      list(curr.spectra, scan.info)
    },curr.df = raw.file)
    sameFileSpectraList
  }, j = file)
  raw.data.setB <- unlist(unflat.raw.datB,recursive=F)
  rm(unflat.raw.datB)
  raw.data.setB = do.call(Map, c(f = rbind, raw.data.setB))
  cat("All data corresponding to the both groups have been read and processed! \n")
  raw.data.setB
}