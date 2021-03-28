getNormalizedMSPeaksA <- function(file, controlPeak, TICnormalization, scaleNormalization, peakElbow)
{
  cat("Processing 'Normal' mass spectrometry format data \n")
  raw.data.set <- lapply(seq_along(file),function(list_index, j){
    cat("Starting processing list element", list_index, "from", j[list_index], "\n")
    curr.spectra <- fread(j[list_index], header = T, verbose = F, stringsAsFactors = F)
    incProgress(amount = 1.0/length(file), message = paste("Starting processing file ", list_index, " from Set 1", sep = ""))
    #print(curr.spectra)
    if(dim(curr.spectra)[2] > 2)
    {
      cat("The wrong file format was selected. This looks like DESI-MSI IMAGING data. Please try the application again and select 'Imaging' file format.\n")
      incProgress(message = "The correct file format was not selected in the input module!")
      stopApp()
      return(NULL)
    }
    curr.spectra <- curr.spectra %>% na.omit()
    if(is.character(curr.spectra[1,1]))
    {
      cat("This looks like DESI-MSI IMAGING data. Please try the application again and select 'Imaging' file format.\n")
      incProgress(message = "The correct file format was not selected in the input module!")
      stopApp()
      return(NULL)
    }
    if(peakElbow=="Yes")
    {
      local_peak_indices = which(diff(sign(diff(as.numeric(unlist(curr.spectra[,2])))))==-2)+1
      curr.spectra <- curr.spectra[local_peak_indices,]
    } 
    if (TICnormalization == "Yes")
    {
      curr.spectra[,2] <- (curr.spectra[,2])/(sum(curr.spectra[,2]))
    }
    if(controlPeak != -1)
    {
      control.peak <- curr.spectra[trunc(curr.spectra[,1]*100)/100 == controlPeak,2]
      if (length(control.peak) > 0)
      {
        curr.spectra[,2] <- curr.spectra[,2]/control.peak
      }
    }
    if (scaleNormalization == "Yes")
    {
      curr.spectra[,2] <- curr.spectra[,2]/max(curr.spectra[,2])
    }
    curr.spectra$Centroid <- 0
    curr.spectra$Disease.State <- 0
    curr.spectra$Index <- list_index
    setcolorder(curr.spectra, c(5,1,2,3,4))
    colnames(curr.spectra) <- c("File.Name", "Mass","Intensity", "Centroid", "Disease.State")
    curr.spectra
  }, j = file)
  #print("in getnormalizedmspeaksA")
  #print(pryr::mem_used())
  cat("All data corresponding to the first group has been read and processed! \n")
  raw.data.set
} 
getNormalizedMSPeaksB <- function(file, controlPeak, numAFiles, TICnormalization, scaleNormalization, peakElbow)
{
  raw.data.set <- lapply(seq_along(file),function(list_index, j){
    cat("Starting processing list element", list_index, "from", j[list_index], "\n")
    curr.spectra <- fread(j[list_index], header = T, verbose = F, stringsAsFactors = F)
    incProgress(amount = 1.0/length(file), message = paste("Starting processing file ", list_index, " from Set 2", sep = ""))
    if(dim(curr.spectra)[2] > 2)
    {
      cat("The wrong file format was selected. This looks like DESI-MSI IMAGING data. Please try the application again and select 'Imaging' file format.\n")
      stopApp()
      return(NULL)
    }
    curr.spectra <- curr.spectra %>% na.omit()
    if(is.character(curr.spectra[1,1]))
    {
      cat("The wrong file format was selected. This looks like DESI-MSI IMAGING data. Please try the application again and select 'Imaging' file format.\n")
      stopApp()
      return(NULL)
    }
    if(peakElbow=="Yes")
    {
      local_peak_indices = which(diff(sign(diff(as.numeric(unlist(curr.spectra[,2])))))==-2)+1
      curr.spectra <- curr.spectra[local_peak_indices,]
    } 
    if (TICnormalization == "Yes")
    {
      curr.spectra[,2] <- (curr.spectra[,2])/(sum(curr.spectra[,2]))
    }
    if(controlPeak != -1)
    {
      control.peak <- curr.spectra[trunc(curr.spectra[,1]*100)/100 == controlPeak,2]
      if (length(control.peak) > 0)
      {
        curr.spectra[,2] <- curr.spectra[,2]/control.peak
      }
    }
    if (scaleNormalization == "Yes")
    {
      curr.spectra[,2] <- curr.spectra[,2]/max(curr.spectra[,2])
    }
    curr.spectra$Centroid <- 0
    curr.spectra$Disease.State <- 1 #main change between A and B
    curr.spectra$Index <- list_index + numAFiles #unique file number across entire dataset
    setcolorder(curr.spectra, c(5,1,2,3,4))
    colnames(curr.spectra) <- c("File.Name", "Mass","Intensity", "Centroid", "Disease.State")
    curr.spectra
  }, j = file)
  #print("in getnormalizedmspeaksB")
  #print(pryr::mem_used())
  cat("All data corresponding to both groups of data have been read and processed! \n")
  raw.data.set
}  