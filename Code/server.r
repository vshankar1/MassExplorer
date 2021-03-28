#Vishnu Shankar
#Mass Spectrometry Online Tool
#July 10, 2020
set.seed(123)
#setwd("~/Dropbox/Zare_Project/Week4_Update/MassSpecTool/")
#options(repos = BiocInstaller::biocinstallRepos())
library(shiny)
library(shinydashboard)
library(ggplot2)
library(plotly)
require(class)
library(data.table)
library(samr)
library(DT)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(glmnet)
source(file = "normalDataProcessing_v2.R", local = T)
source(file = "imagingDataProcessing_v2.R", local = T)
library(devtools)
library(roxygen2)
library(fcluster)
#devtools::load_all(pkg = "../fcluster", recompile = T, quiet = T)
#print("after loading packages")
#print(pryr::mem_used())
#fcluster_path <- inline::getDynLib("fcluster")
#dyn.load(unlist(fcluster_path[2]))

#devtools::load_all(fcluster_path)

options(shiny.maxRequestSize=1000000000*1024^2) 
shinyServer(function(input, output, session) {
  #Display Table in Instructions Panel  
  session$onSessionEnded(stopApp)

  sample.csv.data <- data.frame(mass = c(50, 51, 52, 53, 54, 55, 57, 59), Intensity = c(0, 2.566398,
                                  5.055928, 5.678105, 5.522893, 4.351526, 3.11324, 0))
  colnames(sample.csv.data)[1] <- "m/z"
  sample.csv2.data <- data.frame(a= c("Scan", "X", "Y", "Z", "m/z", 50.13077, 50.130837), 
                                 b= c(360, 3.8, 2.4, 1, "Intensity", 443.581421, 1167.427856),
                                 c= c("Scan", "X", "Y", "Z", "m/z", 50.120516, 50.120583),
                                 d= c(361, 4, 2.4, 1, "Intensity", 553.307312, 1084.705078))
  names(sample.csv2.data) <- NULL
  output$tbl <- renderTable({sample.csv.data}, spacing = 'xs')
  output$tbl2 <- renderTable({sample.csv2.data}, spacing = 'xs')
  
  globalPath<-paste0("~/Downloads/MassExplorer_Output/")
  dir.create("~/Downloads/MassExplorer_Output")
  sink(paste0(globalPath, "log.txt"))
  
  #Pairwise Comparison Module: What are the inputted file choices for A and B.
  output$choicesFileA <- renderUI({
    if (is.null(input$file1)) 
    {
      # User has not uploaded a file yet
      return(NULL)
    }
    if(input$fileType == "Normal")
    {
      selectInput("fileAInput", "Inputted Spectra",
                  input$file1[1])
    } else {
      fileAdata <- data1()
      fileANames <- paste("Scan:", fileAdata[[2]]$Scan, ";  X:", fileAdata[[2]]$X, 
            ";  Y:", fileAdata[[2]]$Y, ";  Patient:", fileAdata[[2]]$Patient, sep = "")
      selectInput("fileAScan", "Inputted Spectra", fileANames)
    }
  })
  output$choicesFileB <- renderUI({
    if (is.null(input$file2)) 
    {
      # User has not uploaded a file yet
      return(NULL)
    }
    if (input$fileType == "Normal")
    {
      selectInput("fileBInput", "Inputted Spectra",
                  input$file2[1])
    } else {
      fileBdata <- data2()
      fileBNames <- paste("Scan:", fileBdata[[2]]$Scan, ";  X:", fileBdata[[2]]$X, 
                          ";  Y:", fileBdata[[2]]$Y, ";  Patient:", fileBdata[[2]]$Patient, sep = "")
      selectInput(inputId = "fileBScan",label = "Inputted Spectra", 
                  choices = fileBNames, selected = fileBNames[2])
    }
  })
  
  
  #Display plots based on inputted files
  #Data1 and data2 load and read in the files
  data1 <- reactive({
    withProgress(message = "Processing Samples from the First Group", {
      #To load the data, we use the following method
      #We get the list of files and the inputted control peak from the user
      #If the user does not upload any files, we do not proceed.
      #We then check if the format is normal or imaging, which then selects the appropriate dataProcessing method
      #After processing the data, we label the list elements by fileName, so it is clear what element corresponds to what file
      file1 = input$file1
      controlPeak = input$controlPeak
      TICnormalization = input$ticNormalization
      scaleNormalization = input$scaleNormalization
      peakElbow = input$peakElbow
      if (is.null(file1))
        return(NULL) #user has not uploaded files yet
      data <- c()
      if (input$fileType=="Normal")
      {
        data <- getNormalizedMSPeaksA(file1$datapath, controlPeak, TICnormalization, scaleNormalization, peakElbow)
        fileNames <- input$file1[1]
        if(is.null(fileNames)){return(NULL)}
        myNames <- fileNames[,1]
        names(data) <- myNames
        data
      } else {
        data <- getNormalizedImagingPeaksA(input$file1, file1$datapath, controlPeak, TICnormalization, scaleNormalization, peakElbow)
        data
      }
    })
  })
  data2 <- reactive({
    withProgress(message = "Processing Samples from the Second Group", {
      file2 = input$file2
      controlPeak = input$controlPeak
      TICnormalization = input$ticNormalization
      scaleNormalization = input$scaleNormalization
      peakElbow = input$peakElbow
      if (is.null(file2))
        return(NULL) #user has not uploaded files yet
      data <- c()
      numAFiles <- dim(input$file1[1])[1]
      if (input$fileType=="Normal")
      {
        #we pass in numAFiles, so that we can set up a unique index for files from set B to define the dim of training.mat
        data <- getNormalizedMSPeaksB(file2$datapath, controlPeak, numAFiles, TICnormalization, scaleNormalization, peakElbow)
        fileNames <- input$file2[1]
        if(is.null(fileNames)){return(NULL)}
        myNames <- fileNames[,1]
        names(data) <- myNames
        data
      } else {
        data <- getNormalizedImagingPeaksB(input$file2, file2$datapath, controlPeak, TICnormalization, scaleNormalization, peakElbow)
        #print(pryr::mem_used())
        data
      }
    })
  })
  #We then load the data, normalize by total ion current, isolate local peaks, then plot
  #using ggplot.
  output$plotSetA <- renderPlotly({
    #This method works as follows
    #1. Use reactive method to load a cache of inputted data as a list
    #2. We then use the drop-down menu to get the selected file name
    #3. Check if the dropdown ui has been rendered, to avoid any error
    #4. We then select the list element that matches with the selected name
    #5. Since ggplot takes in a dataframe, we create a fileAdisplay, which takes a mass, intensity
    #5b. We subset according to how the user has selected the range on the slider to get the modified data frame
    #6. We then use ggplot to create the displayed mass-spectra distribution
    fileAdata <- data1()
    if (input$fileType == "Normal")
    {
      selectedFile <- input$fileAInput
      if(is.null(selectedFile)){return(NULL)}
      selectedIndex <- which(selectedFile == names(fileAdata))
      fileAdisplay <- fileAdata[[selectedIndex]]
      sliderRange <- input$mz_range
      sliderLower <- sliderRange[1]
      sliderUpper <- sliderRange[2]
      fileAdisplay <- data.frame(mass = fileAdisplay[,2], intensity = fileAdisplay[,3], isotope.mass = fileAdisplay[,2], file.num = fileAdisplay[,1])
      selectedDisplay <- fileAdisplay[which(fileAdisplay$mass >= sliderLower & fileAdisplay$mass <= sliderUpper),]
      plot1 <- ggplot(selectedDisplay, aes(x = mass, y = intensity, label = file.num)) + 
        geom_segment(aes(xend=isotope.mass, yend = 0), colour="red") + 
        ggtitle(paste("Mass Spectra from " , input$fileAInput ," (", input$labelSet1, ")", sep = "")) + 
        xlab("m/z") + 
        ylab("Relative Intensity") 
      ggplotly(plot1, tooltip = c("isotope.mass","intensity"))
    } else {
      selectedSpectra <- input$fileAScan
      #get the patient, then the index
      patient_num = extractPatientNum(selectedSpectra)
      selectedPatientIndex = which(fileAdata[[2]]$Patient == patient_num)
      fileAdataNames <- paste("Scan:", fileAdata[[2]]$Scan[selectedPatientIndex], ";  X:", fileAdata[[2]]$X[selectedPatientIndex], 
                              ";  Y:", fileAdata[[2]]$Y[selectedPatientIndex], ";  Patient:", fileAdata[[2]]$Patient[selectedPatientIndex], sep = "")
      #cat("\n----------Patient Num ---------\n")
      #print(patient_num)
      #cat("\n----------Selected Spectra ---------\n")
      #print(selectedSpectra)
      #cat("\n --------- fileAdataNames------------\n")      
      #print(fileAdataNames)
      selectedIndex <- which(selectedSpectra == fileAdataNames)[1]
      tableIndex <- which(fileAdata[[1]]$Scan.Index == selectedIndex)
      fileAdisplay <- data.frame(mass = fileAdata[[1]][tableIndex,1], intensity = fileAdata[[1]][tableIndex,2], isotope.mass = fileAdata[[1]][tableIndex,1])
      colnames(fileAdisplay) <- c("mass" , "intensity", "isotope.mass")
      sliderRange <- input$mz_range
      sliderLower <- sliderRange[1]
      sliderUpper <- sliderRange[2]
      selectedDisplay <- fileAdisplay[which(fileAdisplay$mass >= sliderLower & fileAdisplay$mass <= sliderUpper),]
      plot1 <- ggplot(selectedDisplay, aes(x = mass, y = intensity)) + 
        geom_segment(aes(xend=isotope.mass, yend = 0), colour="red") + 
        ggtitle(paste("Mass Spectra from " , fileAdataNames[selectedIndex], sep = "")) + 
        xlab("m/z") + 
        ylab("Relative Intensity") 
      ggplotly(plot1, tooltip = c("isotope.mass","intensity"))
    }
  })
  output$plotSetB <- renderPlotly({
    fileBdata <- data2()
    if (input$fileType == "Normal")
    {
      selectedFile <- input$fileBInput
      if(is.null(selectedFile)){return(NULL)}
      selectedIndex <- which(selectedFile == names(fileBdata))
      fileBdisplay <- fileBdata[[selectedIndex]]
      fileBdisplay <- data.frame(mass = fileBdisplay[,2], intensity = fileBdisplay[,3], isotope.mass = fileBdisplay[,2], file.num = fileBdisplay[,1])
      sliderRange <- input$mz_range
      sliderLower <- sliderRange[1]
      sliderUpper <- sliderRange[2]
      selectedDisplay <- fileBdisplay[which(fileBdisplay$mass >= sliderLower & fileBdisplay$mass <= sliderUpper),]
      plot1 <- ggplot(selectedDisplay, aes(x = mass, y = intensity, label = file.num)) + 
        geom_segment(aes(xend=isotope.mass, yend = 0), colour="blue") + 
        ggtitle(paste("Mass Spectra from ", input$fileBInput, " (", input$labelSet2, ")", sep = "")) + 
        xlab("m/z") + 
        ylab("Relative Intensity") 
      ggplotly(plot1, tooltip = c("isotope.mass","intensity"))
    } else {
      selectedSpectra <- input$fileBScan
      patient_num = extractPatientNum(selectedSpectra)
      selectedPatientIndex = which(fileBdata[[2]]$Patient == patient_num)
      fileBdataNames <- paste("Scan:", fileBdata[[2]]$Scan[selectedPatientIndex], ";  X:", fileBdata[[2]]$X[selectedPatientIndex], 
                              ";  Y:", fileBdata[[2]]$Y[selectedPatientIndex], ";  Patient:", fileBdata[[2]]$Patient[selectedPatientIndex], sep = "")
      selectedIndex <- which(selectedSpectra == fileBdataNames)[1]
      tableIndex <- which(fileBdata[[1]]$Scan.Index == selectedIndex)
      fileBdisplay <- data.frame(mass = fileBdata[[1]][tableIndex,1], intensity = fileBdata[[1]][tableIndex,2], isotope.mass = fileBdata[[1]][tableIndex,1])
      colnames(fileBdisplay) <- c("mass" , "intensity", "isotope.mass")
      sliderRange <- input$mz_range
      sliderLower <- sliderRange[1]
      sliderUpper <- sliderRange[2]
      selectedDisplay <- fileBdisplay[which(fileBdisplay$mass >= sliderLower & fileBdisplay$mass <= sliderUpper),]
      plot1 <- ggplot(selectedDisplay, aes(x = mass, y = intensity)) + 
        geom_segment(aes(xend=isotope.mass, yend = 0), colour="blue") + 
        ggtitle(paste("Mass Spectra from " , fileBdataNames[selectedIndex], sep = "")) + 
        xlab("m/z") + 
        ylab("Relative Intensity") 
      ggplotly(plot1, tooltip = c("isotope.mass","intensity"))
    }
  }) 
  output$mytableA = DT::renderDT({
    #dataTable, unlike ggplot can handle matrices
    #therefore, we can just subset the matrix in R and return the mass and intensity
    fileAdata <- data1()
    if (input$fileType == "Normal")
    {
      selectedFile <- input$fileAInput
      if(is.null(selectedFile)){return(NULL)}
      selectedIndex <- which(selectedFile == names(fileAdata))
      fileAdisplay <- fileAdata[[selectedIndex]]
      sliderRange <- input$mz_range
      sliderLower <- sliderRange[1]
      sliderUpper <- sliderRange[2]
      selectedData <- fileAdisplay[fileAdisplay[,2] > sliderLower & fileAdisplay[,2] < sliderUpper,]
      #selectedData <- selectedData[order(selectedData[,3], decreasing = T),]
      selectedData[,2:3] #only mass, intensity needed
    } else {
      selectedSpectra <- input$fileAScan
      patient_num = extractPatientNum(selectedSpectra)
      selectedPatientIndex = which(fileAdata[[2]]$Patient == patient_num)
      fileAdataNames <- paste("Scan:", fileAdata[[2]]$Scan[selectedPatientIndex], ";  X:", fileAdata[[2]]$X[selectedPatientIndex], 
                              ";  Y:", fileAdata[[2]]$Y[selectedPatientIndex], ";  Patient:", fileAdata[[2]]$Patient[selectedPatientIndex], sep = "")
      selectedIndex <- which(selectedSpectra == fileAdataNames)[1]
      tableIndex <- which(fileAdata[[1]]$Scan.Index == selectedIndex)
      fileAdisplay <- data.frame(mass = fileAdata[[1]][tableIndex,1], intensity = fileAdata[[1]][tableIndex,2])
      colnames(fileAdisplay) <- c("m/z" , "Intensity")
      sliderRange <- input$mz_range
      sliderLower <- sliderRange[1]
      sliderUpper <- sliderRange[2]
      selectedData <- fileAdisplay[fileAdisplay[,1] > sliderLower & fileAdisplay[,1] < sliderUpper,]
      #selectedData <- selectedData[order(selectedData[,2], decreasing = T),]
      selectedData
    }
  }, rownames = F, options = list(initComplete = JS('function(setting, json) {}')))
  output$mytableB = DT::renderDT({
    fileBdata <- data2()
    if (input$fileType == "Normal")
    {
      selectedFile <- input$fileBInput
      if(is.null(selectedFile)){return(NULL)}
      selectedIndex <- which(selectedFile == names(fileBdata))
      fileBdisplay <- fileBdata[[selectedIndex]]
      sliderRange <- input$mz_range
      sliderLower <- sliderRange[1]
      sliderUpper <- sliderRange[2]
      selectedData <- fileBdisplay[fileBdisplay[,2] > sliderLower & fileBdisplay[,2] < sliderUpper,]
      #selectedData <- selectedData[order(selectedData[,3], decreasing = T),]
      selectedData[,2:3] #only mass, intensity needed
    } else {
      selectedSpectra <- input$fileBScan
      patient_num = extractPatientNum(selectedSpectra)
      selectedPatientIndex = which(fileBdata[[2]]$Patient == patient_num)
      fileBdataNames <- paste("Scan:", fileBdata[[2]]$Scan[selectedPatientIndex], ";  X:", fileBdata[[2]]$X[selectedPatientIndex], 
                              ";  Y:", fileBdata[[2]]$Y[selectedPatientIndex], ";  Patient:", fileBdata[[2]]$Patient[selectedPatientIndex], sep = "")
      selectedIndex <- which(selectedSpectra == fileBdataNames)[1]
      tableIndex <- which(fileBdata[[1]]$Scan.Index == selectedIndex)
      fileBdisplay <- data.frame(mass = fileBdata[[1]][tableIndex,1], intensity = fileBdata[[1]][tableIndex,2])
      colnames(fileBdisplay) <- c("m/z" , "Intensity")
      sliderRange <- input$mz_range
      sliderLower <- sliderRange[1]
      sliderUpper <- sliderRange[2]
      selectedData <- fileBdisplay[fileBdisplay[,1] > sliderLower & fileBdisplay[,1] < sliderUpper,]
      #selectedData <- selectedData[order(selectedData[,2], decreasing = T),]
      selectedData
    }
  }, rownames = F)
  
  processedImagingData <- reactive({
    fileAdata <- data1()
    fileBdata <- data2()
    if (is.null(fileAdata))
      return(NULL)
    if (is.null(fileBdata))
      return(NULL)
    raw.peaks = c(fileAdata[[1]]$Mass, fileBdata[[1]]$Mass)
    saveRDS(fileAdata, file = paste0(globalPath, "fileA_data_raw"))
    saveRDS(fileBdata, file = paste0(globalPath, "fileB_data_raw"))
    cat("Raw processed data successfully saved (fileA_data_raw, fileB_data_raw).!\n")
    cat("\nClustering detected peaks from all samples. \n")
    incProgress(0.2, message = "Clustering Peaks Across All Samples")
    hier_cluster <- fcluster(raw.peaks, method = 'complete')
    cut_tree <- fcutree(raw.peaks, hier_cluster, h=input$error)
    centroids <- data.frame(cut_tree$centroid)
    saveRDS(centroids[,1], file = paste0(globalPath, "1_centroid_pks"))
    cat("Calculated centroids successfully saved (1_centroid_pks)!\n")
    incProgress(0.2, message = "Gathering peaks and initiating calculations")
    cat("Clustering in progress, gathering peaks, and initiating calculations. This is the most expensive part of the computation, since we are currently assigning each peak to the nearest centroid.\n")
    knn.out <- knn1(centroids[,1], data.frame(raw.peaks), 1:length(centroids[,1]))
    closest_centroid = centroids[knn.out,]
    saveRDS(knn.out, file = paste0(globalPath, "1_centroid_assignment"))
    cat("Centroid Assignment successfully saved (1_centroid_assignment)!\n")
    num_rows_set_A = nrow(fileAdata[[1]])
    fileAdata[[1]]$Centroid = closest_centroid[1:num_rows_set_A]
    fileBdata[[1]]$Centroid = closest_centroid[(num_rows_set_A+1):length(closest_centroid)]
    saveRDS(fileAdata, file = paste0(globalPath, "fileA_data_raw"))  
    saveRDS(fileBdata, file = paste0(globalPath, "fileB_data_raw"))
    cat("Raw data with calculated centroids successfully saved (fileA_data_raw, fileB_data_raw)!\n")
    rm(closest_centroid)
    rm(knn.out)
    rm(cut_tree)
    rm(centroids)
    rm(hier_cluster)
    incProgress(0.2, message = "Gathering Peaks from Each Scan")
    
    fileAdatalist=split(fileAdata[[1]], by = "Patient.Index") 
    #saveRDS(fileAdatalist, file = "~/Downloads/MassExplorer_R_Objects/fileAdatalist")
    #fileAdatalist=readRDS(file="~/Downloads/MassExplorer_R_Objects/fileAdatalist")
    fileA_datatable_list=lapply(fileAdatalist, function(x){
      data.table::dcast.data.table(x, formula = Patient.Index + Scan.Index ~ Centroid, 
                       value.var = "Relative.Intensity", fun.aggregate = max, fill = 0)
    })
    fileAdata[[1]] = data.table::rbindlist(fileA_datatable_list, use.names = T, fill = T)
    rm(fileA_datatable_list)
    rm(fileAdatalist)
    
    #saveRDS(fileA_datatable, file = "~/Downloads/MassExplorer_R_Objects/2_fileA_datatable")
    
    #fileAdata[[1]] = dcast.data.table(fileAdata[[1]], formula = Patient.Index + Scan.Index ~ Centroid, 
    #                                  value.var = "Relative.Intensity", fun.aggregate = max, fill = 0)
    saveRDS(fileAdata, file = paste0(globalPath, "2_fileA_training_dat"))
    cat("Training matrix built for fileA data (2_fileA_training_dat)!\n")
    fileBdatalist=split(fileBdata[[1]], by = "Patient.Index") 
    fileB_datatable_list=lapply(fileBdatalist, function(x){
      data.table::dcast.data.table(x, formula = Patient.Index + Scan.Index ~ Centroid, 
                                   value.var = "Relative.Intensity", fun.aggregate = max, fill = 0)
    })
    fileBdata[[1]] = data.table::rbindlist(fileB_datatable_list, use.names = T, fill = T)
    rm(fileB_datatable_list)
    rm(fileBdatalist)
    #fileBdata[[1]] = dcast.data.table(fileBdata[[1]], formula = Patient.Index + Scan.Index ~ Centroid, 
    #                                  value.var = "Relative.Intensity", fun.aggregate = max, fill = 0)
    saveRDS(fileBdata, file = paste0(globalPath, "2_fileB_training_dat"))
    cat("Training matrix built for fileB data (2_fileB_training_dat)!\n")
    
    training.table = rbindlist(list(fileAdata[[1]], fileBdata[[1]]), use.names = T, fill = T)
    descriptor.mat = rbindlist(list(fileAdata[[2]], fileBdata[[2]]), use.names = T, fill = T)
    saveRDS(descriptor.mat, file = paste0(globalPath, "3_descriptor.table"))
    cat("Successfully saved Scan, X, Y, Patient, Disease.State information (3_descriptor.table)!\n")
    
    na.cols = unique(which(is.na(training.table), arr.ind=TRUE)[,2])
    if(length(na.cols) > 0)
    {
      training.table[,(na.cols) := lapply(.SD, function(x){out <- x; out[is.na(out)] <- 0; out}), .SDcols = na.cols]
    }
    training.table[ ,`:=`(Scan.Index = NULL, Patient.Index = NULL)]
    saveRDS(training.table, file = paste0(globalPath, "3_training.table"))
    cat("Successfully saved training matrix (3_training_table)!\n")
    output <- list(training.table, descriptor.mat)
    sink(NULL)
    closeAllConnections()
    output
  })
  
  getDataMat <- reactive({
    #to get the training matrix
    #we get a list of the raw data for inputted set A. This consists of file.num, mass, intensity, centroid, disease.state
    #we get a similar list for inputted set b
    #we combine the list to form a training.matrix
    #Using the hierarchical clustering method, we separate the data into clusters, based on the inputted error
    #We then build a training matrix, with the columns as the peaks and the rows representing each "patient"
    withProgress(message = 'Constructing Training Matrix', {
      fileAdata <- data1()
      fileBdata <- data2()
      if (is.null(fileAdata))
        return(NULL)
      if (is.null(fileBdata))
        return(NULL)
      if (input$fileType == "Normal") 
      {
        num_spectra_A = length(fileAdata)
        num_spectra_B = length(fileBdata)
        training.df <- rbindlist(c(fileAdata, fileBdata))
        rm(fileAdata)
        rm(fileBdata)
        #SAM or t-statistics
        plotType <- "T"
        if (min(num_spectra_A, num_spectra_B) > 2)
          plotType = "S"
        incProgress(0.2, message = "Clustering Peaks Across All Samples")
        cat("Clustering peaks across all samples!\n")
        hier_cluster <- fcluster(training.df$Mass, method = 'complete')
        cut_tree <- fcutree(training.df$Mass, hier_cluster, h=input$error)
        rm(hier_cluster)
        centroids <- data.frame(cut_tree$centroid)
        incProgress(0.2, message = "Clustering completed, gathering peaks, and initiating calculations")
        cat("Clustering in progress, gathering peaks, and initiating calculations. This is the most expensive part of the computation, since we are currently assigning each peak to the nearest centroid.\n")
        knn.out <- knn1(centroids[,1], data.frame(training.df$Mass), 1:length(centroids[,1]))
        training.df$Centroid <- centroids[knn.out, ]
        rm(centroids)
        rm(knn.out)
        output <- list(training.df, plotType, input$fileType)
        #we want to be able to retrieve the following:
        #(1) training data.table
        #(2) whether to do t-statistic or SAM calculation
        #(3) What is the imaging type: normal or imaging
      } else {
        plotType = "S"
        processedList <- processedImagingData()
        training.table<- processedList[[1]]
        descriptor.mat <- processedList[[2]]
        rm(processedList)
        incProgress(0.2, message = "Setting up calculations based on gathered peaks")
        cat(input$threshold)
        if (input$threshold != -1 && input$threshold > 0 && input$threshold < 1)
        {
          #print('removing sparse metabolites (imaging mode)')
          cat("Removing sparse metabolites, based on specifications in input module!\n")
          fraction_nonzero = colSums(training.table != 0)/nrow(training.table)
          selected_ind = as.vector(which(fraction_nonzero >= input$threshold))
          training.table <- training.table[, selected_ind, with = F]
          #print(dim(training.table))
        }
        output <- list(training.table, descriptor.mat, plotType, input$fileType)
      }
    })
  })
  
  getLassoOutput <- reactive({
    withProgress(message = 'Running Lasso', {
      if(input$fileType == "Normal")
      {
        training_X <- getTrainingXMatrix()
        training_Y <- training_X$Disease.State
        training_X <- scale(as.matrix(training_X[,-c(1:2)]), center = T)
        if (length(training_Y) >= 8)
        {
          cat('Leave one patient out cross-validation for normal mass spectrometry data!\n')
          set.seed(123)
          cv.glmmod <- cv.glmnet(x = training_X, y = training_Y, 
                    alpha = 1, foldid = 1:length(training_Y), type.measure = "class",
                    family = "binomial", keep = T)
        } else {
          withProgress(message = 'Too few observations to run the lasso')
          cat('Too few observations to run the lasso!\n')
          
        }
        saveRDS(cv.glmmod, file = paste0(globalPath, "glmnet_cv_model"))
        cv.glmmod
      }
      else
      {
        cat('Building Lasso model for imaging data!\n')
        training.list <- getDataMat()
        training.table <- training.list[[1]]
        descriptor.mat <- training.list[[2]]
        
        training_Y <- descriptor.mat$Disease.State #training y must be 1 or 2
        training_X <- as.matrix(training.table)
        training_X <- scale(training_X, center = T)
        training_folds <- dense_rank(descriptor.mat$Patient)
        rm(training.table)
        rm(descriptor.mat)
        
        if (length(unique(training_folds)) >= 8)
        {
          cat('Using leave one patient out cross-validation!\n')
          set.seed(123)
          cv.glmmod <- cv.glmnet(x = training_X, y=training_Y, alpha=1,
                                 foldid = training_folds,
                                 family = "binomial", type.measure = "class",
                                 keep = T)
        } else {
          cat('Using 10-fold cross-validation!\n')
          set.seed(123)
          cv.glmmod <- cv.glmnet(x = training_X, y=training_Y, alpha=1, 
                                 family = "binomial", type.measure = "class",
                                 keep = T)
        }
        saveRDS(cv.glmmod, file = paste0(globalPath, "glmnet_cv_model"))
        cv.glmmod
      }
    })
  })
  modelPerformance <- function(){
    if (input$fileType == "Imaging")
    {
      cv.glmmod <- getLassoOutput()
      training.list <- getDataMat()
      training.table <- training.list[[1]]
      descriptor.mat <- training.list[[2]]
      tmp_coeffs <- coef(cv.glmmod, s = "lambda.min")
      selected_coefficients <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x,
                                          stringsAsFactors = F)
      perform.df <- data.frame(Description = c("Number of Unique Detected Metabolites",
                                               "Number of Observations (Pixels)",
                                               "Number of Selected LASSO Peaks",
                                               "Method of Model Selection",
                                               "Selected Model Tuning Parameter",
                                               "Accuracy",
                                               "True Positive Rate",
                                               "False Positive Rate",
                                               "True Negative Rate",
                                               "False Negative Rate"), Output = 0, stringsAsFactors = F)
      training_Y <- descriptor.mat$Disease.State
      imin <- which.min(cv.glmmod$cvm)
      yhat <- 1*(cv.glmmod$fit.pre>mean(training_Y))[,imin]
      #print(cv.glmmod$fit.pre[,imin])
      #write.csv(cv.glmmod$fit.pre[,imin], file = "cv.glmmod.fit.pre.csv", col.names = T, row.names = F)
      #print(yhat)
      #write.csv(yhat, file = "cv.glmmod.yhat.csv", col.names = T, row.names = F)
      
      perform.df[1,2] = ncol(training.table)
      perform.df[2,2] = nrow(training.table)
      perform.df[3,2] = nrow(selected_coefficients)-1
      perform.df[4,2] = ifelse(length(unique(descriptor.mat$Patient)) >= 8, "Leave one patient out CV", "10-fold CV")
      perform.df[5,2] = min(cv.glmmod$cvm)
      perform.df[6,2] = round(mean(training_Y == yhat), 3)
      perform.df[7,2] = round(length(which(training_Y == yhat & training_Y == 1))/length(which(training_Y == 1)), 3)
      perform.df[8,2] = round(length(which(training_Y != yhat & training_Y == 1))/length(which(training_Y == 1)), 3)
      perform.df[9,2] = round(length(which(training_Y == yhat & training_Y == 0))/length(which(training_Y == 0)), 3)
      perform.df[10,2] = round(length(which(training_Y != yhat & training_Y == 0))/length(which(training_Y == 0)), 3)
      print(perform.df)
      perform.df
    } else {
      cv.glmmod <- getLassoOutput()
      training_X <- getTrainingXMatrix()
      training_Y <- training_X$Disease.State
      training_X <- as.matrix(training_X[,-c(1:2)])
      #print('normal glmmod performance characteristics')
      print(training_X[1:5,1:5])
      tmp_coeffs <- coef(cv.glmmod, s = "lambda.min")
      selected_coefficients <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x,
                                          stringsAsFactors = F)
      perform.df <- data.frame(Description = c("Number of Unique Detected Metabolites",
                                               "Number of Observations",
                                               "Number of Selected LASSO Peaks",
                                               "Method of Model Selection",
                                               "Selected Model Tuning Parameter",
                                               "Accuracy",
                                               "True Positive Rate",
                                               "False Positive Rate",
                                               "True Negative Rate",
                                               "False Negative Rate"), Output = 0, stringsAsFactors = F)
      imin <- which.min(cv.glmmod$cvm)
      yhat <- 1*(cv.glmmod$fit.pre>mean(training_Y))[,imin]
      perform.df[1,2] = ncol(training_X)
      perform.df[2,2] = nrow(training_X)
      perform.df[3,2] = nrow(selected_coefficients)-1
      perform.df[4,2] = "Leave one patient out CV"
      perform.df[5,2] = min(cv.glmmod$cvm)
      perform.df[6,2] = round(mean(training_Y == yhat), 3)
      perform.df[7,2] = round(length(which(training_Y == yhat & training_Y == 1))/length(which(training_Y == 1)), 3)
      perform.df[8,2] = round(length(which(training_Y != yhat & training_Y == 1))/length(which(training_Y == 1)), 3)
      perform.df[9,2] = round(length(which(training_Y == yhat & training_Y == 0))/length(which(training_Y == 0)), 3)
      perform.df[10,2] = round(length(which(training_Y != yhat & training_Y == 0))/length(which(training_Y == 0)), 3)
      print(perform.df)
      perform.df
      

    }
  }
  output$modelPerformanceDT <- DT::renderDT({
    modelPerformance()
  }, rownames = F)
  selectedLassoPeaksPlotPrep <- function(){
    cv.glmmod <- getLassoOutput()
    #print()
    tmp_coeffs <- coef(cv.glmmod, s = "lambda.min")
    selected_coefficients <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x,
                                        stringsAsFactors = F)
    selected_coefficients <- selected_coefficients[-1,]
    print(selected_coefficients)
    plot(x = as.numeric(selected_coefficients$name), y = selected_coefficients[,2],
         type = 'h', col = 'blue', xlab = 'm/z', ylab = 'Coefficient Weight', cex.lab = 1.2, cex.axis = 1.3)
  }
  crossvalMSEPlotPrep <- function() {
    cv.glmmod <- getLassoOutput()
    plot(cv.glmmod, cex.axis = 1.2, cex.lab = 1.2)
  }
  output$crossvalMSEPlot <- renderPlot({
    crossvalMSEPlotPrep()
  })
  output$selectedLassoPeaks <- renderPlot({
    selectedLassoPeaksPlotPrep()
  })
  output$lassoPeaksTable <- DT::renderDT({
    cv.glmmod <- getLassoOutput()
    tmp_coeffs <- coef(cv.glmmod, s = "lambda.min")
    selected_coefficients <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
    selected_coefficients[order(abs(selected_coefficients$coefficient), decreasing = T),]
  }, rownames = F)
  
  finalSAMPlotNormal <- reactive({
    withProgress(message = 'Completing SAM Calculations', {
      #Takes the existing training data frame and re-formats it such that the columns
      #correspond to m/z peaks and the rows correspond to a unique observation
      #Manually checked this method for 2-3 peaks
      training_X <- getTrainingXMatrix()
      SAM_training_Y <- training_X$Disease.State + 1
      SAM_training_X = t(training_X[,-c(1:2)])
      print(SAM_training_X[1:5,1:5])
      results <- SAM(SAM_training_X, SAM_training_Y, genenames = rownames(SAM_training_X), resp.type = "Two class unpaired",
                     random.seed = 123)
      rm(SAM_training_X)
      rm(SAM_training_Y)
      incProgress(0.5, message = "Calculation Completed, Preparing Plot")
      cat('Completed SAM Calculation, Now Preparing Plot!\n')
      if(!is.null(results$siggenes.table$genes.up) && !is.null(results$siggenes.table$genes.lo))
      {
        print("case 1")
        results.table <- rbind(results$siggenes.table$genes.up, results$siggenes.table$genes.lo)
        print(results.table)
      }
      if(!is.null(results$siggenes.table$genes.lo) && is.null(results$siggenes.table$genes.up))
      {
        print("case 2")
        results.table <- results$siggenes.table$genes.lo
        print(results.table)
      }
      if(is.null(results$siggenes.table$genes.lo) && !is.null(results$siggenes.table$genes.up))
      {
        print("case 3")
        results.table <- results$siggenes.table$genes.up
        print(results.table)
      }
      incProgress(amount = 0.2)
      if (nrow(results.table) > 0) {  
        final.results <- results.table[as.numeric(results.table[,7]) < 5 & as.numeric(results.table[,6])!=Inf,]
        if (!is.null(dim(final.results)[1]) && dim(final.results)[1] > 0)
        {
          final.results <- as.matrix(final.results[order(as.numeric(final.results[,1])),])
          SAM.plot <- data.frame(Peaks = as.numeric(final.results[,1]), Peak.Type = "Up regulated", Fold.Change=as.numeric(final.results[,6]), stringsAsFactors = F)
          print(head(final.results))
        } else {
          final.results <- as.data.frame(final.results)
          SAM.plot <- data.frame(Peaks = final.results[1,], Peak.Type = "Up regulated", Fold.Change=as.numeric(final.results[,6]), stringsAsFactors = F)
        }
        SAM.plot[,2] = ifelse(as.numeric(final.results[,3]) < 0, "Down regulated", "Up regulated")
        list(SAM.plot, final.results)
      } else {
        df<- data.frame()
        df
      }
    })
  })
  finalSAMPlotImaging <- reactive({
    withProgress(message = 'Preparing SAM Calculations', {
      training.list <- getDataMat()
      training.table <- training.list[[1]]
      descriptor.mat <- training.list[[2]]
      SAM_training_Y <- descriptor.mat$Disease.State+1 #training y must be 1 or 2
      SAM_training_X <- t(training.table) #remove the auxillary-non peak variables
      rm(training.table) 
      rm(descriptor.mat)
      cat('Running SAM Calculation for Imaging Mode Data!\n')
      print(SAM_training_X[1:5,1:5])
      incProgress(0.2, message = "Running SAM Calculation")
      sink(paste0(globalPath,"SAM_progress.txt"))
      results <- SAM(SAM_training_X, SAM_training_Y, genenames = rownames(SAM_training_X), resp.type = "Two class unpaired",
                     random.seed = 123)
      saveRDS(results, file = paste0(globalPath, "SAM_results"))
      sink(NULL)
      closeAllConnections()
      results.table <- rbind(results$siggenes.table$genes.up, results$siggenes.table$genes.lo)
      final.results <-results.table[as.numeric(results.table[,7]) < 5 & as.numeric(results.table[,6])!=Inf,]
      final.results <- final.results[order(as.numeric(final.results[,1])),]
      cat('SAM Calculations Completed for Imaging Mode Data!\n')
      
      if (nrow(final.results) > 0) {  
        SAM.plot <- data.frame(Peaks = as.numeric(final.results[,1]), Peak.Change = "Small.Metabolite", Fold.Change=as.numeric(final.results[,6]), stringsAsFactors = F)
        printable.results = data.frame(Peaks = as.numeric(final.results[,1]), Score = as.numeric(final.results[,3]),
                                       Numerator= as.numeric(final.results[,4]), Denominator=as.numeric(final.results[,5]),
                                       Fold.Change = as.numeric(final.results[,6]), Q.Value = as.numeric(final.results[,7]), stringsAsFactors = F)
        printable.results$Change = ifelse(printable.results$Score < 0, "Down regulated", "Up regulated")
        #printable.results$Fold.Change = ifelse(printable.results$Score < 0, -1*printable.results$Fold.Change, printable.results$Fold.Change)
        SAM.plot[which(printable.results$Score < 0),2] = "Down regulated"
        SAM.plot[which(printable.results$Score > 0),2] = "Up regulated"
        list(SAM.plot, printable.results)
      } else {
        df <- data.frame()
        df
      }
    })
  })
  samPlotPrep <- function() {
    #training.list <- getDataMat()
    if (input$fileType == "Imaging")
    {
      plots <- finalSAMPlotImaging()
      final.plot <- plots[[1]]
      rm(plots)
      if (nrow(final.plot) > 0 && !is.na(final.plot[,1])) {
        #print("final plot imaging")
        #final.plot$Fold.Change = ifelse(final.plot$Fold.Change < 1.0, -1.0*log2(1/final.plot$Fold.Change), log2(final.plot$Fold.Change))
        plot(x = final.plot$Peaks, y = log2(final.plot$Fold.Change), 
             col = 'blue', type = 'h', xlab = "m/z", ylab = "Fold Change (Log2)", 
             main = "Peaks with statistically significant differences in relative peak intensity (q < 0.05)",
             cex.lab = 1.3, cex.axis = 1.3)
        abline(h = 0, col = 'black')
      } else {
        plot(x = 50:1000, y = rep(0, length(50:1000)), main = "No statistically significant differences detected",
             type = 'h', col = 'blue', xlab = "m/z", ylab = "Fold Change", cex.lab = 1.3, cex.axis = 1.3)
      }
    } else {
     # print("final plot normal")
      final.plot <- finalSAMPlotNormal()[[1]]
      if (nrow(final.plot) > 0 && !is.na(final.plot[,1])) {
        #final.plot$Fold.Change = ifelse(final.plot$Fold.Change < 1.0, -1.0*log2(1/final.plot$Fold.Change), log2(final.plot$Fold.Change))
        plot(x = final.plot$Peaks, y = log2(final.plot$Fold.Change), 
             col = 'blue', type = 'h', xlab = "m/z", ylab = "Fold Change (Log2)", 
             main = "Peaks with statistically significant differences in relative peak intensity (q < 0.05)",
             cex.lab = 1.3, cex.axis = 1.3)
        abline(h = 0, col = 'black')
      } else {
        #print("nothing panned out!")
        plot(x = 50:1000, y = rep(0, length(50:1000)), main = "No statistically significant differences detected",
             type = 'h', col = 'blue', xlab = "m/z", ylab = "Fold Change", cex.lab = 1.3, cex.axis = 1.3)
      }
    }
  }
  output$clusterPlot <- renderPlot({
      samPlotPrep()
  })
  output$interestingClusters = DT::renderDT({
    if (input$fileType == "Imaging")
    {
      all.plots <- finalSAMPlotImaging()
      final.plot <- all.plots[[1]]
    }else if (input$fileType == "Normal") {
      all.plots <- finalSAMPlotNormal()
      final.plot <- all.plots[[1]]
    } else {
      return(NULL)
    }
    final.plot <- final.plot[order(abs(final.plot$Fold.Change), decreasing = T),]
    colnames(final.plot) <- c("Peak", "Direction of Change", "Fold Change")
    print(final.plot)
    #final.plot[,c(1,3)]
  }, rownames=F)
  getTrainingXMatrix <- reactive({
    training.list <- getDataMat()
    training.df <- as.data.table(training.list[[1]])
    if (input$fileType == "Normal")
    {
      training.df <- dcast.data.table(training.df, formula = File.Name + Disease.State ~ Centroid, 
                                     fun.aggregate = max, value.var = "Intensity", fill = 0)
      if (input$threshold != -1 && input$threshold > 0 && input$threshold < 1)
      {
        cat('Removing sparse metabolites (normal mode) \n')
        descriptor_info = which(colnames(training.df) %in% c("File.Name", "Disease.State"))
        fraction_nonzero = colSums(training.df != 0)/nrow(training.df)
        selected_ind = which(fraction_nonzero >= input$threshold)
        print(dim(training.df))
        training.df <- training.df[,unique(c(descriptor_info, selected_ind)), with = F]
        print(dim(training.df))
      }
      sink(NULL)
      saveRDS(training.df, file = paste0(globalPath, "normal_data_training_df"))
      training.df
    } 
  })
  
  output$dropDownPeaks <- renderUI({
    withProgress(message = 'Setting Up Display', {
      if (input$fileType == "Normal") #in normal file format, disease.state is the first column which is trainingY
      {
        training_X <- getTrainingXMatrix()
        training.list <- getDataMat()
        inputType <- training.list[[3]]
        unique.pks = as.numeric(names(training_X[,-c(1:2)]))
      } else { #in imaging mode
        training.list <- getDataMat()
        training.table <- training.list[[1]]
        unique.pks = as.numeric(names(training.table[,-c(1:2)]))
      }
      incProgress(amount = 0.5, message = "Collected Peaks, Now Setting up Boxplot")
      selectizeInput(inputId = "peakChoice", label = "Selected Peak", unique.pks, multiple = F, 
                  selected = unique.pks[2], options = list(maxItems = length(unique.pks)))
    })
  })
  averageSpectraPlotPrep <- function(){
    if (input$fileType == "Normal") {
      training_X <- getTrainingXMatrix()
      peaks = colnames(training_X)[-c(1:2)]
      healthy_state = min(training_X$Disease.State)
      healthy_ind = which(training_X$Disease.State == healthy_state)
      plot(x = as.numeric(peaks), y = colMeans(training_X[healthy_ind,-c(1:2)]), xlab = "m/z", ylab = "Relative Intensity",
           type = 'h', col = 'green', cex.axis = 1.3, cex.lab = 1.3)
      lines(x = as.numeric(peaks), y = colMeans(training_X[-healthy_ind,-c(1:2)]), xlab = "m/z", ylab = "Relative Intensity",
            type = 'h', col = 'red', cex.axis = 1.3, cex.lab = 1.3)
      legend("topright", legend=c(input$labelSet1, input$labelSet2),
             col=c("green", "red"), lty=1, text.font = c(15))
    } else {
      training.list = getDataMat()
      training.table = training.list[[1]]
      descriptor.table = training.list[[2]]
      print(head(descriptor.table))
      print(training.table[1:5,1:5])
      peaks = colnames(training.table)
      healthy_state = min(descriptor.table$Disease.State)
      healthy_ind = which(descriptor.table$Disease.State == healthy_state)
      plot(x = as.numeric(peaks), y = colMeans(training.table[healthy_ind,]), xlab = "m/z", ylab = "Relative Intensity",
           type = 'h', col = 'green', cex.axis = 1.3, cex.lab = 1.3)
      lines(x = as.numeric(peaks), y = colMeans(training.table[-healthy_ind,]), xlab = "m/z", ylab = "Relative Intensity",
            type = 'h', col = 'red', cex.axis = 1.3, cex.lab = 1.3)
      legend("topright", legend=c(input$labelSet1, input$labelSet2),
             col=c("green", "red"), lty=1, text.font = c(15))
    }
  }
  output$spatialHealthyMetaboliteDist <- renderPlotly({
    if(input$fileType == "Normal")
    {
      df <- data.frame()
      ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
    } else {
      selectedPeak <- input$peakChoice
      training.list = getDataMat()
      training.table = training.list[[1]]
      descriptor.table = training.list[[2]]
      #print(descriptor.table)
      descriptor.table$Disease.State = descriptor.table$Disease.State + 1
      selectedPeakIndex = which(colnames(training.table) %in% selectedPeak)
      final.table = bind_cols(descriptor.table, training.table[,selectedPeakIndex, with = F])
      peak.name = trunc(as.numeric(colnames(final.table)[6])*1000)/1000
      colnames(final.table)[6] <- "Relative.Intensity"
      #Can't assign 6 names
      healthy_state = min(final.table$Disease.State)
      final.table = final.table[which(final.table$Disease.State == healthy_state),]
      processed.table = as.data.frame(final.table %>% group_by(X, Y) %>% summarise(Relative.Intensity = mean(Relative.Intensity)))
      print(dim(processed.table))
      ggplot(data = processed.table, mapping = aes(x = X, y = Y)) + 
        geom_tile(aes(fill = Relative.Intensity), color = 'black') + 
        theme_void() + ggtitle(label = input$labelSet1) +
        scale_fill_gradient2(low = "white", high = "black") + xlab(label = "X") + ylab(label = "Y") 
    }
  })
  output$spatialDiseaseMetaboliteDist <- renderPlotly({
    if(input$fileType == "Normal")
    {
      df <- data.frame()
      ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
    } else {
      selectedPeak <- input$peakChoice
      training.list = getDataMat()
      training.table = training.list[[1]]
      descriptor.table = training.list[[2]]
      descriptor.table$Disease.State = descriptor.table$Disease.State + 1
      selectedPeakIndex = which(colnames(training.table) %in% selectedPeak)
      final.table = bind_cols(descriptor.table, training.table[,selectedPeakIndex, with = F])
      peak.name = trunc(as.numeric(colnames(final.table)[6])*1000)/1000
      colnames(final.table)[6] <- "Relative.Intensity"
      #Can't assign 6 names
      healthy_state = min(final.table$Disease.State)
      final.table = final.table[which(final.table$Disease.State != healthy_state),]
      processed.table = as.data.frame(final.table %>% group_by(X, Y) %>% summarise(Relative.Intensity = mean(Relative.Intensity)))
      #print(processed.table)
      print(dim(processed.table))
      ggplot(data = processed.table, mapping = aes(x = X, y = Y, fill = Relative.Intensity)) +
        geom_tile(color = 'black') +
        xlab(label = "X") + ylab(label = "Y") + scale_fill_gradient2(low = "white", high = "black") + 
        theme_void() + ggtitle(label = input$labelSet2)
    }
  })
  
  output$averageSpectraPlot <- renderPlot({
    withProgress(message = 'Preparing Avg. Spectra Plot', {
      averageSpectraPlotPrep()
    })
  })
  peakVarBoxPlotPrep <- reactive({
    withProgress(message = 'Preparing Boxplot', {
      selectedPeak <- input$peakChoice
      if (input$fileType == "Normal") #in normal file format, disease.state is the first column which is trainingY
      {
        training_X <- getTrainingXMatrix()
        fileAdata <- data1()
        fileBdata <- data2()
        fileNames <- c(names(fileAdata), names(fileBdata))
        rm(fileAdata)
        rm(fileBdata)
        #we need the file names to access the meta-data that accompanies every relative intensity
        training_X.final = training_X[,c(1:2, which(names(training_X) %in% selectedPeak)), with = F]
        training_X.final[,2] = training_X.final[,2]+1
        training_X.final[, "File.Name" := c(fileNames)] 
        peak.name = trunc(as.numeric(colnames(training_X.final)[3])*1000)/1000
        colnames(training_X.final)[3] <- "Relative.Intensity"
        healthy_state = min(training_X.final$Disease.State)
        training_X.final$Disease.State = as.factor(ifelse(training_X.final$Disease.State == healthy_state, input$labelSet1, input$labelSet2))
        training_X.final = as.data.frame(training_X.final)
        p <- ggboxplot(data = training_X.final, x = "Disease.State", y = "Relative.Intensity",
                       palette = "aaas", add = "jitter", fill="Disease.State", order = c(input$labelSet1, input$labelSet2),
                       font.tickslab = c(18)) 
        p <- p + grids(linetype = "dashed") + xlab("") + ylab("Relative Intensity") + theme(legend.position = "none")
      } else { #in imaging mode, disease.state is the 5th column which is trainingY
        training.list = getDataMat()
        training.table = training.list[[1]]
        descriptor.table = training.list[[2]]
        descriptor.table$Disease.State = descriptor.table$Disease.State + 1
        selectedPeakIndex = which(colnames(training.table) %in% selectedPeak)
        final.table = bind_cols(descriptor.table, training.table[,selectedPeakIndex, with = F])
        peak.name = trunc(as.numeric(colnames(final.table)[6])*1000)/1000
        colnames(final.table)[6] <- "Relative.Intensity"
        #Can't assign 6 names
        healthy_state = min(final.table$Disease.State)
        final.table$Disease.State = as.factor(ifelse(final.table$Disease.State == healthy_state, input$labelSet1, input$labelSet2))
        final.table = as.data.frame(final.table)
        p <- ggboxplot(data = final.table, x = "Disease.State", y = "Relative.Intensity",
                       palette = "aaas", add = "jitter", fill="Disease.State", order = c(input$labelSet1, input$labelSet2),
                       font.tickslab = c(18)) 
        #outlier setting: geom_boxplot(outlier.colour="red", outlier.shape=8,
        #             outlier.size=4)
        p<- p + xlab("")  + 
          grids(linetype = "dashed") + ylab("Relative Intensity") + theme(legend.position = "none")
      }
      ggpar(p, font.y = c(14, "bold"))
      
    })
  })
  
  output$peakvarBoxPlot <- renderPlot({
    print(peakVarBoxPlotPrep())
  })
  
  output$download1 <- downloadHandler(
    filename = function() {
      paste("training_matrix.csv", sep = "")
    },
    content = function(file) {
      if (input$fileType == "Normal")
      {
        training_X <- getTrainingXMatrix()
        withProgress(message = 'Writing Training Matrix', {  
          setProgress(value = 0.5)
          write.csv(training_X, file, row.names = FALSE)
        })
      } else {
        withProgress(message = 'Constructing Training Matrix', {
          training.list <- getDataMat()
          training.table <- training.list[[1]]
          descriptor.mat <- training.list[[2]]
          final.out<- bind_cols(descriptor.mat, training.table)
          incProgress(amount = 0.5, message = "Training Matrix Constructed, Writing Out Matrix")
          fwrite(final.out, file, append = F, row.names = F, col.names = T)
          #write.csv(final.out, file, row.names = FALSE)
        })
      }
    }
  )
  output$downloadMetabolitePlot <- downloadHandler(
    filename = function() { 
      paste(input$peakChoice, '_metabolite_distribution.png', sep='') 
    },
    content = function(file) 
    {
      final.plot <- peakVarBoxPlotPrep()
      ggsave(file, final.plot, dpi=400, height=4, width=5)
    }
  )
  output$downloadlassoPeaks <- downloadHandler(
   filename = paste0('selected_lasso_peaks.pdf'),
   content = function(file){
     pdf(file = file, height = 6, width = 11)
     selectedLassoPeaksPlotPrep()
     dev.off()
   }
  )
  output$downloadAvgSpectraPlot <- downloadHandler(
    filename = paste0('average_spectra_comparison.pdf'),
    content = function(file){
      pdf(file=file, height = 5, width = 7.5)
      averageSpectraPlotPrep()
      dev.off()
    }
  )
  output$downloadcrossvalPlot <- downloadHandler(
    filename = paste0('Lasso_CV_Performance.pdf'),
    content = function(file){
      pdf(file=file, height = 6, width = 11)
      crossvalMSEPlotPrep()
      dev.off()
      
    }
  )
  output$downloadSAMPlot <- downloadHandler(
    filename = paste0('SAM_Peak_Comparisons.pdf'),
    content = function(file) {
      pdf(file = file, height = 6, width = 11)
      samPlotPrep()
      dev.off()
    }
  )
  
  output$downloadModel <- downloadHandler(
    filename= function(){
      paste('selected_lasso_coefficients.csv', sep = "")
    },
    content = function(file)
    {
      cv.glmmod <- getLassoOutput()
      tmp_coeffs <- coef(cv.glmmod, s = "lambda.min")
      selected_coefficients <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x,
                                          stringsAsFactors = F)
      write.csv(selected_coefficients, file, row.names = F, col.names = T)
    }
  )
  output$downloadPerformance <- downloadHandler(
    filename = function(){
      paste('lasso_selected_cv_model_performance.csv', sep = "")
    },
    content = function(file)
    {
      performance.df <- modelPerformance()
      write.csv(performance.df, file, col.names = T, row.names = F)
    }
  )
  output$downloadSAMOutput <- downloadHandler(
    filename = function(){
      paste("sam_output.csv", sep = "")
    },
    content = function(file) 
    {
      if (input$fileType == "Normal")  
      {
        all.plots <- finalSAMPlotNormal()
        final.plot <- all.plots[[2]]
        colnames(final.plot)[1:2]<-c("Peaks", "Peak ID")
        final.plot = cbind(final.plot, Direction = ifelse(as.numeric(final.plot[,3]) < 0, "Down regulated", "Up regulated"))
        print(final.plot)
      } else if(input$fileType == "Imaging") {
        all.plots <- finalSAMPlotImaging()
        final.plot <- all.plots[[2]]
      } else {
        return(NULL)
      }
      write.csv(final.plot, file, row.names = F, col.names = T)
    }
  )
})
