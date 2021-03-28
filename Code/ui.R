#Vishnu Shankar
#Mass Spectrometry Online Tool
#July 10, 2020
library(shiny)
library(shinydashboard)
library(ggplot2)
library(DT)
library(plotly)

dashboardPage(title = "MassExplorer",
              skin = "black",
              header <- dashboardHeader(
                title = tags$b("MassExplorer"),
                titleWidth = 270
              ),
              sidebar <- dashboardSidebar(
                sidebarMenu(
                  id="tabs",
                  menuItem(tags$b("Input Data"), tabName = "setup", icon = icon("sitemap"), selected = T),
                  menuItem(tags$b("Pairwise Comparisons"),tabName = "crudeOutput", icon = icon("chevron-circle-right")),
                  menuItem(tags$b("Get Training Matrix"),tabName = "trainingMat", icon = icon("chevron-circle-right")),
                  menuItem(tags$b("Find Specific Differences"),tabName = "specificOutput", icon = icon("chevron-circle-right")),
                  menuItem(tags$b("Build Model"),tabName = "buildModel", icon = icon("chevron-circle-right")),
                  menuItem(tags$b("Understand Peak Variability"),tabName = "exploreData", icon = icon("chevron-circle-right")),
                  #menuItem(tags$b("Evaluate Model"),tabName = "evalModel", icon = icon("chevron-circle-right")),
                  menuItem(tags$b("About MassExplore"),tabName = "about", icon = icon("info"))
                )    
              ),
              body <- dashboardBody(
                tabItems(
                  tabItem("setup",
                          fluidRow(box(width = 12, title = tags$b("Welcome to MassExplorer!"), solidHeader = T, status = 'primary', collapsible = F, 
                                       tags$h4("This application helps identify differences in the abundances of detected molecular species between two groups of mass spectrometry data. It is uniquely suited
                                        for processing Desorption Electrospray Ionization Mass Spectrometry Imaging (DESI-MSI) data."),
                                       tags$br(),
                                       tags$h4("It consists of the following modules:"),
                                       tags$h4("(1) Pairwise comparisons: Visualize mass spectra and compare individual spectra between groups"),
                                       #tags$br(),
                                       tags$h4("(2) Get Training Matrix: Download average spectra and processed data for further analysis"),
                                       #tags$br(),
                                       tags$h4("(3) Find Specific Differences: Identify statistically significant differences in peak abundances between groups"),
                                       #tags$br(),
                                       tags$h4("(4) Build Model: Train a predictive model to distinguish groups using the LASSO via cross-validation"),
                                       #tags$br(),
                                       tags$h4("(5) Understand Peak Variability: Study metabolite variability and spatial distribution"),
                                       #tags$br(),
                                       #tags$h4("(6) Evaluate Model: Evaluate a pre-trained model on inputted mass spectrometry data"),
                                       tags$br(),
                                       tags$h6("Please contact the developer Vishnu Shankar (vishnus1@cs.stanford.edu), if you have any feedback or questions.", style ="color:red")
                                   )),
                          fluidRow(
                            box(width = 6,title = tags$b('How should the data be processed?'), solidHeader = T, status = 'primary', collapsible = T,
                                radioButtons("ticNormalization", label="Normalize by Total Ion Current?", choices=c("Yes", "No"), 
                                             selected = "Yes", inline = T),
                                radioButtons("scaleNormalization", label="Scale Maximum Ion Abundances to 1.0?", choices=c("Yes", "No"), 
                                             selected = "Yes", inline = T),
                                radioButtons("peakElbow", label = "Exclude peak elbows during data processing?", choices=c("Yes","No"), selected="Yes", inline = T),
                                numericInput(inputId = "controlPeak", label = "Specify the internal standard m/z peak to two decimal places (e.g. 514.28). The default (-1) corresponds to no internal standard.", value = -1),
                                numericInput(inputId = "error", label = "Specify the tolerance (m/z units) for peak shifting (how much can a peak shift between samples and be considered the same?)", value = 0.05),
                                numericInput(inputId = "threshold", label = "Peaks included in the analysis should be found in at least what fraction of samples? The default (-1) includes all detected metabolites in the analysis. Values should be between (0,1).", value = -1)
                            ),
                            box(width = 6, title = tags$b("Label Sample Sets"), solidHeader = T, status = 'primary', collapsible = T, 
                                HTML("<br>"),
                                textInput("labelSet1", tags$h4(tags$b("Label the first group of samples (e.g. Healthy)")), "Healthy"), 
                                HTML("<br><br><br><br><br>"),
                                textInput("labelSet2", tags$h4(tags$b("Label the second group of samples (e.g. Disease)")), "Disease"),
                                HTML("<br><br><br><br><br><br><br>")
                            )),
                          fluidRow(
                            box(width = 6, title = tags$b('Select File Format'), solidHeader = T, status = 'danger', 
                                helpText("For 'Normal' file formats, multiple spectra can be inputted as separate files with the format below."),
                                helpText("For 'Imaging' file formats, multiple spectra can be inputted in the same and separate files. Imaging filenames should have a unique number for sample identification
                                         (e.g. 1A.csv/Patient_1.csv, 2A/Patient_2.csv). Please ensure that no two files from distinct samples or patients in the same group have the same identification number."),
                                fluidRow(column(3, tags$b("Select file format:"), radioButtons("fileType", " ",choices = c("Normal", "Imaging"), selected = "Normal")
                                ),
                                column(3, tags$b("Sample File:"),
                                       conditionalPanel(condition = "input.fileType == 'Normal'", tableOutput('tbl')),
                                       conditionalPanel(condition = "input.fileType == 'Imaging'", tableOutput('tbl2')))
                                )
                            ),
                            box(width = 6,title = tags$b('Input Dataset'), solidHeader = T, status = 'danger', collapsible = F,
                                helpText("Please upload both files for analysis in csv format. Specify the appropriate format in 'Format Input Data' panel."),
                                fileInput("file1","Upload the first group of samples (e.g. healthy samples)", multiple=TRUE),
                                fileInput("file2","Upload the second group of samples (e.g. diseased samples)", multiple=TRUE)
                            )
                          )
                  ),
                  tabItem("crudeOutput",
                          fluidRow(box(width = 12, title = tags$b("Instructions"), solidHeader = T, status = 'primary',
                                       helpText("Use the slider below to adjust the m/z range for both selected mass spectra. The plot and table below
                                    reflect relative intensity, based on the data processing settings selected in the Input module."),
                                       sliderInput("mz_range", "Select m/z Range of Interest:", min = 0, max = 2000, value = c(0,2000), width = 1200)
                          )),
                          fluidRow(
                            column(6, uiOutput("choicesFileA")),
                            column(6, uiOutput("choicesFileB"))
                          ),
                          fluidRow(
                            tags$style(type="text/css",
                                       ".shiny-output-error { visibility: hidden; }",
                                       ".shiny-output-error:before { visibility: hidden; }"
                            ),
                            column(6, plotlyOutput("plotSetA")),
                            column(6, plotlyOutput("plotSetB"))
                          ),
                          fluidRow(
                            column(6, 
                                   box(width = 20, title = tags$b('Selected m/z Spectra for Set A'), solidHeader = T, status = 'primary', collapsible = T,
                                       DTOutput('mytableA'))),
                            column(6, 
                                   box(width = 20, title = tags$b('Selected m/z Spectra for Set B'), solidHeader = T, status = 'primary', collapsible = T,
                                       DTOutput('mytableB')))
                          )
                  ),
                  tabItem("specificOutput",
                          fluidRow(box(width = 12, title = tags$b("What does this module do?"), solidHeader = T, status = 'primary',
                                       helpText("This module calculates the statistically significant differences in relative abundance of detected molecular species between groups.
                                        The calculation is done using significance analysis of microarrays (SAM) (PNAS 2001 98: 5116-5121), which uses repeated permutations
                                        to account for the false discovery rate and determine if the changes in relative abundances of analytes is significantly related to the outcome.")
                                       )
                          ),
                          fluidRow(box(width = 12, title = tags$b('Statistically Significant Differences in Peaks'), 
                                       solidHeader = T, status = 'primary', 
                                       fluidRow(column(12, plotOutput("clusterPlot"))), 
                                       downloadButton('downloadSAMPlot', 'Download Plot'))
                                   ),
                          fluidRow(column(12, box(width = 20, title = tags$b('How do you interpret the output below?'), solidHeader = T,
                                                  status = 'primary',
                                                helpText("If a particular metabolite m/z 89.02 is 'Up Regulated' with a Fold Change of '4.1', it means in 
                                                the second set of samples the mean relative intensity of m/z 89.02 is 4.1 times higher than the mean in the first set of samples.
                                                If m/z 89.02 is 'Down Regulated' with Fold Change of 0.8, it means in the second set of samples the mean relative intensity is 0.8 times
                                                the mean relative intensity of the first set of samples. Note the fold change table of significant peaks shown below is not log transformed, unlike the plot above.")))),
                          fluidRow(column(12, box(width = 20, title = tags$b('Table of Statistically Significant Peaks'), solidHeader = T, 
                                                  DTOutput('interestingClusters'), status = 'primary', collapsible = T),
                                          downloadButton("downloadSAMOutput", "Download SAM Output"))
                          )
                  ), 
                  tabItem("buildModel",
                          fluidRow(box(width = 12, title = tags$b("What does this module do?"), solidHeader = T, status = 'primary',
                                       helpText("This module builds a predictive binary logistic regression pixel-based classifier using the LASSO (J Roy Stat Soc Ser B 1996; 58:267–88.). 
                                                  As the LASSO is known to return a sparse model, it can help reduce model overfitting and focus on the most informative
                                                molecular species that can help distinguish both groups.
                                                The model is trained on the relative abundances of detected metabolites and is selected via cross-validation. This module requires
                                                8 or more mass spectra to fit the model. If there are more than 8 patients, a leave-one patient out cross-validation scheme is used (See https://doi.org/10.1002/ijc.32843 for more details). 
                                                Otherwise, this module arbitrarily splits the data into 10 folds and optimizes the model accordingly. "))
                         ),
                         fluidRow(box(width = 12, title = tags$b('Cross-validation Performance: Selection of Tuning Parameter'),
                                      solidHeader = T, status = 'primary',
                                      helpText('The plot below shows the cross-validation classification performance as a function
                                               of the tuning parameter (lambda), selected via cross-validation. The vertical
                                               dotted line corresponds to the selected model. The numbers on top of the plot
                                               indicate the number of non-zero peaks that correspond to the model 
                                               with the particular tuning parameter and classification error.'),
                                      fluidRow(column(12, plotOutput('crossvalMSEPlot'))),
                                      downloadButton('downloadcrossvalPlot', 'Download Cross Validation Plot')
                                      )),
                         
                          fluidRow(box(width = 12, title = tags$b('Selected LASSO Peaks'), 
                                       solidHeader = T, status = 'primary', 
                                       fluidRow(column(12, plotOutput('selectedLassoPeaks'))),
                                       downloadButton('downloadlassoPeaks', 'Download Plot'),
                                       downloadButton('downloadModel', 'Download Model'),
                                       DTOutput('lassoPeaksTable')
                                      )),
                         fluidRow(box(width = 12, title = tags$b('Model Performance'),
                                      solidHeader = T, status = 'primary',
                                      fluidRow(column(12, DTOutput('modelPerformanceDT'))),
                                      downloadButton('downloadPerformance', 'Download Performance Statistics')
                                      )
                                  )
                  ),
                  tabItem("exploreData", 
                        fluidRow(box(width = 12, title = tags$b('Understand Peak Variability'), solidHeader = T, status = 'primary', 
                                       fluidRow(column(12, uiOutput("dropDownPeaks"))),
                                       fluidRow(column(12, plotOutput("peakvarBoxPlot"))),
                                       downloadButton('downloadMetabolitePlot', 'Download Plot'),
                                       fluidRow(column(6, plotlyOutput('spatialHealthyMetaboliteDist')),
                                                column(6, plotlyOutput('spatialDiseaseMetaboliteDist')))
                          ))
                          #box(width = 12, title = tags$b('Visualize Spatial Distribution'), solidHeader = T, status = 'primary', 
                          #    fluidRow(column(6, plotlyOutput('spatialHealthyMetaboliteDist')),
                          #             column(6, plotlyOutput('spatialDiseaseMetaboliteDist')))
                              #downloadButton('downloadSpatialPlot', 'Download Plot')
                          #))
                  ), 
                  tabItem("trainingMat",
                          fluidRow(box(width = 12, title = tags$b('Compare Average Spectra'), solidHeader = T, status = 'primary',
                                       fluidRow(column(12, plotOutput('averageSpectraPlot'))),
                                       downloadButton('downloadAvgSpectraPlot', 'Download Plot'))),
                          fluidRow(box(width = 12, title = tags$b('Download Training Matrix'), solidHeader = T, status = 'primary',
                                       helpText("The processed data (training matrix) is formatted as follows: The matrix has n x p 
                                    dimensions, with n = number of unique spectra and
                                    p = the processed detected analytes. The matrix is populated by the intensities
                                    for each patient with 0 indicating that the analyte was not detected for the particular
                                    spectra. The column corresponding to 'Disease.State' indicates what class each observation
                                    corresponds to. In the case of 'Normal' file type, an observation (row) corresponds to each patient/file. 
                                    In the case of 'Imaging' file type, we consider each scan, X, Y from each patient as a new row"),
                                       downloadButton('download1',"Download Training Matrix")))
                  ),
                  tabItem("about",
                          fluidRow(box(width = 12, title = tags$b('About MassExplore'), solidHeader = T, status = 'primary',
                                       helpText("This application was developed by Vishnu Shankar (vishnus1@stanford.edu), with the guidance of 
                                    Prof. Richard N. Zare and Prof. Robert Tibshirani at Stanford University, along with important feedback from members of the Zare lab. Please send an email to 
                                    vishnus1@stanford.edu if you have any comments, questions, or issues running the program."))) 
                  )
                )
              )
)
              