library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyWidgets)
library(shinyjs)
library(shinythemes)
library(ggplot2)
library(knitr)
library(rmarkdown)
library(randomForest)
library(caret)
library(rmarkdown)

load("ABMR_RF_model_3000.RData")
load("ABMR_Outcome.RData")
load("ordered_gene_data.RData")
load("TCMR_RF_model_1000.RData")
load("TCMR_Outcome.RData")
load("ordered_gene_data_tcmr.RData")

# Define UI for application
ui <- navbarPage(theme = shinytheme("yeti"),
                 title = "GenePath",
                 tabPanel(title = "Predictor",
                          fluidPage(
                            useShinyjs(),
                            useShinydashboard(),
                            
                            # Application title
                            titlePanel("GenePath"),
                            
                            # Sidebar with a slider input for number of bins 
                            sidebarLayout(
                              sidebarPanel(
                                style = "height:170px;",
                                fileInput("file", "Please Upload Your File.", accept = ".csv"),
                                
                                actionButton("Result", "Result")
                              ),
                              mainPanel(
                                valueBoxOutput("Result_A", width = 6),
                                valueBoxOutput("Result_T", width = 6)
                              )
                            ),
                            
                            sidebarLayout(
                              sidebarPanel(
                                checkboxInput("child", "Child Patient", FALSE),
                                downloadButton("download")
                              ),
                              
                              mainPanel(
                                htmlOutput("text")
                              )
                            ) 
                          )),
                 
                 tabPanel(title = "Top Genes",
                          fluidRow(
                            box(
                              title = "ABMR",
                              width = 12, 
                              status = "primary",
                              fluidRow(
                                splitLayout(cellWidths = c("20%", "20%", "20%", "20%", "20%"), 
                                            plotOutput("boxa1"), 
                                            plotOutput("boxa2"), 
                                            plotOutput("boxa3"),
                                            plotOutput("boxa4"),
                                            plotOutput("boxa5"))
                              )
                            ), 
                            
                            box(
                              title = "TCMR",
                              width = 12, 
                              status = "primary",
                              fluidRow(
                                splitLayout(cellWidths = c("20%", "20%", "20%", "20%", "20%"), 
                                            plotOutput("boxt1"), 
                                            plotOutput("boxt2"), 
                                            plotOutput("boxt3"),
                                            plotOutput("boxt4"),
                                            plotOutput("boxt5"))
                              )
                            )
                          ))
)

# Define server logic
server <- function(input, output, session) {
  
  csv_data <- reactive({
    inFile <- input$file
    if (is.null(inFile)) return(NULL)
    read.csv(inFile$datapath)
  })
  
  Result <- reactiveValues(abmr = -1, tcmr = -1, R = "None", RT = "None", M = "None", G = "None", plota = NULL, plott = NULL)
  
  observeEvent(input$Result, {
    data <- csv_data()  # Get the uploaded data
    
    # Get the variable names from the loaded model
    top_ABMR <- ordered_gene_data[1:5,]
    top_TCMR <- ordered_gene_data_tcmr[1:5,]
    
    var_names_abmr <- rownames(varImp(final_rf_model_3000))
    var_names_tcmr <- rownames(varImp(final_rf_model_tcmr))
    
    common_genes_abmr <- intersect(data$V1, var_names_abmr)
    common_genes_tcmr <- intersect(data$V1, var_names_tcmr)
    
    topgene_abmr <- intersect(data$V1, rownames(top_ABMR))
    topgene_tcmr <- intersect(data$V1, rownames(top_TCMR))
    
    rownames(data) = data$V1
    subset_abmr <- data[common_genes_abmr, ]
    subset_tcmr <- data[common_genes_tcmr, ]
    
    top_csv_a <- data[topgene_abmr, ]
    top_csv_t <- data[topgene_tcmr, ]
    # print(top_csv_a)
    # print(rownames(top_ABMR))
    abmr_rn <- rownames(top_ABMR)
    tcmr_rn <- rownames(top_TCMR)
    
    Result$plota <- top_csv_a[abmr_rn, ]
    Result$plott <- top_csv_t[tcmr_rn, ]
    
    # Make sure the new data is in the same order as the variable names in the model
    ordered_abmr <- subset_abmr[var_names_abmr, ]
    ordered_tcmr <- subset_tcmr[var_names_tcmr, ]
    
    # Transpose the new data
    t_abmr <- t(ordered_abmr)
    t_tcmr <- t(ordered_tcmr)
    # Use the final_svm_model to predict class probabilities
    # Make sure to set the probability argument to TRUE
    Result$abmr <- predict(final_rf_model_3000, t_abmr[-1, ],type = "prob")[1]  # Predict ABMR risk
    Result$tcmr <- predict(final_rf_model_tcmr, t_tcmr[-1, ],type = "prob")[2]  # Predict TCMR risk
  })
  
  # Percentage output
  output$Result_A <- renderValueBox({
    (if (Result$abmr == -1) {
      valueBox(
        paste0("--%"), "Risk of ABMR", icon = icon("percent", lib = "font-awesome"),
        color = "green"
      )
    } else if (Result$abmr >= 0 & Result$abmr <= 0.3) {
      valueBox(
        paste0(round(Result$abmr * 100, 2), "%"), "Risk of ABMR", icon = icon("percent", lib = "font-awesome"),
        color = "green"
      )
    } else if (Result$abmr > 0.3 & Result$abmr <= 0.6) {
      valueBox(
        paste0(round(Result$abmr * 100, 2), "%"), "Risk of ABMR", icon = icon("percent", lib = "font-awesome"),
        color = "orange"
      )
    } else if (Result$abmr > 0.6 & Result$abmr <= 100) {
      valueBox(
        paste0(round(Result$abmr * 100, 2), "%"), "Risk of ABMR", icon = icon("percent", lib = "font-awesome"),
        color = "red"
      )
    })
  })
  
  output$Result_T <- renderValueBox({
    (if (Result$tcmr == -1) {
      valueBox(
        paste0("--%"), "Risk of TCMR", icon = icon("percent", lib = "font-awesome"),
        color = "green"
      )
    } else if (Result$tcmr >= 0 & Result$tcmr <= 0.3) {
      valueBox(
        paste0(round(Result$tcmr * 100, 2), "%"), "Risk of TCMR", icon = icon("percent", lib = "font-awesome"),
        color = "green"
      )
    } else if (Result$tcmr > 0.3 & Result$tcmr <= 0.6) {
      valueBox(
        paste0(round(Result$tcmr * 100, 2), "%"), "Risk of TCMR", icon = icon("percent", lib = "font-awesome"),
        color = "orange"
      )
    } else if (Result$tcmr > 0.6 & Result$tcmr <= 100) {
      valueBox(
        paste0(round(Result$tcmr * 100, 2), "%"), "Risk of TCMR", icon = icon("percent", lib = "font-awesome"),
        color = "red"
      )
    })
  })

  # Text output
  output$text <- renderText({
    Result$G <- paste("<p style='font-family:monaco;'> This tool allows for further confirmation of the diagnosis made by the histologist using a 
                      microscope. </p>", 
                      "<p style='font-family:monaco;'> Using genetic analysis means that the likelihood of another histologist coming to a
                      different conclusion about your biopsy is lessened, as this model studies the genetic presentation of ABMR and TCMR, rather
                      than only the physical presentations. </p>", 
                      sep = "\n" )
    (if (input$child == FALSE) {
      Result$M <- paste("<p style='font-family:monaco;'> Immunosuppressant medication reduces your immune system's response to your kidney being 
                          foreign object. This prevents your immune system from damaging your kidney. </p>", 
                        "<p style='font-family:monaco;'> No matter the outcome of your genetic analysis, continuing to take your immunosuppressant 
                          medication in the way that your medical team has advised is extremely important. </p>", 
                        "<p style='font-family:monaco;'> Even if no signs of rejection are perceived, this is BECAUSE of your adherence to
                        medication and not an indication to stop taking it. </p>", 
                        "<p style='font-family:monaco;'> If you are having uncomfortable side effects, please talk to your doctors and specialists, 
                          so that they can be managed. </p>", 
                        sep = "\n")
      if (Result$abmr == -1 && Result$tcmr == -1) {
        paste("<h3 style='font-family:optima;'> Upload a patient's data and click 'Result' button to predict. </h3>",
              "<h3 style='font-family:optima;'> Please wait until it says 'Upload complete' before clicking the 'Result' button. </h3>",
              "<h3 style='font-family:optima;'> Tick 'Child Patient' box to show result for child patient's parents. </h3>",
              "<h3 style='font-family:optima;'> Click 'Download' button to get a pdf report for patient. </h3>",
              sep = "\n")
      } else if (Result$abmr >= 0 && Result$abmr <= 0.3 && Result$tcmr >= 0 && Result$tcmr <= 0.3) {
        Result$R <- paste("<h1 style='font-family:optima;'> Non-Rejection </h1>", 
                          "<p style='font-family:monaco;'> Your immune system is functioning well with your new kidney! </p>", 
                           "<p style='font-family:monaco;'> Functioning is adequate and no significant signs of ABMR or TCMR were observed. </p>",
                           sep = "\n")
        paste(
          "<h1 style='font-family:optima;'> Non-Rejection </h1>", 
          "<p style='font-family:monaco;'> Your immune system is functioning well with your new kidney! </p>", 
          "<p style='font-family:monaco;'> Functioning is adequate and no significant signs of ABMR or TCMR were observed. </p>", 
          
          "<h1 style='font-family:optima;'> Medication Importance </h3>", 
          "<p style='font-family:monaco;'> Immunosuppressant medication reduces your immune system's response to your kidney being 
                          foreign object. This prevents your immune system from damaging your kidney. </p>", 
          "<p style='font-family:monaco;'> No matter the outcome of your genetic analysis, continuing to take your immunosuppressant 
          medication in the way that your medical team has advised is extremely important. </p>", 
          "<p style='font-family:monaco;'> Even if no signs of rejection are perceived, this is BECAUSE of your adherence to medication 
          and not an indication to stop taking it. </p>", 
          "<p style='font-family:monaco;'> If you are having uncomfortable side effects, please talk to your doctors and specialists, 
          so that they can be managed. </p>", 
          
          "<h1 style='font-family:optima;'> GenePath function and purpose </h3>", 
          "<p style='font-family:monaco;'> This tool allows for further confirmation of the diagnosis made by the histologist using a 
          microscope. </p>", 
          "<p style='font-family:monaco;'> Using genetic analysis means that the likelihood of another histologist coming to a different 
          conclusion about your biopsy is lessened, as this model studies the genetic presentation of ABMR and TCMR, rather than only the 
          physical presentations. </p>", 
          
          sep = "\n" 
        )
      } else if (Result$abmr > 0.3 && Result$tcmr > 0.3) {
        Result$R <- paste("<h1 style='font-family:optima;'> Mixed Rejection </h1>", 
                          "<p style='font-family:monaco;'> Your immune system has recognised the kidney as foreign and has activated the T-cells
                          and B-cells in your immune system to target your donor's genetic material. </p>", 
                          "<p style='font-family:monaco;'> This means that your immune system is directing T-cells to your kidney and producing
                          donor-specific antibodies, leading to inflammation in your kidney. </p>", 
                          "<p style='font-family:monaco;'> This results in damage on a cellular and structural level, which has serious
                          implications for graft functioning. </p>", sep = "\n")
        paste(
          "<h1 style='font-family:optima;'> Mixed Rejection </h1>", 
          "<p style='font-family:monaco;'> Your immune system has recognised the kidney as foreign and has activated the T-cells and B-cells 
          in your immune system to target your donor's genetic material. </p>", 
          "<p style='font-family:monaco;'> This means that your immune system is directing T-cells to your kidney and producing donor-specific 
          antibodies, leading to inflammation in your kidney. </p>", 
          "<p style='font-family:monaco;'> This results in damage on a cellular and structural level, which has serious implications for 
          graft functioning. </p>", 
          
          "<h1 style='font-family:optima;'> Medication Importance </h3>", 
          "<p style='font-family:monaco;'> Immunosuppressant medication reduces your immune system's response to your kidney being 
          foreign object. This prevents your immune system from damaging your kidney. </p>", 
          "<p style='font-family:monaco;'> No matter the outcome of your genetic analysis, continuing to take your immunosuppressant 
          medication in the way that your medical team has advised is extremely important. </p>", 
          "<p style='font-family:monaco;'> Even if no signs of rejection are perceived, this is BECAUSE of your adherence to medication 
          and not an indication to stop taking it. </p>", 
          "<p style='font-family:monaco;'> If you are having uncomfortable side effects, please talk to your doctors and specialists, 
          so that they can be managed. </p>", 
          
          "<h1 style='font-family:optima;'> GenePath function and purpose </h3>", 
          "<p style='font-family:monaco;'> This tool allows for further confirmation of the diagnosis made by the histologist using a 
          microscope. </p>", 
          "<p style='font-family:monaco;'> Using genetic analysis means that the likelihood of another histologist coming to a different 
          conclusion about your biopsy is lessened, as this model studies the genetic presentation of ABMR and TCMR, rather than only the 
          physical presentations. </p>", 
          
          sep = "\n" 
        )
      } else if (Result$abmr > 0.3 && Result$tcmr <= 0.3) {
        Result$R <- paste("<h1 style='font-family:optima;'> ABMR </h1>", 
                          "<p style='font-family:monaco;'> Your immune system has recognised the kidney as foreign and has produced antibodies which target 
                          your donor's genetic material. </p>", 
                          "<p style='font-family:monaco;'> This process is mediated by B-cells. This causes damage to the tissues and blood vessels of your 
                          kidney. </p>", 
                          "<p style='font-family:monaco;'> When the walls of the blood vessels are thickened, blood flow, and therefore kidney function, is 
                          impaired, meaning it has serious implications for graft loss. </p>", sep = "\n")
        paste(
          "<h1 style='font-family:optima;'> ABMR </h1>", 
          "<p style='font-family:monaco;'> Your immune system has recognised the kidney as foreign and has produced antibodies which target 
          your donor's genetic material. </p>", 
          "<p style='font-family:monaco;'> This process is mediated by B-cells. This causes damage to the tissues and blood vessels of your 
          kidney. </p>", 
          "<p style='font-family:monaco;'> When the walls of the blood vessels are thickened, blood flow, and therefore kidney function, is 
          impaired, meaning it has serious implications for graft loss. </p>", 
          
          "<h1 style='font-family:optima;'> Medication Importance </h3>", 
          "<p style='font-family:monaco;'> Immunosuppressant medication reduces your immune system's response to your kidney being 
                          foreign object. This prevents your immune system from damaging your kidney. </p>", 
          "<p style='font-family:monaco;'> No matter the outcome of your genetic analysis, continuing to take your immunosuppressant 
          medication in the way that your medical team has advised is extremely important. </p>", 
          "<p style='font-family:monaco;'> Even if no signs of rejection are perceived, this is BECAUSE of your adherence to medication 
          and not an indication to stop taking it. </p>", 
          "<p style='font-family:monaco;'> If you are having uncomfortable side effects, please talk to your doctors and specialists, 
          so that they can be managed. </p>", 
          
          "<h1 style='font-family:optima;'> GenePath function and purpose </h3>", 
          "<p style='font-family:monaco;'> This tool allows for further confirmation of the diagnosis made by the histologist using a 
          microscope. </p>", 
          "<p style='font-family:monaco;'> Using genetic analysis means that the likelihood of another histologist coming to a different 
          conclusion about your biopsy is lessened, as this model studies the genetic presentation of ABMR and TCMR, rather than only the 
          physical presentations. </p>", 
          
          sep = "\n" 
        )
      } else if (Result$abmr <= 0.3 && Result$tcmr > 0.3) {
        Result$R <- paste("<h1 style='font-family:optima;'> TCMR </h1>", 
                          "<p style='font-family:monaco;'> Your immune system has recognised the kidney as foreign and has activated the T-cells in your immune 
                          system to target your donor's genetic material. </p>", 
                          "<p style='font-family:monaco;'> This immune system response causes inflammation and cellular injury. This injury can extend to blood 
                          vessels and inhibit blood flow. </p>", 
                          "<p style='font-family:monaco;'> Prolonged TCMR can lead to tissue scarring, which is known as fibrosis, and disrupt the structure 
                          and function of the kidney. </p>", 
                          "<p style='font-family:monaco;'> The extent of the damage is dependent on the severity of the TCMR. </p>", sep = "\n")
        paste(
          "<h1 style='font-family:optima;'> TCMR </h1>", 
          "<p style='font-family:monaco;'> Your immune system has recognised the kidney as foreign and has activated the T-cells in your immune 
          system to target your donor's genetic material. </p>", 
          "<p style='font-family:monaco;'> This immune system response causes inflammation and cellular injury. This injury can extend to blood 
          vessels and inhibit blood flow. </p>", 
          "<p style='font-family:monaco;'> Prolonged TCMR can lead to tissue scarring, which is known as fibrosis, and disrupt the structure 
          and function of the kidney. </p>", 
          "<p style='font-family:monaco;'> The extent of the damage is dependent on the severity of the TCMR. </p>", 
          
          "<h1 style='font-family:optima;'> Medication Importance </h3>", 
          "<p style='font-family:monaco;'> Immunosuppressant medication reduces your immune system's response to your kidney being 
                          foreign object. This prevents your immune system from damaging your kidney. </p>", 
          "<p style='font-family:monaco;'> No matter the outcome of your genetic analysis, continuing to take your immunosuppressant 
          medication in the way that your medical team has advised is extremely important. </p>", 
          "<p style='font-family:monaco;'> Even if no signs of rejection are perceived, this is BECAUSE of your adherence to medication 
          and not an indication to stop taking it. </p>", 
          "<p style='font-family:monaco;'> If you are having uncomfortable side effects, please talk to your doctors and specialists, 
          so that they can be managed. </p>", 
          
          "<h1 style='font-family:optima;'> GenePath function and purpose </h3>", 
          "<p style='font-family:monaco;'> This tool allows for further confirmation of the diagnosis made by the histologist using a 
          microscope. </p>", 
          "<p style='font-family:monaco;'> Using genetic analysis means that the likelihood of another histologist coming to a different 
          conclusion about your biopsy is lessened, as this model studies the genetic presentation of ABMR and TCMR, rather than only the 
          physical presentations. </p>", 
          
          sep = "\n" 
        )
      }
    } else {
      Result$M <- paste("<p style='font-family:monaco;'> Immunosuppressant medication reduces your child's immune system's response to their
                          kidney being foreign object. This prevents their immune system from damaging their kidney. </p>", 
                        "<p style='font-family:monaco;'> No matter the outcome of your child's genetic analysis, continuing their immunosuppressant medication 
                          regime in the way that the medical team has advised is extremely important. </p>", 
                        "<p style='font-family:monaco;'> Even if no signs of rejection are perceived, this is BECAUSE of your child's adherence to medication 
                          and not an indication to stop taking it. </p>", 
                        "<p style='font-family:monaco;'> If your child is having uncomfortable side effects, please talk to their doctors and specialists, so 
                          that they can be managed. </p>", 
                        "<p style='font-family:monaco;'> Keep an open conversation with them so that you can facilitate this for them. </p>", 
                        sep = "\n")
      if (Result$abmr == -1 && Result$tcmr == -1) {
        paste("<h3 style='font-family:optima;'> Upload a patient's data and click 'Result' button to predict. </h3>",
              "<h3 style='font-family:optima;'> Please wait until it says 'Upload complete' before clicking the 'Result' button. </h3>",
              "<h3 style='font-family:optima;'> Tick 'Child Patient' box to show result for child patient's parents. </h3>",
              "<h3 style='font-family:optima;'> Click 'Download' button to get a pdf report for patient. </h3>",
              sep = "\n")
      } else if (Result$abmr >= 0 && Result$abmr <= 0.3 && Result$tcmr >= 0 && Result$tcmr <= 0.3) {
        Result$R <- paste("<h1 style='font-family:optima;'> Non-Rejection </h1>", 
                          "<p style='font-family:monaco;'> Your child's immune system is functioning well with their new kidney! </p>", 
                          "<p style='font-family:monaco;'> Functioning is adequate and no significant signs of ABMR or TCMR were observed. </p>", sep = "\n")
        paste(
          "<h1 style='font-family:optima;'> Non-Rejection </h1>", 
          "<p style='font-family:monaco;'> Your child's immune system is functioning well with their new kidney! </p>", 
          "<p style='font-family:monaco;'> Functioning is adequate and no significant signs of ABMR or TCMR were observed. </p>", 
          
          "<h1 style='font-family:optima;'> Medication Importance </h3>", 
          "<p style='font-family:monaco;'> Immunosuppressant medication reduces your child's immune system's response to their kidney being foreign 
          object. This prevents their immune system from damaging their kidney. </p>", 
          "<p style='font-family:monaco;'> No matter the outcome of your child's genetic analysis, continuing their immunosuppressant medication 
          regime in the way that the medical team has advised is extremely important. </p>", 
          "<p style='font-family:monaco;'> Even if no signs of rejection are perceived, this is BECAUSE of your child's adherence to medication 
          and not an indication to stop taking it. </p>", 
          "<p style='font-family:monaco;'> If your child is having uncomfortable side effects, please talk to their doctors and specialists, so 
          that they can be managed. </p>", 
          "<p style='font-family:monaco;'> Keep an open conversation with them so that you can facilitate this for them. </p>", 
          
          "<h1 style='font-family:optima;'> GenePath function and purpose </h3>", 
          "<p style='font-family:monaco;'> This tool allows for further confirmation of the diagnosis made by the histologist using a 
          microscope. </p>", 
          "<p style='font-family:monaco;'> Using genetic analysis means that the likelihood of another histologist coming to a different 
          conclusion about your biopsy is lessened, as this model studies the genetic presentation of ABMR and TCMR, rather than only the 
          physical presentations. </p>", 
          
          sep = "\n" 
        )
      } else if (Result$abmr > 0.3 && Result$tcmr > 0.3) {
        Result$R <- paste("<h1 style='font-family:optima;'> Mixed Rejection </h1>", 
                          "<p style='font-family:monaco;'> Your child's immune system has recognised the kidney as foreign and has activated the T-cells and 
                          B-cells in their immune system to target your child's donor's genetic material. </p>", 
                          "<p style='font-family:monaco;'> This means that their immune system is directing T-cells to your kidney and producing donor-specific 
                          antibodies, leading to inflammation in their kidney. </p>", 
                          "<p style='font-family:monaco;'> This results in damage on a cellular and structural level, which has serious implications for graft 
                          functioning. </p>", sep = "\n")
        paste(
          "<h1 style='font-family:optima;'> Mixed Rejection </h1>", 
          "<p style='font-family:monaco;'> Your child's immune system has recognised the kidney as foreign and has activated the T-cells and 
          B-cells in their immune system to target your child's donor's genetic material. </p>", 
          "<p style='font-family:monaco;'> This means that their immune system is directing T-cells to your kidney and producing donor-specific 
          antibodies, leading to inflammation in their kidney. </p>", 
          "<p style='font-family:monaco;'> This results in damage on a cellular and structural level, which has serious implications for graft 
          functioning. </p>", 
          
          "<h1 style='font-family:optima;'> Medication Importance </h3>", 
          "<p style='font-family:monaco;'> Immunosuppressant medication reduces your child's immune system's response to their kidney being foreign 
          object. This prevents their immune system from damaging their kidney. </p>", 
          "<p style='font-family:monaco;'> No matter the outcome of your child's genetic analysis, continuing their immunosuppressant medication 
          regime in the way that the medical team has advised is extremely important. </p>", 
          "<p style='font-family:monaco;'> Even if no signs of rejection are perceived, this is BECAUSE of your child's adherence to medication 
          and not an indication to stop taking it. </p>", 
          "<p style='font-family:monaco;'> If your child is having uncomfortable side effects, please talk to their doctors and specialists, so 
          that they can be managed. </p>", 
          "<p style='font-family:monaco;'> Keep an open conversation with them so that you can facilitate this for them. </p>", 
          
          "<h1 style='font-family:optima;'> GenePath function and purpose </h3>", 
          "<p style='font-family:monaco;'> This tool allows for further confirmation of the diagnosis made by the histologist using a 
          microscope. </p>", 
          "<p style='font-family:monaco;'> Using genetic analysis means that the likelihood of another histologist coming to a different 
          conclusion about your biopsy is lessened, as this model studies the genetic presentation of ABMR and TCMR, rather than only the 
          physical presentations. </p>", 
          
          sep = "\n" 
        )
      } else if (Result$abmr > 0.3 && Result$tcmr <= 0.3) {
        Result$R <- paste("<h1 style='font-family:optima;'> ABMR </h1>", 
                          "<p style='font-family:monaco;'> Your child's immune system has recognised the kidney as foreign and has produced antibodies which 
                          target your child's donor's genetic material. </p>", 
                          "<p style='font-family:monaco;'> This process is mediated by B-cells. This causes damage to the tissues and blood vessels of their 
                          kidney. </p>", 
                          "<p style='font-family:monaco;'> When the walls of the blood vessels are thickened, blood flow, and therefore kidney function, is 
                          impaired, meaning it has serious implications for graft loss. </p>", sep = "\n")
        paste(
          "<h1 style='font-family:optima;'> ABMR </h1>", 
          "<p style='font-family:monaco;'> Your child's immune system has recognised the kidney as foreign and has produced antibodies which 
          target your child's donor's genetic material. </p>", 
          "<p style='font-family:monaco;'> This process is mediated by B-cells. This causes damage to the tissues and blood vessels of their 
          kidney. </p>", 
          "<p style='font-family:monaco;'> When the walls of the blood vessels are thickened, blood flow, and therefore kidney function, is 
          impaired, meaning it has serious implications for graft loss. </p>", 
          
          "<h1 style='font-family:optima;'> Medication Importance </h3>", 
          "<p style='font-family:monaco;'> Immunosuppressant medication reduces your child's immune system's response to their kidney being foreign 
          object. This prevents their immune system from damaging their kidney. </p>", 
          "<p style='font-family:monaco;'> No matter the outcome of your child's genetic analysis, continuing their immunosuppressant medication 
          regime in the way that the medical team has advised is extremely important. </p>", 
          "<p style='font-family:monaco;'> Even if no signs of rejection are perceived, this is BECAUSE of your child's adherence to medication 
          and not an indication to stop taking it. </p>", 
          "<p style='font-family:monaco;'> If your child is having uncomfortable side effects, please talk to their doctors and specialists, so 
          that they can be managed. </p>", 
          "<p style='font-family:monaco;'> Keep an open conversation with them so that you can facilitate this for them. </p>", 
          
          "<h1 style='font-family:optima;'> GenePath function and purpose </h3>", 
          "<p style='font-family:monaco;'> This tool allows for further confirmation of the diagnosis made by the histologist using a 
          microscope. </p>", 
          "<p style='font-family:monaco;'> Using genetic analysis means that the likelihood of another histologist coming to a different 
          conclusion about your biopsy is lessened, as this model studies the genetic presentation of ABMR and TCMR, rather than only the 
          physical presentations. </p>", 
          
          sep = "\n" 
        )
      } else if (Result$abmr <= 0.3 && Result$tcmr > 0.3) {
        Result$R <- paste("<h1 style='font-family:optima;'> TCMR </h1>", 
                          "<p style='font-family:monaco;'> Your child's immune system has recognised the kidney as foreign and has activated the T-cells in your 
                          immune system to target your child's donor's genetic material. </p>", 
                          "<p style='font-family:monaco;'> This immune system response causes inflammation and cellular injury. This injury can extend to blood 
                          vessels and inhibit blood flow. </p>", 
                          "<p style='font-family:monaco;'> Prolonged TCMR can lead to tissue scarring, which is known as fibrosis, and disrupt the structure and 
                          function of the kidney. </p>", 
                          "<p style='font-family:monaco;'> The extent of the damage is dependent on the severity of the TCMR. </p>", sep = "\n")
        paste(
          "<h1 style='font-family:optima;'> TCMR </h1>", 
          "<p style='font-family:monaco;'> Your child's immune system has recognised the kidney as foreign and has activated the T-cells in your 
          immune system to target your child's donor's genetic material. </p>", 
          "<p style='font-family:monaco;'> This immune system response causes inflammation and cellular injury. This injury can extend to blood 
          vessels and inhibit blood flow. </p>", 
          "<p style='font-family:monaco;'> Prolonged TCMR can lead to tissue scarring, which is known as fibrosis, and disrupt the structure and 
          function of the kidney. </p>", 
          "<p style='font-family:monaco;'> The extent of the damage is dependent on the severity of the TCMR. </p>", 
          
          "<h1 style='font-family:optima;'> Medication Importance </h3>", 
          "<p style='font-family:monaco;'> Immunosuppressant medication reduces your child's immune system's response to their kidney being foreign 
          object. This prevents their immune system from damaging their kidney. </p>", 
          "<p style='font-family:monaco;'> No matter the outcome of your child's genetic analysis, continuing their immunosuppressant medication 
          regime in the way that the medical team has advised is extremely important. </p>", 
          "<p style='font-family:monaco;'> Even if no signs of rejection are perceived, this is BECAUSE of your child's adherence to medication 
          and not an indication to stop taking it. </p>", 
          "<p style='font-family:monaco;'> If your child is having uncomfortable side effects, please talk to their doctors and specialists, so 
          that they can be managed. </p>", 
          "<p style='font-family:monaco;'> Keep an open conversation with them so that you can facilitate this for them. </p>", 
          
          "<h1 style='font-family:optima;'> GenePath function and purpose </h3>", 
          "<p style='font-family:monaco;'> This tool allows for further confirmation of the diagnosis made by the histologist using a 
          microscope. </p>", 
          "<p style='font-family:monaco;'> Using genetic analysis means that the likelihood of another histologist coming to a different 
          conclusion about your biopsy is lessened, as this model studies the genetic presentation of ABMR and TCMR, rather than only the 
          physical presentations. </p>", 
          
          sep = "\n" 
        )
      }
    })
  })
  
  # Boxplot for Top5 Genes
  top5_ABMR <- ordered_gene_data[1:5,]
  top5_TCMR <- ordered_gene_data_tcmr[1:5,]
  
  colors = c("tomato2", "limegreen")
  
  abmr1 <- data.frame(gene_expression = top5_ABMR[1,],
                      ABMR_Outcome = factor(ABMR_Outcome))
  abmr2 <- data.frame(gene_expression = top5_ABMR[2,],
                      ABMR_Outcome = factor(ABMR_Outcome))
  abmr3 <- data.frame(gene_expression = top5_ABMR[3,],
                      ABMR_Outcome = factor(ABMR_Outcome))
  abmr4 <- data.frame(gene_expression = top5_ABMR[4,],
                      ABMR_Outcome = factor(ABMR_Outcome))
  abmr5 <- data.frame(gene_expression = top5_ABMR[5,],
                      ABMR_Outcome = factor(ABMR_Outcome))
  
  tcmr1 <- data.frame(gene_expression = top5_TCMR[1,],
                      TCMR_Outcome = factor(TCMR_Outcome, levels = c("TCMR", "Non_rejection")))
  tcmr2 <- data.frame(gene_expression = top5_TCMR[2,],
                      TCMR_Outcome = factor(TCMR_Outcome, levels = c("TCMR", "Non_rejection")))
  tcmr3 <- data.frame(gene_expression = top5_TCMR[3,],
                      TCMR_Outcome = factor(TCMR_Outcome, levels = c("TCMR", "Non_rejection")))
  tcmr4 <- data.frame(gene_expression = top5_TCMR[4,],
                      TCMR_Outcome = factor(TCMR_Outcome, levels = c("TCMR", "Non_rejection")))
  tcmr5 <- data.frame(gene_expression = top5_TCMR[5,],
                      TCMR_Outcome = factor(TCMR_Outcome, levels = c("TCMR", "Non_rejection")))
  
  output$boxa1 <- renderPlot({
    if (is.null(Result$plota[1, -1])) {
      ggplot(abmr1, aes(x = ABMR_Outcome, y = gene_expression, fill = ABMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_ABMR)[1]) + 
        ylab("Gene Expressions") 
    } else {
      ggplot(abmr1, aes(x = ABMR_Outcome, y = gene_expression, fill = ABMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        geom_point(aes(y = Result$plota[1, -1]), color = "blue", size = 3) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_ABMR)[1]) + 
        ylab("Gene Expressions") 
    }
  })
  
  output$boxa2 <- renderPlot({
    if (is.null(Result$plota[2, -1])) {
      ggplot(abmr2, aes(x = ABMR_Outcome, y = gene_expression, fill = ABMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_ABMR)[2]) + 
        ylab("Gene Expressions") 
    } else {
      ggplot(abmr2, aes(x = ABMR_Outcome, y = gene_expression, fill = ABMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        geom_point(aes(y = Result$plota[2, -1]), color = "blue", size = 3) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_ABMR)[2]) + 
        ylab("Gene Expressions") 
    }
  })
  
  output$boxa3 <- renderPlot({
    if (is.null(Result$plota[3, -1])) {
      ggplot(abmr3, aes(x = ABMR_Outcome, y = gene_expression, fill = ABMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_ABMR)[3]) + 
        ylab("Gene Expressions") 
    } else {
      ggplot(abmr3, aes(x = ABMR_Outcome, y = gene_expression, fill = ABMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        geom_point(aes(y = Result$plota[3, -1]), color = "blue", size = 3) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_ABMR)[3]) + 
        ylab("Gene Expressions") 
    }
  })
  
  output$boxa4 <- renderPlot({
    if (is.null(Result$plota[4, -1])) {
      ggplot(abmr4, aes(x = ABMR_Outcome, y = gene_expression, fill = ABMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_ABMR)[4]) + 
        ylab("Gene Expressions") 
    } else {
      ggplot(abmr4, aes(x = ABMR_Outcome, y = gene_expression, fill = ABMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        geom_point(aes(y = Result$plota[4, -1]), color = "blue", size = 3) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_ABMR)[4]) + 
        ylab("Gene Expressions") 
    }
  })
  
  output$boxa5 <- renderPlot({
    if (is.null(Result$plota[5, -1])) {
      ggplot(abmr5, aes(x = ABMR_Outcome, y = gene_expression, fill = ABMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_ABMR)[5]) + 
        ylab("Gene Expressions") 
    } else {
      ggplot(abmr5, aes(x = ABMR_Outcome, y = gene_expression, fill = ABMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        geom_point(aes(y = Result$plota[5, -1]), color = "blue", size = 3) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_ABMR)[5]) + 
        ylab("Gene Expressions") 
    }
  })
  
  output$boxt1 <- renderPlot({
    if (is.null(Result$plott[1, -1])) {
      ggplot(tcmr1, aes(x = TCMR_Outcome, y = gene_expression, fill = TCMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_TCMR)[1]) + 
        ylab("Gene Expressions") 
    } else {
      ggplot(tcmr1, aes(x = TCMR_Outcome, y = gene_expression, fill = TCMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        geom_point(aes(y = Result$plott[1, -1]), color = "blue", size = 3) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_TCMR)[1]) + 
        ylab("Gene Expressions") 
    }
  })
  
  output$boxt2 <- renderPlot({
    if (is.null(Result$plott[2, -1])) {
      ggplot(tcmr2, aes(x = TCMR_Outcome, y = gene_expression, fill = TCMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_TCMR)[2]) + 
        ylab("Gene Expressions") 
    } else {
      ggplot(tcmr2, aes(x = TCMR_Outcome, y = gene_expression, fill = TCMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        geom_point(aes(y = Result$plott[2, -1]), color = "blue", size = 3) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_TCMR)[2]) + 
        ylab("Gene Expressions") 
    }
  })
  
  output$boxt3 <- renderPlot({
    if (is.null(Result$plott[3, -1])) {
      ggplot(tcmr3, aes(x = TCMR_Outcome, y = gene_expression, fill = TCMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_TCMR)[3]) + 
        ylab("Gene Expressions") 
    } else {
      ggplot(tcmr3, aes(x = TCMR_Outcome, y = gene_expression, fill = TCMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        geom_point(aes(y = Result$plott[3, -1]), color = "blue", size = 3) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_TCMR)[3]) + 
        ylab("Gene Expressions") 
    }
  })
  
  output$boxt4 <- renderPlot({
    if (is.null(Result$plott[4, -1])) {
      ggplot(tcmr4, aes(x = TCMR_Outcome, y = gene_expression, fill = TCMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_TCMR)[4]) + 
        ylab("Gene Expressions") 
    } else {
      ggplot(tcmr4, aes(x = TCMR_Outcome, y = gene_expression, fill = TCMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        geom_point(aes(y = Result$plott[4, -1]), color = "blue", size = 3) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_TCMR)[4]) + 
        ylab("Gene Expressions") 
    }
  })
  
  output$boxt5 <- renderPlot({
    if (is.null(Result$plott[5, -1])) {
      ggplot(tcmr5, aes(x = TCMR_Outcome, y = gene_expression, fill = TCMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_TCMR)[5]) + 
        ylab("Gene Expressions") 
    } else {
      ggplot(tcmr5, aes(x = TCMR_Outcome, y = gene_expression, fill = TCMR_Outcome)) + 
        geom_boxplot(notch = TRUE) + 
        geom_point(aes(y = Result$plott[5, -1]), color = "blue", size = 3) + 
        scale_fill_manual(values = colors, guide = "none") + 
        xlab(rownames(top5_TCMR)[5]) + 
        ylab("Gene Expressions") 
    }
  })
  
  # download a pdf report
  output$download <- downloadHandler(
    filename = "GenePath_Report.pdf",
    content = function(file) {
      rmarkdown::render("report.Rmd", output_file = file, params =  list(abmr = Result$abmr, tcmr = Result$tcmr, r = Result$R, m = Result$M, g = Result$G))
    }
  )
}

# Run the application 
shinyApp(ui, server) 