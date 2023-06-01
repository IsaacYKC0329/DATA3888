# Gene Path

## Starting the app 

1. Our app has been published through shiny.io. You can access our webpage using this link: https://kidneya12.shinyapps.io/DATA3888_Kidney_A12/
2. If the link is not working, you can open the app.R file in this repository. Click 'Run App' button on the top right side. 
   - Please make sure you have **installed all required packages**! 
   - If the model is not working, please **run these in your R console** before run the app: 
   
   `load("ABMR_RF_model_3000.RData")`
   
   `load("ABMR_Outcome.RData")`
   
   `load("ordered_gene_data.RData")`
   
   `load("TCMR_RF_model_1000.RData")`
   
   `load("TCMR_Outcome.RData")`
   
   `load("ordered_gene_data_tcmr.RData")`

## Functions

### Predictor Page

The 'Predictor' Page can determine the probability of patients experiencing rejection.

##### Predictor

Upload a patient's data and click 'Result' button to predict. Please **wait until it says 'Upload complete'** before clicking the 'Result' button. 

Sample CSV files are in the 'Test csv files' repository. 

##### Child Patient

Select 'Child Patient' check box to transform the resulting output into suggestions for the patient’s patients (or guardians).

##### Download

After the prediction has been made, click 'Download' to generate a PDF report for the patients, including the displayed result. 

### Top Genes Page

The 'Top Genes' page default shows the top 5 gene expression of ABMR and TCMR against Non-rejection in boxplots. 

After the prediction has been made, The blue points represent the genetic signals of the patient. 
