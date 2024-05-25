#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

######################
#       IMPORTS      #
######################

library(shiny)
library(DT)
library(matrixStats)
library(tidyverse)
library(pheatmap)
library(colourpicker)
library(fgsea)
library(beeswarm)
library(DESeq2)
library(data.table)

# Change this number in order to change the max upload size for files
options(shiny.maxRequestSize = 200 * 1024^2)  # Set max upload size to 100 MB

###################################################
#                 FRONT END                       #
###################################################

# This section of the code contains all of the UI for the code
# Rshiny code is written in a similar format to HTML, but it's all stored in the UI variable.

####################
# Frontend Variables
####################
# Variables that are used in the dropdowns in the UI and need to be defined beforehand

deseq_choices <-
  c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

####################
# UI
####################

ui <- fluidPage(
  # Title 
  titlePanel("Plotting of Normalized Counts Data"),
  
  # Beginning of all tabs
  tabsetPanel(
    ####################################################
    # SAMPLES TAB: where you look at metadata
    ####################################################
    tabPanel("Samples",
      ####################################################
      # Samples Sidebar
      ####################################################
       sidebarPanel(
         fileInput("samples_file", "Choose a CSV file for sample information exploration:"),
         
         selectInput("numeric_columns", "Select Numeric Column", choices =c(
           "age_of_death","AvgSpotLen","mrna.seq_reads","pmi","RIN", "age_of_onset", "cag", "duration", "h.v_cortical_score", "h.v_striatal_score", "vonsattel_grade")
         ), # TODO: Make these choices not hardcoded, have the same UI_output as below 
         
         uiOutput("samples_columns"), # dropdown menu that shows all of the columns in the csv
         
         submitButton(
           text = "Submit",
           icon = icon("car-crash"),
           width = "100%"
         ) # Submit button
         
        ),
      ####################################################
      # Samples Main Panel
      ####################################################
       mainPanel(
         tabsetPanel(
           
           # Table tab has a Dataframe output of the samples table given
           tabPanel("Table", 
                    DTOutput("samples_table")),
           
           # Select columns gives the table with the columns that you selected 
           tabPanel("Select Columns", 
                    DTOutput("filtered_samples_table")),
           
           # Summary tab has a summary of all the selected columns, their type (int, string, etc.) and most common entry
           tabPanel("Summary", 
                    tableOutput("samples_summary")),
           
           # Plots numeric column vs selected column
           # TODO: Currently the numeric columns are hardcoded so you can't select anything else
           tabPanel("Plots", 
                    plotOutput("selected_plot"))
         )
       )
    ),
    ####################################################
    # End of Samples Tab
    ####################################################
    
    ####################################################
    # Start of Counts Tab: Tab for filtering Counts Matrix
    ####################################################
    tabPanel("Counts", 
       ####################################################
       # Counts Sidebar
       ####################################################
       sidebarPanel(
         
         # CSV Input
         fileInput("counts_file", "Choose a CSV file for counts matrix exploration:"),
         
         # Choose how much variance you want the row to have to be selected
         sliderInput("counts_variance", "Select the percentile of variance:", min = 0, max = 100, value = 0),
         
         # Choose how many non zero samples you want the row to have to be selected
         sliderInput("counts_non_zero", "Select the number of samples that are non-zero to be included:", min = 0, max = 100, value = 0),
         
         # Submit Button
         submitButton(
           text = "Submit",
           icon = icon("car-crash"),
           width = "100%"
         )
       ),
       ####################################################
       # Counts Main Panel
       ####################################################
       mainPanel(
         tabsetPanel(
           
           # Shows the results of filtering from the sidebar
           # TODO: Show DTOutput of the filtered table underneath the counts filtering table
           tabPanel("Filtering Effect", tableOutput("counts_filtering_table")),
           
           # Shows how many samples actually pass the filters and where they are in terms of variance and # of zeros
           tabPanel("Diagnostic Scatter Plots", plotOutput("plot_variance"), plotOutput("plot_zeros")),
           
           # Shows Samples vs Genes heatmap
           tabPanel("Clustered Heatmap", plotOutput("counts_heatmap")),
           
           # Beeswarm of which components are most important from a PCA
           tabPanel("PCA", sliderInput("num_components", "Choose number of components for PCA",min = 0, max = 69, value = 2 ),
                    plotOutput("counts_pca"))
         )
       )
    ),
    
    ####################################################
    # End of Counts Tab
    ####################################################
    
    ####################################################
    # Start of DE Tab
    ####################################################
    
    # This tab will either let you insert differential expression results from a csv, or take the sample and counts
    # matrices from before and do differential expression.
    
    # TODO: Change the tab structure so Samples and Counts are under "DE Preprocessing" tab, 
    # and DE tab becomes "DE Analysis".
    
    tabPanel("DE", 
       ####################################################
       # DE Analysis Sidebar
       ####################################################
       sidebarPanel(
         # Insert DE CSV file
         # TODO: Make it so that there isn't an error when this is empty
         fileInput("DE_file", "Choose a CSV file for differential expression:"),
         
         # Radio buttons where the user can choose whether they want to use the new data or the old data
         # TODO: Automatically use old data unless they want to use new data
         # Possibly hide file input unless no is checked on this radio button?
         radioButtons(
           inputId = "deseq",
           label = "Run Deseq2 on the filtered counts data from previous data?",
           choices = c("yes", "no"),
           selected = "no"
         ),
         
         # The following are choices for the volcano plots
         # deseq_choices are at the top 
         # TODO: Make hovering the points give you the name of the gene
         # TODO: Make clicking the points take you to the uniprot page for the gene
         radioButtons(
           inputId = "x_axis",
           label = "Choose the column for the x-axis",
           choices = deseq_choices,
           selected = "log2FoldChange"
         ),
         radioButtons(
           inputId = "y_axis",
           label = "Choose the column for the y-axis",
           choices = deseq_choices,
           selected = "padj"
         ),
         colourInput(
           inputId = "base",
           label = "Base point color",
           value = "#22577A",
           closeOnClick = T
         ),
         colourInput(
           inputId = "highlight",
           label = "Highlight point color",
           value = "#FFCF56",
           closeOnClick = T
         ),
         sliderInput(
           "slider",
           "Select the magnitude of the p adjusted coloring:",
           min = -50,
           max = 0,
           value = -5,
           step = 1
         ),
         
         # Submit Button
         submitButton(
           text = "Plot",
           icon = icon("car-crash"),
           width = "100%"
         )
       ),
       
       ####################################################
       # DE Analysis Main Panel
       ####################################################
       mainPanel(
         tabsetPanel(
           # Show the DE Table
           tabPanel("Sortable Table", DTOutput("DE_table")),
           
           # Show the volcano plot and the corresponding table``
           tabPanel("DE Analysis",
                    tabsetPanel(tabPanel("Volcano Plots", plotOutput("volcano")), tabPanel("Table", DTOutput("table"))))
         )
       )
    ),
    ####################################################
    # End of DE Tab
    ####################################################
    
    ####################################################
    # Start of EA Tab
    ####################################################
    tabPanel("Enrichment Analysis", 
       tabsetPanel(
         ####################################################
         # EA Table
         ####################################################
         # This tab creates a table of the enrichment analysis that you can sort through
         # You can also select which NES pathways to select.
         
         tabPanel("Sortable Data Table", 
            ####################################################
            # EA Table Sidebar
            ####################################################
            sidebarPanel(
              # Insert Differential Expression CSV
              # TODO: Have tab automatically detect previous differential expression tab
              fileInput("EA_file", "Choose a CSV file for differential expression:"),
              
              # Slider to select max p adjusted for NES pathways
              sliderInput("EA_p_adj_table", "p adjusted to select pathways by:", min =-30, max = 0, value = 0),
              
              # Select whether you want negative, positive or all pathways
              radioButtons("EA_NES_select", "Choose which NES pathways to select:",
                           c("positive", "negative", "all")),
              
            
              submitButton(
                text = "Update Table",
                icon = icon("car-crash"),
                width = "100%"
              ),
              
              #Download button
              # TODO: Add download button to all other tabs with tables, see if possible with graphs
              p("Download the current filtered table:"),
              downloadButton("EA_download", "Download Table")
              
            ),
            ####################################################
            # EA Table Main Panel
            ####################################################
            mainPanel(
              DTOutput("EA_Table")
            )
         ),
         
         ####################################################
         # EA Top Pathways
         ####################################################
         tabPanel("Top Pathways", 
            ####################################################
            # EA Top Pathways Sidebar
            ####################################################
            sidebarPanel(
              fileInput("EA_file", "Choose a CSV file for differential expression:"),
              sliderInput("EA_p_adj", "p adjusted to select pathways by:", min = -30, max = 0, value = -10),
              submitButton(
                text = "Update Table",
                icon = icon("car-crash"),
                width = "100%"
              )
            ),
            ####################################################
            # EA Top Pathways Main Panel
            ####################################################
            mainPanel(
              plotOutput("EA_barplot")
            ),
                  
         ),
         ####################################################
         # EA Scatter
         ####################################################
         tabPanel("Scatter Plot of NES", 
            ####################################################
            # EA Scatter Sidebar
            ####################################################
            sidebarPanel(
              fileInput("EA_file", "Choose a CSV file for differential expression:"),
              sliderInput("EA_p_adj_plot", "p adjusted to select pathways by:", min = -30, max = 0, value = -10),
              submitButton(
                text = "Update Table",
                icon = icon("car-crash"),
                width = "100%"
              )
            ),
            ####################################################
            # EA Scatter Main Panel
            ####################################################
            mainPanel(
              plotOutput("EA_scatter_plot")
            )
         )
       )
    )
  )
)


###################################################
#                 BACK END                        #
###################################################

# All of the functions that are used in the frontend 
server <- function(input, output) {
  
  ####################################################
  # SAMPLES FUNCTIONS
  ####################################################
  
  # Functions used to create the outputs for the Samples tab
  
  # PURPOSE: Function to read the file input into "samples_file" variable of input.
  # INPUTS: samples_file is the variable in the frontend that accepts a csv
  # REACTIVE: As soon as a csv file is inserted into the variable, it is automatically read
  
  load_samples_data <- reactive({
    new_file <- fread(input$samples_file$datapath)
    return(as.tibble(new_file))
  })
  
  
  # PURPOSE: This function is used for the "select columns" tab of the samples tab. 
  # Selects columns from the full 
  # REACTIVE:
  # TODO: Not currently reactive, only changes when the submit button is clicked
  # Possibly because it references the load_samples_data()? Try to put load_samples_data() into 
  
  column_filter <- reactive({
    selected_data <- load_samples_data() %>%
      select(input$samples_columns)
    return(selected_data)
  })
  
  # Creates the summary for the samples tab
  samples_summary <-
    function(dataf) {
      newdf <- tibble(
        Columns = colnames(dataf),
        type = sapply(dataf, typeof),
        most_common = sapply(dataf, function(col) {
          if (typeof(col) %in% c("character", "factor")) {
            names(table(col))[which.max(table(col))]
          } else {
            mean(na.omit(col))
          }
        })
      )
      return(newdf)
    }
  
  
  ####################################################
  # SAMPLES OUTPUTS
  ####################################################
  
  # Outputs for the 
  
  output$samples_columns <- renderUI({
    # Get the column names of the data frame
    choices <- colnames(load_samples_data())
    # Create a selectInput with dynamically set choices
    selectInput("samples_columns", "Select Columns", choices, multiple = TRUE)
  })
  
  
  
  output$samples_summary <- renderTable(
    samples_summary(load_samples_data())
  ) 
  
  # Load the samples summary into an output object for the UI
  output$samples_table <- renderDT(
    datatable(load_samples_data(), options = list(ordering = TRUE))
  )
  
  # For the radio buttons to choose which things to plot
  output$filtered_samples_table <- renderDT({
    # Subset the data based on the selected columns
    # Render the table
    datatable(column_filter(), options = list(dom = 't',pageLength = nrow(column_filter())))
  })
  
  
  # Render the plot or display for the selected numeric column
  output$selected_plot <- renderPlot({
    # Example: Create a histogram for the selected numeric column
    data <- load_samples_data()
    vector <- input$numeric_columns
    # print(vector)
    hist(data[[vector]], main = paste("Histogram of", vector))
  })
  
  
  ####################################################
  # COUNTS FUNCTIONS
  ####################################################
  
  output$counts <- renderDT(
    datatable(load_counts_data(), options = list(ordering = TRUE))
  )
  
  # Make a table that summarizes the effect of filtering
  # Return a table with the number of samples, number of genes, number and % of genes passing filter and not passing filter
  counts_filtering_table <- function(dataf, variance, zero) {
    filtered_counts <- counts_filtering(variance, zero)
    
    # Create the summary table
    table <- tibble(
      `Number of samples` = ncol(dataf),
      `Total Number of Genes` = nrow(dataf),
      `Number of Genes passing the filter` = nrow(filtered_counts),
      `Percent of Genes passing the filter` = nrow(filtered_counts)/nrow(dataf) * 100,
      `Number of Genes not passing filter` = nrow(dataf) - nrow(filtered_counts),
      `Percent of Genes not passing the filter` =  ((nrow(dataf) - nrow(filtered_counts))/nrow(dataf))*100
    )
    
    # Return the summary table
    return(table)
  }
  
  output$counts_filtering_table <- renderTable({
    counts_filtering_table(load_counts_data(), input$counts_variance, input$counts_non_zero)
  })
  
  #Make diagnostic scatter plots
  # Genes passing in a darker color, genes not passing in a lighter color
  # Do median count vs variance and median count vs number of zeros
  counts_scatter_plot <- function(variance, non_zero) {
    dataf <- load_counts_data()
    filtered_counts <- counts_filtering(variance, non_zero)
    df_var <- rowVars(as.matrix(dataf[, -1]), na.rm = TRUE)
    df_zeros <- rowSums(dataf == 0)
    
    diagnostic_data <- cbind(dataf, df_var, df_zeros)
    diagnostic_data$filtered <- ifelse(df_var >= variance & df_zeros>=non_zero, "Pass", "Fail")
    
    # Create scatter plot for median count vs variance
    plot_variance <- ggplot(diagnostic_data, aes(x = rank(rowMeans(dataf[, -1], na.rm = TRUE)), y = log2(df_var))) +
      geom_point(aes(color = filtered), alpha = 0.7) +
      labs(title = "Median Count vs Variance",
           x = "rank Median Count",
           y = "Log2 Variance") +
      theme_minimal()
    
    # Create scatter plot for median count vs number of zeros
    plot_zeros <- ggplot(diagnostic_data, aes(x = rank(rowMeans(dataf[, -1], na.rm = TRUE)), y = df_zeros)) +
      geom_point(aes(color = filtered), alpha = 0.7) +
      labs(title = "Median Count vs Number of Zeros",
           x = "rank Median Count",
           y = "Number of Zeros") +
      theme_minimal()
    
    # Return a list of plots
    return(list(plot_variance, plot_zeros))
  }
  
  # Display the plots
  output$plot_variance <- renderPlot({
    plots <- counts_scatter_plot(input$counts_variance, input$counts_non_zero)
    print(plots[[1]])
  })
  
  output$plot_zeros <- renderPlot({
    plots <- counts_scatter_plot(input$counts_variance, input$counts_non_zero)
    print(plots[[2]])
  })
  
  
  
  
  
  
  
}



###################################################
#                 RUN APP                         #
###################################################


# Run the application 
shinyApp(ui = ui, server = server)

