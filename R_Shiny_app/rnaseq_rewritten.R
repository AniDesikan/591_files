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

ui <- fluidPage(
  titlePanel("Plotting of Normalized Counts Data"),
)
  
  
# Run the application 
shinyApp(ui = ui, server = server)

