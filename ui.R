#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

# Define UI for application that draws a histogram
fluidPage(

    # Application title
    titlePanel("MiniMonsterPlex"),

    textInput("output.folder", "Enter output folder here", value = "output",placeholder = "Must not already exist"),
    textInput("metadata.file", "Enter metadata_file here", value = "metadata.csv"),
    textInput("input.folder", "Enter input folder here", value = "fastq"),
    textInput("isolate.list", "Enter a space seperated list of isolates you want included here", value = ""),
    textInput("isolate.file", "Enter a file containing a new line list of isolates you want included here", value = ""),
    textInput("host.list", "Enter a space seperated list of hosts you want included here", value = ""),
    textInput("host.file", "Enter a file containing a new line list of hosts you want included here", value = ""),
    actionButton("myButton", "Run MiniMonsterPlex")
)
