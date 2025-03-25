#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(glue)

#server logic to run mini monsterplex

runMiniMonsterPlex <- function(output.folder,metadata.file,input.folder,isolate.list,
                               isolate.file,host.list,host.file){
  # If the isolate list is empty, set it to NULL; otherwise, split it into a list
  if (isolate.list == "") {
  	isolate.list = NULL
  } else {
  	isolate.list = unlist(strsplit(isolate.list, " "))
  }
  
  # If the isolate file is empty, set it to NULL
  if (isolate.file == "") {
    isolate.file = NULL
  }
  
  # If the host list is empty, set it to NULL; otherwise, split it into a list
  if (host.list == "") {
    host.list = NULL
  } else {
    host.list = unlist(strsplit(host.list, " "))
  }
  
  # If the host file is empty, set it to NULL
  if (host.file == "") {
    host.file = NULL
  }
  
  # Construct the command options based on the provided isolate and host lists/files
  if (is.null(isolate.list)) {
    isolate.list.command = " "
  } else {
    isolate.list.command = glue(' -i {isolate.list}')
  }
  if (is.null(isolate.file)) {
    isolate.file.command = " "
  } else {
    isolate.file.command = glue(' -il {isolate.file}')
  }
  if (is.null(host.list)) {
    host.list.command = " "
  } else {
    host.list.command = glue(' -hf {host.list}')
  }
  if (is.null(host.file)) {
    host.file.command = " "
  } else {
    host.file.command = glue(' -hfl {host.file}')
  }
  
  # Construct the final command to run the MiniMonsterPlex.py script with the provided options
  command = glue('python MiniMonsterPlex.py -o {output.folder} -m {metadata.file} -f {input.folder} {isolate.list.command} {isolate.file.command} {host.list.command} {host.file.command}')
  
  # Execute the constructed command
  print(command)
  
}



# Define server logic required to draw a histogram
function(input, output, session) {
    observeEvent(input$myButton, {
      runMiniMonsterPlex(input$output.folder,input$metadata.file,input$input.folder,input$isolate.list,input$isolate.file,input$host.list,input$host.file)
    })
    output$distPlot <- renderPlot({

        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')

    })

}
