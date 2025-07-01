# app.R
library(shiny)
library(glue)
library(DT)
library(reticulate)
library(shinydashboard)
library(promises)
library(future)
library(shinyWidgets)
library(spsComps)

plan(multisession)

# ---- Python Imports ----
os        <- import("os")
glob      <- import("glob")
re        <- import("re")
subprocess<- import("subprocess")
argparse  <- import("argparse")
multi     <- import("multiprocessing")
shutil    <- import("shutil")

# Ensure Python path
reticulate::py_run_string(
  "import sys; import os; sys.path.append(f'{os.getcwd()}')"
)

Sys.setenv(PYTHONUNBUFFERED = "1")


options(shiny.maxRequestSize = 100*1024^2)


# ---- UI ----
ui <- dashboardPage(
  dashboardHeader(
    title      = HTML("<em>Pyricularia oryzae</em><br>MiniMonsterPlex Tool"),
    titleWidth = 350
  ),
  dashboardSidebar(
    width = 350,
    sidebarMenu(
      id = "tabs",
      menuItem("Inputs",            tabName = "inputs",  icon = icon("upload")),
      menuItem("Alignment Summary", tabName = "summary", icon = icon("table")),
      menuItem("Phylogenetic Tree", tabName = "tree",    icon = icon("tree", lib = "font-awesome"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .form-group { margin-bottom: 15px; }
      .btn-primary { background-color: #007bff; border-color: #007bff; }
      .btn-primary:hover { background-color: #0056b3; border-color: #004085; }
    "))),
    tabItems(
      # ---- Inputs tab ----
      tabItem(tabName = "inputs",
              fluidRow(
                # Left: placeholder for uploads
                box(
                  width = 6, title = "Upload Data", status = "primary", solidHeader = TRUE,
                  uiOutput("file_upload_ui")
                ),
                # Right: project selector + metadata filename
                box(
                  width = 6, title = "Settings", status = "primary", solidHeader = TRUE,
                  selectizeInput(
                    "project_id", "Project ID",
                    choices = NULL, selected = "",
                    options = list(create = TRUE, placeholder = "Enter Project ID")
                  ),
                  textInput("metadata.file", "Metadata filename", value = "metadata.csv")
                )
              ),
              fluidRow(
                column(12, align = "center",
                       actionButton(
                         "myButton", "Run MiniMonsterPlex",
                         class = "btn-primary",
                         style = "width:250px; margin-top:20px;"
                       )
                )
              )
      ),
      
      # ---- Alignment Summary tab ----
      tabItem(tabName = "summary",
              fluidRow(
                box(
                  width = 12, title = "Alignment Summary Table", status = "info", solidHeader = TRUE,
                  DTOutput("alignment_table")
                )
              )
      ),
      
      # ---- Phylogenetic Tree tab ----
      tabItem(tabName = "tree",
              fluidRow(
                box(
                  width = 12, title = "Phylogenetic Tree", status = "info", solidHeader = TRUE,
                  uiOutput("phylogenetic_tree_display")
                )
              )
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  
  
  setwd(getwd())
  base_dir     <- getwd()
  projects_dir <- file.path(base_dir, "Projects")
  dir.create(projects_dir, showWarnings = FALSE)
  # Create www directory for serving files like the PDF tree
  dir.create("www", showWarnings = FALSE)
  
  # Reactive value to store the final vector of selected files
  final_file_vector <- reactiveVal(character(0))
  # Reactive value to store the URL of the current tree PDF
  current_tree_url <- reactiveVal(NULL)
  
  # Create necessary directories as soon as project ID is selected
  observeEvent(input$project_id, {
    req(nzchar(input$project_id))  # Ensure project_id is not empty
    
    project_id      <- input$project_id
    project_dir     <- file.path(projects_dir, project_id)
    fastq_dir       <- file.path(project_dir, "newFastq")
    output_dir      <- file.path(project_dir, "output")
    metadata_dir    <- file.path(project_dir, "metadata")
    
    dir.create(fastq_dir,    recursive = TRUE, showWarnings = FALSE)
    dir.create(output_dir,   recursive = TRUE, showWarnings = FALSE)
    dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
    
    message("Project directories created for: ", project_id)
    
    # Check for an existing tree PDF and update the display if it's found
    tree_filename <- paste0(project_id, "_tree.pdf")
    source_pdf_path <- file.path(output_dir, "tree_out", tree_filename)
    
    if (file.exists(source_pdf_path)) {
      dest_pdf_path <- file.path("www", tree_filename)
      file.copy(source_pdf_path, dest_pdf_path, overwrite = TRUE)
      current_tree_url(tree_filename)
    } else {
      current_tree_url(NULL)
    }
  })
  
  
  # Populate project dropdown
  observe({
    existing <- list.dirs(projects_dir, full.names = FALSE, recursive = FALSE)
    updateSelectizeInput(session, "project_id",
                         choices = existing,
                         selected = "",
                         server = TRUE)
  })
  
  # Dynamically render fileInputs only when project_id is non-empty
  output$file_upload_ui <- renderUI({
    req(nzchar(input$project_id))
    tagList(
      fileInput("fastq_files",     "Upload FASTQ files:", multiple = TRUE,
                accept = c(".gz", ".fastq", ".fq")),
      fileInput("metadata_upload", "Upload metadata file", accept = ".csv")
    )
  })
  
  # Copy FASTQ files
  observe({
    req(nzchar(input$project_id), input$fastq_files)
    fastq_dir <- file.path(projects_dir, input$project_id, "newFastq")
    dir.create(fastq_dir, recursive = TRUE, showWarnings = FALSE)
    for (i in seq_along(input$fastq_files$datapath)) {
      file.copy(
        input$fastq_files$datapath[i],
        file.path(fastq_dir, input$fastq_files$name[i]),
        overwrite = TRUE
      )
    }
    message("FASTQ files copied to: ", fastq_dir)
  })
  
  # Copy metadata file
  observe({
    req(nzchar(input$project_id), input$metadata_upload)
    meta_dir <- file.path(projects_dir, input$project_id, "metadata")
    dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(
      input$metadata_upload$datapath,
      file.path(meta_dir, input$metadata_upload$name),
      overwrite = TRUE
    )
    message("Metadata copied to: ", file.path(meta_dir, input$metadata_upload$name))
  })
  
  # Run analysis
  observeEvent(input$myButton, {
    withProgress( 
      message = 'MiniMonsterPlex is running...',
      value = 0,
      min = 0,
      max = 10,
      {
        req(nzchar(input$project_id), input$metadata.file)
        
        project_id <<- input$project_id
        
        setProgress(1, message = 'Sourcing Python script...')
        reticulate::source_python("MiniMonsterPlex_shiny.py")
        
        setProgress(3, message = 'Starting main analysis...')
        
        tryCatch({
          # Call the main python function
          file_name_list <- main(project_id, input$metadata.file)
          
          setProgress(4, message = 'Analysis complete. Displaying filter selection.')
          
          # Show the modal dialog with checkboxes for file selection
          showModal(modalDialog(
            title = "Select Files",
            p("Select the files you want to include in the final anaylsis"),
            # Checkbox group with the list of files from the python script
            checkboxGroupInput("selected_files_checkbox", "Returned Files:",
                               choices = file_name_list,
                               selected = file_name_list), # All selected by default
            footer = tagList(
              # This button confirms the selection and triggers the next observer
              actionButton("confirm_selection", "Create Vector", class = "btn-primary"),
              modalButton("Cancel")
            ),
            easyClose = TRUE
          ))
          
          setProgress(6, message = 'Awaiting file selection...')
          
        },
        error = function(e) {
          # This block now catches the error raised from Python
          error_message <- e$message
          
          # Clean up the message from reticulate if needed
          error_message <- gsub(".*RuntimeError: ", "", error_message)
          
          # Show the detailed error in a modal dialog
          showModal(modalDialog(
            title = "A Python Error Has Occurred",
            tags$h4("The analysis failed with the following error:"),
            # Use a <pre> tag to preserve formatting of the error message
            tags$pre(style = "white-space: pre-wrap; word-wrap: break-word; color: #a94442; background-color: #f2dede; border: 1px solid #ebccd1; padding: 15px; border-radius: 4px;", 
                     error_message),
            easyClose = TRUE,
            footer = NULL
          ))
        }
        )
      })
    
    # Trigger table refresh regardless of success or failure
    output$alignment_table <- renderDT({
      req(nzchar(input$project_id))
      
      project_id <- input$project_id
      summary_path <- file.path(projects_dir, project_id, "output/alignment_summary.csv")
      
      if (file.exists(summary_path)) {
        df <- read.csv(summary_path)
        colnames(df) <- c("Sample ID", "# Reads", "# Aligned", "Percent")
        # Remove duplicates based on "sample ID", keeping the latest occurrence
        df <- df[!duplicated(df[["Sample ID"]], fromLast = TRUE), ]
        df <- df[order(df[["Sample ID"]]), ]
        datatable(
          df, extensions = 'Buttons', rownames = FALSE,
          options = list(dom='Bfrtip', buttons=c('copy','csv','excel','pdf','print'))
        )
      } else {
        datatable(
          data.frame(Message = "No alignment summary found. Run the analysis or check for errors."),
          rownames = FALSE
        )
      }
    })
  })
  
  # Observer for the modal dialog confirmation button
  observeEvent(input$confirm_selection, {
    # Store the selected files from the checkbox input into the reactive value
    final_file_vector(input$selected_files_checkbox)
    
    # Close the pop-up window
    removeModal()
    
    filtered <- final_file_vector()
    
    reticulate::source_python("MiniMonsterPlex_shiny_raxml.py")
    
    setProgress(8, message = 'Building Tree')
    tryCatch({
      # Call the main python function
      main(project_id, filtered)
      
      # After the python script runs, find the PDF and prepare it for display
      tree_filename <- paste0(project_id, "_tree.pdf")
      source_pdf_path <- file.path(projects_dir, project_id, "output", "tree_out", tree_filename)
      
      if (file.exists(source_pdf_path)) {
        dest_pdf_path <- file.path("www", tree_filename)
        file.copy(source_pdf_path, dest_pdf_path, overwrite = TRUE)
        current_tree_url(tree_filename) # Update the reactive value with the file name
        updateTabItems(session, "tabs", "tree") # Switch to the tree tab
      } else {
        current_tree_url(NULL) # Set to NULL if tree not found
      }
      
      setProgress(10, message = 'Finished')
      # Show a final success message to the user
      showModal(modalDialog(
        title = "Success!",
        paste("Analysis for project", input$project_id, "completed."),
        easyClose = TRUE
      ))
      
    },
    error = function(e) {
      # This block now catches the error raised from Python
      error_message <- e$message
      
      # Clean up the message from reticulate if needed
      error_message <- gsub(".*RuntimeError: ", "", error_message)
      
      # Show the detailed error in a modal dialog
      showModal(modalDialog(
        title = "A Python Error Has Occurred",
        tags$h4("The analysis failed with the following error:"),
        # Use a <pre> tag to preserve formatting of the error message
        tags$pre(style = "white-space: pre-wrap; word-wrap: break-word; color: #a94442; background-color: #f2dede; border: 1px solid #ebccd1; padding: 15px; border-radius: 4px;", 
                 error_message),
        easyClose = TRUE,
        footer = NULL
      ))
    }
    )
  })
  
  # Initial render of the alignment summary table
  output$alignment_table <- renderDT({
    datatable(
      data.frame(Message = "Select a project and run the analysis to see results here."),
      rownames = FALSE
    )
  })
  
  # Render for the phylogenetic tree UI. This will show the PDF or a message.
  output$phylogenetic_tree_display <- renderUI({
    url <- current_tree_url()
    if (!is.null(url)) {
      # Use an iframe to embed the PDF from the 'www' directory
      tags$iframe(style = "height: 800px; width: 100%; border: none;", src = url)
    } else {
      p("Phylogenetic tree not available. Please run the full analysis or select a project with a completed analysis.")
    }
  })
}

shinyApp(ui, server)
