library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(shinycssloaders)
library(shinyjs)
library(reticulate)
library(DT)
library(readr)
library(shinyBS)
library(ggpubr)
library(shinyWidgets)

# Load optional packages with graceful fallback
if (!requireNamespace("fenr", quietly = TRUE)) {
  warning("Package 'fenr' not available. Some enrichment functions may not work.")
} else {
  library(fenr)
}

if (!requireNamespace("shinydisconnect", quietly = TRUE)) {
  warning("Package 'shinydisconnect' not available. Disconnect handling disabled.")
} else {
  library(shinydisconnect)
}

library(stringr)
library(jsonlite)

# Configure Python environment for reticulate
tryCatch({
  # Check if we're in a conda environment
  if (Sys.getenv("CONDA_DEFAULT_ENV") != "") {
    # Use conda environment
    reticulate::use_condaenv("fibrosis_shiny", required = TRUE)
  } else {
    # Fallback to virtualenv
    reticulate::use_virtualenv("fibrosis_shiny")
  }
  
  # Import Python modules
  sc <- reticulate::import("scanpy")
  # Configure scanpy figures for high resolution export (600 DPI)
  # This significantly improves the quality of images shown in the app and saved by users
  sc$set_figure_params(dpi = 100, dpi_save = 600, format = 'png')
  dc <- reticulate::import("decoupler")
  pydeseq2_dds <- reticulate::import("pydeseq2.dds")
  pydeseq2_ds <- reticulate::import("pydeseq2.ds")
  
  cat("✅ Python environment and packages loaded successfully\n")
}, error = function(e) {
  cat("⚠️ Python environment setup failed:", e$message, "\n")
  cat("📦 Application will run with limited Python functionality\n")
  
  # Create placeholder objects to prevent errors
  sc <- NULL
  dc <- NULL
  pydeseq2_dds <- NULL
  pydeseq2_ds <- NULL
})

# Load datasets configuration  
datasets_config <- jsonlite::fromJSON("config/datasets_config.json")

# Define NULL-coalescing operator for cleaner code
`%||%` <- function(x, y) if (is.null(x)) y else x

# Helper function to get organism choices from config
get_organism_choices <- function(config) {
  choices <- list()
  for (organism in names(config)) {
    organism_data <- config[[organism]]
    status <- organism_data[["Status"]] %||% "Unknown"
    description <- organism_data[["Description"]] %||% "No description"
    
    # Create a display name with status indicator
    if (status == "Available") {
      display_name <- paste0("✅ ", organism)
    } else if (status == "Optional") {
      display_name <- paste0("⚪ ", organism, " (Optional)")
    } else {
      display_name <- paste0("⏳ ", organism, " (", status, ")")
    }
    
    choices[[display_name]] <- organism
  }
  return(choices)
}

# Helper function to get organism data safely
get_organism_data <- function(config, organism_key) {
  if (!organism_key %in% names(config)) {
    return(list(
      datasets = character(0),
      status = "Not Found",
      description = "Organism not found in configuration",
      message = NULL
    ))
  }
  
  organism_data <- config[[organism_key]]
  return(list(
    datasets = organism_data[["Datasets"]] %||% character(0),
    status = organism_data[["Status"]] %||% "Unknown",
    description = organism_data[["Description"]] %||% "No description available",
    message = organism_data[["Message"]] %||% NULL,
    download_command = organism_data[["DownloadCommand"]] %||% NULL
  ))
}

# Function to reload datasets configuration if needed
reload_datasets_config <- function() {
  tryCatch({
    jsonlite::fromJSON("config/datasets_config.json")
  }, error = function(e) {
    warning("Failed to reload datasets configuration: ", e$message)
    return(datasets_config) # Return original config if reload fails
  })
}

# Function to validate dataset file existence
validate_dataset_path <- function(organism, dataset_id) {
  base_path <- paste0("datasets/", organism, "/", dataset_id, ".h5ad")
  file_exists <- file.exists(base_path)
  file_size <- if (file_exists) file.size(base_path) else 0
  
  return(list(
    path = base_path,
    exists = file_exists,
    size = file_size
  ))
}

# Generate organism choices from config
organism_choices <- get_organism_choices(datasets_config)

ui <- fluidPage(
  # Add disconnect message only if shinydisconnect is available
  if (requireNamespace("shinydisconnect", quietly = TRUE)) {
    shinydisconnect::disconnectMessage(
      text = "An error occurred. Please refresh the page and try again.",
      refresh = "Refresh",
      background = "#FFFFFF",
      colour = "#444444",
      refreshColour = "#2c3e50",
      overlayColour = "#000000",
      overlayOpacity = 0.6,
      width = 450,
      top = 50,
      size = 22,
      css = ""
    )
  },
  navbarPage(
    title = "Multi-species scRNA-seq Atlas of MASLD",
    theme = bs_theme(preset = 'flatly', base_font = 'Lato', code_font = 'Lato', heading_font = 'Lato'),
    fluid = TRUE,
    windowTitle = "MASLDatlas",
    id = "navbar",
    selected = "tab_explore_datasets",
    tags$head(
      tags$link(rel = "icon", type = "image/png", sizes = "32x32", href = "tabicon.PNG"),
      tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css"),
      tags$link(rel = "stylesheet", type = "text/css", href = "custom-modern.css"),
      useShinyjs(),
      tags$style(HTML("
  /* Variables CSS pour les couleurs */
  :root {
    --primary-color: #2c3e50;
    --secondary-color: #3498db;
  }

  /* Navbar full width et améliorée */
  .navbar {
    background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%) !important;
    border: none !important;
    width: 100vw !important;
    left: 0 !important;
    margin-left: calc(-50vw + 50%) !important;
    margin-right: calc(-50vw + 50%) !important;
    position: relative !important;
    z-index: 1000 !important;
  }

  /* Force navbar container to break out of parent containers */
  .navbar > .container-fluid {
    width: 100vw !important;
    max-width: 100vw !important;
    padding-left: 15px !important;
    padding-right: 15px !important;
    margin: 0 !important;
  }

  /* Remove any padding/margin from body that could affect navbar */
  body .container-fluid {
    padding-left: 0 !important;
    padding-right: 0 !important;
  }

  /* Specific override for navbar container only */
  body > .container-fluid {
    padding: 0 !important;
    margin: 0 !important;
  }

  .navbar-brand {
    font-weight: bold !important;
    color: white !important;
  }

  .navbar-nav .nav-link {
    color: rgba(255,255,255,0.9) !important;
    font-weight: 500 !important;
  }

  .navbar-nav .nav-link:hover,
  .navbar-nav .nav-link.active {
    color: white !important;
    background-color: rgba(255,255,255,0.1) !important;
    border-radius: 4px !important;
  }

  /* Cache le workflow tab */
  li[data-value='tab_workflow'] {
    display: none !important;
  }

  .nav-link[data-value='tab_workflow'] {
    display: none !important;
  }

  /* Container full width */
  .container-fluid {
    width: 100% !important;
    max-width: 100% !important;
    padding-left: 15px !important;
    padding-right: 15px !important;
    margin: 0 !important;
  }

  /* Force body and html to use full width */
  body, html {
    width: 100% !important;
    margin: 0 !important;
    padding: 0 !important;
  }

  /* Bootstrap container overrides */
  .container, .container-fluid {
    width: 100% !important;
    max-width: 100% !important;
    padding-left: 15px !important;
    padding-right: 15px !important;
  }

  /* Tab content full width */
  .tab-content {
    width: 100% !important;
    padding: 0 !important;
  }

  .tab-pane {
    width: 100% !important;
    padding: 15px !important;
  }

  /* Navbar container full width */
  .navbar > .container-fluid {
    width: 100% !important;
    max-width: 100% !important;
    padding-left: 15px !important;
    padding-right: 15px !important;
  }

  #imageoutput_violin_expression_clusters img,
  #imageoutput_violin_expression_clusters_groups img,
  #imageoutput_dge_ranks img,
  #imageoutput_dge_violin img,
  #de_enrichment_table,
  #pca_pseudo_bulk img,
  #collectri_top img,
  #collectri_volcano img,
  #collectri_network img,
  #progeny_top img,
  #progeny_targets img,
  #msigdb_top img,
  #msigdb_running img{
    width: 100% !important;
    height: auto !important;
    object-fit: contain;
    display: block;
  }
  
  @media (min-width: 1300px) and (max-width: 1600px) {
    .umap-left.col-sm-6 {
      width: 60%;
      float: left;
      padding-left: 0;
      padding-right: 0;
    }
    .umap-right.col-sm-6 {
      width: 40%;
      float: left;
      padding-left: 0;
      padding-right: 0;
    }
  }
  
  @media (max-width: 1299px) {
    .umap-left.col-sm-6 {
      width: 70%;
      float: left;
      padding-left: 0;
      padding-right: 0;
    }
    .umap-right.col-sm-6 {
      width: 30%;
      float: left;
      padding-left: 0;
      padding-right: 0;
    }
  }
"))
    ),
    
    # Premier onglet - Explore & Analyze Datasets (par défaut)
    tabPanel(
      title = div(
        span("📊", style = "margin-right: 5px;"),
        "Explore & Analyze Datasets"
      ),
      value = "tab_explore_datasets",      
      # Ajout d'une description de l'onglet
      div(class = "tab-description", style = "background: #f8f9fa; padding: 15px; margin: 10px 0; border-radius: 8px; border-left: 4px solid #3498db;",
        h4("🔬 Single-cell RNA-seq Data Analysis Pipeline", style = "margin-top: 0; color: #2c3e50;"),
        p("Explore multi-species scRNA-seq datasets, perform differential expression analysis, and visualize cellular heterogeneity. Follow the step-by-step workflow below:")
      ),
      
      tabsetPanel(
        type = "pills",
        id = "main_analysis_tabs",
        
        # Step 1: Import Dataset
        tabPanel(
          title = div(
            span("📥", style = "margin-right: 5px;"),
            "1. Import Dataset"
          ),
          value = "subtab_import_dataset",
          
          div(class = "step-header", style = "background: linear-gradient(135deg, #e8f5e8, #f0f8ff); padding: 10px; margin-bottom: 15px; border-radius: 6px;",
            h5("Step 1: Load and Visualize Your Dataset", style = "margin: 0; color: #27ae60;"),
            p("Select an organism and dataset to begin your analysis. View UMAP projections and explore cellular composition.", 
              style = "margin: 5px 0 0 0; font-size: 0.9rem; color: #666;")
          ),
          
          hr(),
          sidebarLayout(
            sidebarPanel(
              width = 2,
              class = "analysis-sidebar",
              
              div(class = "sidebar-section",
                selectInput("selection_organism", "Select Organism",
                           choices = organism_choices)
              ),
              
              div(class = "sidebar-section",
                uiOutput("dataset_selection_list")
              ),
              
              # Add dataset size selection for large datasets
              conditionalPanel(
                condition = "input.selection_dataset && input.selection_dataset.includes('Fibrotic')",
                div(class = "sidebar-section",
                  selectInput("dataset_size_option", "Dataset Size",
                             choices = c(
                               "Full dataset (9.2GB)" = "full",
                               "Large subset (20k cells)" = "sub20k", 
                               "Medium subset (10k cells)" = "sub10k",
                               "Small subset (5k cells)" = "sub5k"
                             ),
                             selected = "sub10k"),
                  div(class = "help-text", 
                      style = "font-size: 11px; color: #666; background: #fff3cd; padding: 8px; border-radius: 4px; margin-top: 8px;",
                      "⚠️ Full dataset may take 30+ minutes to load")
                )
              ),
              
              div(class = "sidebar-action",
                actionButton("import_dataset", 
                           class = "btn-primary btn-block", 
                           "Load Dataset", 
                           width = '100%')
              )
            ),
            mainPanel(width = 8,
              fluidRow(
                column(width = 6, textOutput("textoutput_data_numbers")),
                column(width = 6, uiOutput("umap_col_selection"))
              ),
              fluidRow(
                column(width = 6, class = "umap-left",
                       shinycssloaders::withSpinner(
                         imageOutput("imageoutput_UMAP", width = "100%", height = "auto") #EREN
                       )
                       ),
                column(width = 6, class = "umap-right", 
                       div(
                         style = "width: 100%; height: 100%; padding: 0; overflow: visible;",
                         shinycssloaders::withSpinner(
                           imageOutput("imageoutput_UMAP_selection", width = "100%", height = "auto")
                         )
                )
                )
              ),
              br(),
              hr(),
              br(),
              fluidRow(
                column(width = 6, uiOutput("ranks_col_selection")),
                column(width = 6)
              ),
              fluidRow(
                column(width = 6, 
                       downloadButton("download_markers", "📥 Download Markers (CSV)", 
                                    icon = icon("download"), class = "btn-primary btn-sm mb-2"),
                       br(),
                       shinycssloaders::withSpinner(DTOutput("table_cluster_markers"))
                ),
                column(width = 6, shinycssloaders::withSpinner(imageOutput("imageoutput_CellType_groups")))
              )
            )
          )
        ),
        # Step 2: Cluster Selection
        tabPanel(
          title = div(
            span("🎯", style = "margin-right: 5px;"),
            "2. Cluster Selection"
          ),
          value = "subtab_cluster_selection",
          div(class = "step-header", style = "background: linear-gradient(135deg, #fff3e0, #e8f5e8); padding: 10px; margin-bottom: 15px; border-radius: 6px;",
            h5("Step 2: Select and Analyze Cell Clusters", style = "margin: 0; color: #f39c12;"),
            p("Choose specific cell clusters for detailed analysis. Visualize gene expression patterns and co-expression relationships.", 
              style = "margin: 5px 0 0 0; font-size: 0.9rem; color: #666;")
          ),
          
          hr(),
          sidebarLayout(
            sidebarPanel(width = 2,
                         bsTooltip(id = "selected_culusters",
                                   title = "This is your tooltip text",
                                   placement = "right",
                                   trigger = "hover"),
                         uiOutput("cluster_selection_list"),
                         uiOutput("button_cluster_filter"),
                         hr(),
                         selectInput("cluster_selection_visualization_type", "Select Visualization Method",
                                     choices = c("Visualize Expression of Gene",
                                                 "Visualize Expression of Geneset",
                                                 "Calculate Co-Expression")),
                         hr(),
                         conditionalPanel(condition = "input.cluster_selection_visualization_type == 'Visualize Expression of Gene'",
                                          uiOutput("gene_list_expression")),
                         conditionalPanel(condition = "input.cluster_selection_visualization_type == 'Visualize Expression of Geneset'",
                                          textInput("name_first_geneset", "Name First Geneset", placeholder = "Enter Name", width = '100%'),
                                          textAreaInput("first_geneset_list", "Enter Gene Names", placeholder = "Enter Names", width = '100%', height = '100%')

                                          ),
                         
                         uiOutput("gene_list_coexpresion_first"),
                         uiOutput("gene_list_coexpresion_second"),
                         uiOutput("button_visualize_cluster_selection"),
                         conditionalPanel(condition = "input.cluster_selection_visualization_type == 'Visualize Expression of Geneset'",
                                          br(),
                                          hr(),
                                          br(),
                                          uiOutput("geneset_name_input_second"),
                                          uiOutput("geneset_list_input_second"),
                                          uiOutput("geneset_list_go_second"))
            ),
            mainPanel(
              conditionalPanel(
                condition = "input.cluster_selection_visualization_type == 'Visualize Expression of Gene'",
                fluidRow(
                  column(width = 6, textOutput("info_cluster_selection")),
                  column(width = 6, uiOutput("umap_col_selection_clusters"))
                ),
                fluidRow(
                  column(width = 6, class = "umap-left", shinycssloaders::withSpinner(imageOutput("imageoutput_UMAP_expression_clusters"))),
                  column(width = 6, class = "umap-right", shinycssloaders::withSpinner(imageOutput("imageoutput_UMAP_selection_clusters")))
                ),
                br(),
                hr(),
                br(),
                fluidRow(
                  column(width = 6, 
                         div(
                          style = "width: 100%; height: 100%; padding: 0;",
                         shinycssloaders::withSpinner(imageOutput("imageoutput_violin_expression_clusters_groups", width = "100%", height = "100%"))
                         )
                         ),
                  column(width = 6, 
                         div(
                           style = "width: 100%; height: 100%; padding: 0;",
                         shinycssloaders::withSpinner(imageOutput("imageoutput_violin_expression_clusters", width = "100%", height = "100%"))
                         )
                         )
                )
              ),
              conditionalPanel(condition = "input.cluster_selection_visualization_type == 'Visualize Expression of Geneset'  && input.visualize_cluster_selection != 0",
                               fluidRow(column(width = 6), column(width = 6,
                                                                  selectInput("first_geneset_enrichment_violin_selection", label = NULL, choices = c("Clusters", "Groups")))),
                               fluidRow(column(width = 6,
                                               div(
                                                style = "width: 100%; height: 100%; padding: 0;",
                                               shinycssloaders::withSpinner(imageOutput("image_output_enrichment_first_set", width = "100%", height = "100%")))
                                        ), 
                                        column(width = 6,
                                               div(
                                                style = "width: 100%; height: 100%; padding: 0;",
                                               conditionalPanel(condition = "input.first_geneset_enrichment_violin_selection == 'Clusters'",
                                                                shinycssloaders::withSpinner(imageOutput("image_output_enrichment_first_violin_clusters", width = "100%", height = "100%"))
                                                                ),
                                               conditionalPanel(condition = "input.first_geneset_enrichment_violin_selection == 'Groups'",
                                                                shinycssloaders::withSpinner(imageOutput("image_output_enrichment_first_violin_groups", width = "100%", height = "100%"))
                                                                )
                                               )
                               )
                                               ),
                               br(),
                               hr(),
                               br(),
                               fluidRow(column(width = 6), column(width = 6,
                                                                  conditionalPanel(condition = "input.second_geneset_run != 0",
                                                                                   selectInput("second_geneset_enrichment_violin_selection", label = NULL, choices = c("Clusters", "Groups")))
                                                                  
                                                                  
                                                                  )),
                               fluidRow(column(width = 6,
                                               shinycssloaders::withSpinner(imageOutput("image_output_enrichment_second_set"))
                                               ), column(width = 6,
                                                         conditionalPanel(condition = "input.second_geneset_enrichment_violin_selection == 'Clusters'",
                                                                          shinycssloaders::withSpinner(imageOutput("image_output_enrichment_second_violin_clusters"))
                                                         ),
                                                         conditionalPanel(condition = "input.second_geneset_enrichment_violin_selection == 'Groups'",
                                                                          shinycssloaders::withSpinner(imageOutput("image_output_enrichment_second_violin_groups"))
                                                         )
                                                         
                                                         )))

              ,
              conditionalPanel(
                condition = "input.cluster_selection_visualization_type == 'Calculate Co-Expression' && input.visualize_cluster_selection != 0",
                fluidRow(
                  column(width = 6, shinycssloaders::withSpinner(imageOutput("imageoutput_UMAP_coexpression_first"))),
                  column(width = 6, shinycssloaders::withSpinner(imageOutput("imageoutput_UMAP_coexpression_second")))
                ),
                br(),
                hr(),
                br(),
                fluidRow(
                  column(width = 6,
                         radioGroupButtons(
                           inputId = "test_choice",
                           label = NULL,
                           choices = c("Spearman", "Pearson"),
                           justified = TRUE, width = '100%'
                         ),
                         bslib::input_switch("remove_zero_counts", "Remove Zero Counts", width = '100%')
                  ),
                  column(width = 6, uiOutput("top_correlated_genes"))
                ),
                fluidRow(
                  column(width = 6, shinycssloaders::withSpinner(plotOutput("stats_cor_plot"))),
                  column(width = 6,
                         downloadButton("download_correlation", "📥 Download Correlation Results (CSV)", 
                                      icon = icon("download"), class = "btn-primary btn-sm mb-2"),
                         br(),
                         fluidRow(
                           column(width = 6, shinycssloaders::withSpinner(DTOutput("first_gene_correlation_table"))),
                           column(width = 6, shinycssloaders::withSpinner(DTOutput("second_gene_correlation_table")))
                         )
                  )
                )
              )
            )
          )
        ),
        # Step 3: Differential Expression
        tabPanel(
          title = div(
            span("📊", style = "margin-right: 5px;"),
            "3. Differential Expression"
          ),
          value = "subtab_differential_expression",          
          div(class = "step-header", style = "background: linear-gradient(135deg, #fce4ec, #f3e5f5); padding: 10px; margin-bottom: 15px; border-radius: 6px;",
            h5("Step 3: Compare Gene Expression Between Conditions", style = "margin: 0; color: #e91e63;"),
            p("Perform differential gene expression analysis between clusters or groups. Generate volcano plots and identify key biomarkers.", 
              style = "margin: 5px 0 0 0; font-size: 0.9rem; color: #666;")
          ),
          
          hr(),
          sidebarLayout(
            sidebarPanel(width = 2,
                         radioGroupButtons(
                           inputId = "de_data",
                           label = NULL,
                           choices = c("All Data", "Filtered Data"),
                           justified = TRUE, width = '100%'
                         ),
                         radioGroupButtons(
                           inputId = "de_type",
                           label = NULL,
                           choices = c("Clusters", "Groups"),
                           justified = TRUE, width = '100%'
                         ),
                         hr(),
                         uiOutput("de_ident_selection_name_first"),
                         uiOutput("de_ident_selection_first"),
                         uiOutput("de_ident_selection_name_second"),
                         uiOutput("de_ident_selection_second"),
                         uiOutput("de_method_selection"),
                         uiOutput("de_run_analysis")
            ),
            mainPanel(
              conditionalPanel(
                condition = "input.de_run_dge > 0",
                fluidRow(
                  column(width = 6, 
                         div(
                          style = "width: 100%; height: 100%; padding: 0;",
                         shinycssloaders::withSpinner(imageOutput("imageoutput_dge_ranks", width = "100%", height = "100%")))
                  ),
                  column(width = 6, 
                         downloadButton("download_dge", "📥 Download DGE Results (CSV)", 
                                      icon = icon("download"), class = "btn-primary btn-sm mb-2"),
                         br(),
                         uiOutput("dge_table_filter"),
                         shinycssloaders::withSpinner(DTOutput("dge_dt"))
                  )
                ),
                br(),
                hr(),
                br(),
                fluidRow(
                  column(width = 6, 
                         div(
                           style = "width: 100%; height: 100%; padding: 0;",
                         shinycssloaders::withSpinner(imageOutput("imageoutput_dge_violin", width = "100%", height = "100%")))
                  ),
                  column(width = 6,
                         radioGroupButtons(
                           inputId = "de_radiobutton",
                           label = NULL,
                           choices = c("Violin Plot", "Enrichment"),
                           justified = TRUE, width = '100%'
                         ),
                         conditionalPanel(
                           condition = "input.de_radiobutton == 'Violin Plot'",
                           div(
                             id = "violin_container",
                             tags$div(
                               id = "violin_fallback",
                               style = "display: none; padding: 20px; color: #666; text-align: center;",
                               "Please select a row above to generate the violin plot."
                             ),
                             shinycssloaders::withSpinner(imageOutput("violin_dge")),
                             tags$script(HTML("
      function monitorViolinImage() {
        var $fallback = $('#violin_fallback');
        var attempts = 0;
        var maxAttempts = 50; // Try for up to 5 seconds

        // Start by showing the message
        $fallback.show();

        var interval = setInterval(function() {
          var $img = $('#violin_container img');

          if ($img.length > 0 && $img.attr('src')) {
            // Listen for successful load
            $img.off('load error'); // prevent duplicate bindings

            $img.on('load', function() {
              if (this.naturalWidth > 0) {
                $fallback.hide();
                clearInterval(interval);
              }
            });

            $img.on('error', function() {
              $fallback.show();
              clearInterval(interval);
            });

            // If image is already loaded (quick response)
            if ($img[0].complete && $img[0].naturalWidth > 0) {
              $fallback.hide();
              clearInterval(interval);
            }
          }

          attempts++;
          if (attempts >= maxAttempts) {
            clearInterval(interval);
          }
        }, 100);
      }

      // Run the function on each reactive output change
      $(document).on('shiny:value', monitorViolinImage);
    "))
                           )
                         ),
                         conditionalPanel(
                           condition = "input.de_radiobutton == 'Enrichment'",
                           selectInput("de_enrichment_type", label = NULL, choices = c("GO", "BP", "KEGG", "Reactome", "WikiPathways")),
                           downloadButton("download_enrichment", "📥 Download Enrichment Results (CSV)", 
                                        icon = icon("download"), class = "btn-success btn-sm mb-2"),
                           br(),
                           div(
                             style = "width: 100%; height: 100%; padding: 0;", 
                           shinycssloaders::withSpinner(DTOutput("de_enrichment_table")))
                         )
                  )
                )
              )
            )
          )
        ),
        # Step 4: Pseudo Bulk Analysis
        tabPanel(
          title = div(
            span("📦", style = "margin-right: 5px;"),
            "4. Pseudo Bulk Analysis"
          ),
          value = "subtab_pseudo_bulk",          
          div(class = "step-header", style = "background: linear-gradient(135deg, #e0f2f1, #e8f5e8); padding: 10px; margin-bottom: 15px; border-radius: 6px;",
            h5("Step 4: Aggregate Single Cells for Bulk-like Analysis", style = "margin: 0; color: #00796b;"),
            p("Convert single-cell data to pseudo-bulk samples for robust differential expression analysis using DESeq2.", 
              style = "margin: 5px 0 0 0; font-size: 0.9rem; color: #666;")
          ),
          
          hr(),
          sidebarLayout(
            sidebarPanel(width = 2,
                         radioGroupButtons(inputId = "pseudo_bulk_data",
                                           label = NULL,
                                           choices = c("All Data", "Filtered Data"),
                                           justified = TRUE, width = '100%'
                         ),
                         uiOutput("pseudo_ident_selection_first"),
                         uiOutput("pseudo_ident_selection_second"),
                         uiOutput("run_pseudo_bulk"),
                         hr(),
                         uiOutput("pdata_clusters"),
                         uiOutput("run_deseq2_ui")
            ),
            mainPanel(
              conditionalPanel(
                condition = "input.run_pseudo_bulk_act > 0",
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("pseudo_pca_selection_ui")
                  ),
                  column(width = 6,
                         uiOutput("selection_pseudo_volcano_table_ui"))
                ),
                fluidRow(
                  column(width = 6, 
                         div(
                           style = "width: 100%; height: 100%; display: flex; align-items: center; justify-content: center; padding: 0;",
                         shinycssloaders::withSpinner(imageOutput("pca_pseudo_bulk", width = "100%", height = "100%")))),
                  column(width = 6,
                         conditionalPanel(condition = "input.selection_pseudo_volcano_table == 'Volcano'",
                                          div(
                                            style = "width: 100%; height: 100%; padding: 0;",
                                          shinycssloaders::withSpinner(imageOutput("pca_pseudo_bulk_volcano", width = "100%", height = "100%")))),
                         conditionalPanel(condition = "input.selection_pseudo_volcano_table == 'Result Table'",
                                          downloadButton("download_pseudobulk", "📥 Download Pseudo-bulk Results (CSV)", 
                                                       icon = icon("download"), class = "btn-primary btn-sm mb-2"),
                                          br(),
                                          shinycssloaders::withSpinner(DTOutput("pca_pseudo_bulk_results_table"))))
                ),
                br(),
                hr(),
                br(),
                fluidRow(
                  column(width = 6, uiOutput("pseudo_ora_selection_ui"),
                         
                         conditionalPanel(condition = "input.pseudo_ora_selection == 'CollecTRI'",
                                          uiOutput("collectri_outputs")),
                         conditionalPanel(condition = "input.pseudo_ora_selection == 'PROGENy'",
                                          uiOutput("progeny_outputs")),
                         conditionalPanel(condition = "input.pseudo_ora_selection == 'MSigDB Hallmark'",
                                          uiOutput("msigdb_outputs"))
                         ),
                  column(width = 6, uiOutput("enrichment_1_pseudo"))
                ),
                fluidRow(
                  column(width = 6,
                         conditionalPanel(condition = "input.pseudo_ora_selection == 'CollecTRI'",
                                          conditionalPanel(condition = "input.collectri_plot_selection == 'Top 25'",
                                                           div(
                                                             style = "width: 100%; height: 100%; display: flex; align-items: center; justify-content: center; padding: 0;",
                                                           shinycssloaders::withSpinner(imageOutput("collectri_top", width = "100%", height = "100%")))),
                                          conditionalPanel(condition = "input.collectri_plot_selection == 'Volcano'",
                                                           div(
                                                             style = "width: 100%; height: 100%; display: flex; align-items: center; justify-content: center; padding: 0;",
                                                           shinycssloaders::withSpinner(imageOutput("collectri_volcano", width = "100%", height = "100%")))),
                                          conditionalPanel(condition = "input.collectri_plot_selection == 'Network'",
                                                           div(
                                                             style = "overflow-x: auto; overflow-y: auto;",
                                                             shinycssloaders::withSpinner(imageOutput("collectri_network", width = "100%", height = "100%"))
                                                           ))),
                         conditionalPanel(condition = "input.pseudo_ora_selection == 'PROGENy'",
                                          conditionalPanel(condition = "input.progeny_plot_selection == 'Scores'",
                                                           shinycssloaders::withSpinner(imageOutput("progeny_top"))),
                                          conditionalPanel(condition = "input.progeny_plot_selection == 'Targets'",
                                                           shinycssloaders::withSpinner(imageOutput("progeny_targets")))
                                          
                                          ),
                         conditionalPanel(condition = "input.pseudo_ora_selection == 'MSigDB Hallmark'",
                                          conditionalPanel(condition = "input.msigdb_plot_selection == 'Top Signatures'",
                                                           shinycssloaders::withSpinner(imageOutput("msigdb_top"))),
                                          conditionalPanel(condition = "input.msigdb_plot_selection == 'Running Score'",
                                                           shinycssloaders::withSpinner(imageOutput("msigdb_running")))
                                          
                         )
                         ),
                  column(width = 6, 
                         downloadButton("download_pseudo_enrichment", "📥 Download Pseudo-bulk Enrichment (CSV)", 
                                      icon = icon("download"), class = "btn-success btn-sm mb-2"),
                         br(),
                         shinycssloaders::withSpinner(DTOutput("pseudo_enrichment_table")))
                )
              )
            )
          )
        )
        # ,
        # tabPanel(
        #   title = "Cell-Cell Interaction Analysis",
        #   value = "subtab_cell_cell_interaction_analysis_1",
        #   icon = icon("phone-volume")
        # )
      )
    ),
    
    # Documentation Tab
    tabPanel(
      title = div(
        span("📚", style = "margin-right: 5px;"),
        "Documentation"
      ),
      value = "tab_documentation",
      
      div(class = "container-fluid", style = "padding: 20px; max-width: 1200px; margin: 0 auto;",
        
        # Header
        div(class = "doc-header", style = "text-align: center; margin-bottom: 40px; padding: 30px; background: linear-gradient(135deg, #2c3e50 0%, #3498db 100%); color: white; border-radius: 10px;",
          h1("MASLDatlas Documentation", style = "margin: 0; font-size: 2.5em;"),
          p("Complete guide to using the Multi-species scRNA-seq Atlas of MASLD", style = "margin-top: 10px; font-size: 1.2em; opacity: 0.9;")
        ),
        
        # Table of Contents
        div(class = "toc-card", style = "background: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 30px; border-left: 4px solid #3498db;",
          h3("Table of Contents", style = "margin-top: 0; color: #2c3e50;"),
          tags$ul(style = "list-style-type: none; padding-left: 0;",
            tags$li(style = "margin: 10px 0;", 
              tags$a(href = "#overview", "1. Overview", style = "color: #3498db; text-decoration: none; font-weight: 500;")
            ),
            tags$li(style = "margin: 10px 0;", 
              tags$a(href = "#getting-started", "2. Getting Started", style = "color: #3498db; text-decoration: none; font-weight: 500;")
            ),
            tags$li(style = "margin: 10px 0;", 
              tags$a(href = "#workflow", "3. Analysis Workflow", style = "color: #3498db; text-decoration: none; font-weight: 500;")
            ),
            tags$li(style = "margin: 10px 0;", 
              tags$a(href = "#features", "4. Key Features", style = "color: #3498db; text-decoration: none; font-weight: 500;")
            ),
            tags$li(style = "margin: 10px 0;", 
              tags$a(href = "#export", "5. Exporting Results", style = "color: #3498db; text-decoration: none; font-weight: 500;")
            ),
            tags$li(style = "margin: 10px 0;", 
              tags$a(href = "#troubleshooting", "6. Troubleshooting", style = "color: #3498db; text-decoration: none; font-weight: 500;")
            )
          )
        ),
        
        # Section 1: Overview
        div(id = "overview", class = "doc-section", style = "margin-bottom: 40px;",
          h2("1. Overview", style = "color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px;"),
          div(style = "background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
            p("MASLDatlas is an interactive web application for exploring and analyzing single-cell RNA sequencing (scRNA-seq) data from multiple species affected by Metabolic dysfunction-Associated Steatotic Liver Disease (MASLD)."),
            h4("Key Capabilities:"),
            tags$ul(
              tags$li("Interactive visualization of scRNA-seq datasets (UMAP, Violin plots, Heatmaps)"),
              tags$li("Differential gene expression analysis"),
              tags$li("Gene correlation and co-expression studies"),
              tags$li("Functional enrichment analysis (GO, KEGG, Reactome, WikiPathways)"),
              tags$li("Pseudo-bulk DESeq2 analysis"),
              tags$li("CSV export functionality for all results")
            ),
            h4("Supported Species:"),
            tags$ul(
              tags$li("Human"),
              tags$li("Mouse"),
              tags$li("Zebrafish"),
              tags$li("Integrated datasets")
            )
          )
        ),
        
        # Section 2: Getting Started
        div(id = "getting-started", class = "doc-section", style = "margin-bottom: 40px;",
          h2("2. Getting Started", style = "color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px;"),
          div(style = "background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
            h4("Step 1: Select a Dataset"),
            p("Navigate to the", strong("Import Dataset"), "tab and:"),
            tags$ol(
              tags$li("Choose an organism from the dropdown menu"),
              tags$li("Select a dataset from the available options"),
              tags$li("Wait for the data to load (progress bar will appear)")
            ),
            
            h4("Step 2: Explore the Data", style = "margin-top: 20px;"),
            p("Once loaded, you can:"),
            tags$ul(
              tags$li("Visualize cell types on UMAP plots"),
              tags$li("Filter data by specific cell types or conditions"),
              tags$li("View gene expression patterns")
            ),
            
            div(class = "alert alert-info", style = "background: #e3f2fd; border-left: 4px solid #2196f3; padding: 15px; margin-top: 20px;",
              strong("Tip:"), " Start with a smaller dataset if you're new to the tool to familiarize yourself with the interface."
            )
          )
        ),
        
        # Section 3: Analysis Workflow
        div(id = "workflow", class = "doc-section", style = "margin-bottom: 40px;",
          h2("3. Analysis Workflow", style = "color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px;"),
          div(style = "background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
            
            h4("Step 1: Import Dataset"),
            p("Select your organism and dataset of interest."),
            
            h4("Step 2: Visualize Gene Expression", style = "margin-top: 20px;"),
            tags$ul(
              tags$li(strong("UMAP plots:"), " Visualize cell populations"),
              tags$li(strong("Violin plots:"), " Compare gene expression across cell types"),
              tags$li(strong("Gene filters:"), " Search and select genes of interest")
            ),
            
            h4("Step 3: Gene Correlation Analysis", style = "margin-top: 20px;"),
            tags$ul(
              tags$li("Select two genes to analyze their co-expression"),
              tags$li("Choose statistical test (Spearman or Pearson)"),
              tags$li("View scatter plots with correlation coefficients"),
              tags$li("Download results as CSV")
            ),
            
            h4("Step 4: Differential Gene Expression (DGE)", style = "margin-top: 20px;"),
            tags$ul(
              tags$li("Define two groups to compare"),
              tags$li("Choose statistical method (Wilcoxon, t-test, etc.)"),
              tags$li("View ranked gene lists with statistics"),
              tags$li("Filter by log fold-change, p-value, or score"),
              tags$li("Export results to CSV")
            ),
            
            h4("Step 5: Enrichment Analysis", style = "margin-top: 20px;"),
            tags$ul(
              tags$li("Select genes from DGE results"),
              tags$li("Choose enrichment database (GO, KEGG, Reactome, etc.)"),
              tags$li("Visualize enriched pathways and terms"),
              tags$li("Download enrichment tables")
            ),
            
            h4("Step 6: Pseudo-bulk Analysis", style = "margin-top: 20px;"),
            tags$ul(
              tags$li("Aggregate cells by sample"),
              tags$li("Run DESeq2 for robust differential expression"),
              tags$li("View volcano plots and results tables"),
              tags$li("Perform enrichment on pseudo-bulk results")
            )
          )
        ),
        
        # Section 4: Key Features
        div(id = "features", class = "doc-section", style = "margin-bottom: 40px;",
          h2("4. Key Features", style = "color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px;"),
          
          div(style = "background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
            
            # Feature cards
            div(class = "row",
              div(class = "col-md-6", style = "margin-bottom: 20px;",
                div(style = "background: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;",
                  h5("Interactive Visualizations"),
                  tags$ul(
                    tags$li("UMAP/t-SNE projections"),
                    tags$li("Violin and box plots"),
                    tags$li("Heatmaps and dot plots"),
                    tags$li("Volcano plots for DGE"),
                    tags$li("Customizable color schemes")
                  )
                )
              ),
              div(class = "col-md-6", style = "margin-bottom: 20px;",
                div(style = "background: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;",
                  h5("Statistical Analysis"),
                  tags$ul(
                    tags$li("Wilcoxon rank-sum test"),
                    tags$li("Student's t-test"),
                    tags$li("Spearman/Pearson correlation"),
                    tags$li("DESeq2 for pseudo-bulk"),
                    tags$li("Multiple testing correction")
                  )
                )
              ),
              div(class = "col-md-6", style = "margin-bottom: 20px;",
                div(style = "background: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #2196f3;",
                  h5("Enrichment Databases"),
                  tags$ul(
                    tags$li("Gene Ontology (GO)"),
                    tags$li("Biological Processes (BP)"),
                    tags$li("KEGG Pathways"),
                    tags$li("Reactome"),
                    tags$li("WikiPathways")
                  )
                )
              ),
              div(class = "col-md-6", style = "margin-bottom: 20px;",
                div(style = "background: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;",
                  h5("Advanced Options"),
                  tags$ul(
                    tags$li("Cell type filtering"),
                    tags$li("Gene set enrichment analysis"),
                    tags$li("Custom gene lists"),
                    tags$li("Batch effect visualization"),
                    tags$li("Quality control metrics")
                  )
                )
              )
            )
          )
        ),
        
        # Section 5: Exporting Results
        div(id = "export", class = "doc-section", style = "margin-bottom: 40px;",
          h2("5. Exporting Results", style = "color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px;"),
          div(style = "background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
            p("All analysis results can be exported as CSV files for further analysis in Excel, R, Python, or other tools."),
            
            h4("Available Exports:"),
            tags$ul(
              tags$li(strong("Cell Type Markers:"), " Download marker genes for selected clusters"),
              tags$li(strong("Correlation Results:"), " Export gene-gene correlation data"),
              tags$li(strong("DGE Results:"), " Save differential expression statistics"),
              tags$li(strong("Enrichment Analysis:"), " Download pathway enrichment tables"),
              tags$li(strong("Pseudo-bulk Results:"), " Export DESeq2 results"),
              tags$li(strong("Pseudo-bulk Enrichment:"), " Save enrichment for pseudo-bulk data")
            ),
            
            h4("How to Export:", style = "margin-top: 20px;"),
            tags$ol(
              tags$li("Complete your analysis in the respective tab"),
              tags$li("Look for the ", tags$code("Download ... (CSV)"), " button below each results table"),
              tags$li("Click the button to download the file"),
              tags$li("Files are named with the analysis type and current date")
            ),
            
            div(class = "alert alert-success", style = "background: #e8f5e9; border-left: 4px solid #4caf50; padding: 15px; margin-top: 20px;",
              strong("Pro Tip:"), " Export your results regularly to keep track of different analyses and comparisons."
            )
          )
        ),
        
        # Section 6: Troubleshooting
        div(id = "troubleshooting", class = "doc-section", style = "margin-bottom: 40px;",
          h2("6. Troubleshooting", style = "color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px;"),
          div(style = "background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
            
            h4("Common Issues and Solutions:"),
            
            div(style = "background: #fff3cd; border-left: 4px solid #ffc107; padding: 15px; margin: 15px 0;",
              h5("Dataset not loading"),
              tags$ul(
                tags$li("Check your internet connection"),
                tags$li("Try refreshing the page"),
                tags$li("Select a different dataset to verify the issue"),
                tags$li("Clear browser cache if problem persists")
              )
            ),
            
            div(style = "background: #fff3cd; border-left: 4px solid #ffc107; padding: 15px; margin: 15px 0;",
              h5("Export button returns error 500"),
              tags$ul(
                tags$li("Ensure you have completed the analysis before exporting"),
                tags$li("Check that the results table contains data"),
                tags$li("Try re-running the analysis"),
                tags$li("If error persists, the data may be too large - try filtering first")
              )
            ),
            
            div(style = "background: #fff3cd; border-left: 4px solid #ffc107; padding: 15px; margin: 15px 0;",
              h5("Enrichment analysis not working"),
              tags$ul(
                tags$li("Verify that genes are selected in the DGE table"),
                tags$li("Check that the organism matches your dataset"),
                tags$li("Ensure the fenr package is installed"),
                tags$li("Try with a different enrichment database")
              )
            ),
            
            div(style = "background: #fff3cd; border-left: 4px solid #ffc107; padding: 15px; margin: 15px 0;",
              h5("Plots not displaying"),
              tags$ul(
                tags$li("Wait for analysis to complete (check for spinner)"),
                tags$li("Ensure required inputs are selected"),
                tags$li("Try zooming in/out in your browser"),
                tags$li("Check browser console for JavaScript errors")
              )
            ),
            
            h4("Need More Help?", style = "margin-top: 30px;"),
            p("If you encounter issues not covered here:"),
            tags$ul(
              tags$li("Check the GitHub repository for known issues"),
              tags$li("Open an issue with a detailed description of the problem"),
              tags$li("Include screenshots and error messages when possible")
            )
          )
        ),
        
        # Footer
        div(class = "doc-footer", style = "text-align: center; margin-top: 50px; padding: 30px; background: #f8f9fa; border-radius: 8px;",
          h4("Citation", style = "color: #2c3e50;"),
          p("If you use MASLDatlas in your research, please cite our publication:"),
          div(style = "background: white; padding: 15px; border-radius: 4px; margin: 15px 0; font-family: monospace; text-align: left;",
            "MASLDatlas: Multi-species single-cell RNA-seq atlas of Metabolic dysfunction-Associated Steatotic Liver Disease. [Journal Name], [Year]."
          ),
          hr(),
          p(style = "margin-top: 20px; color: #7f8c8d;",
            "MASLDatlas v1.0 | ", 
            tags$a(href = "https://github.com/BioMAs/MASLDatlas", target = "_blank", "GitHub Repository"),
            " | Last updated: October 2025"
          )
        )
      )
    )
    # ,
    # tabPanel(
    #   title = "Datasets",
    #   value = "tab_datasets",
    #   icon = icon("database")
    # ),
    # tabPanel(
    #   title = "About",
    #   value = "tab_about",
    #   icon = icon("info")
    # )
  )
)




server <- function(input, output,session) {
  
  # Navigation depuis la page d'accueil
  observeEvent(input$go_to_explore, {
    updateTabsetPanel(session, "navbar", selected = "tab_explore_datasets")
  })
  
  observeEvent(input$go_to_correlations, {
    updateTabsetPanel(session, "navbar", selected = "tab_correlation")
  })
  
  ##################################### Import Dataset
  
  
  # Import Dataset
  ##ui actions
  
  output$dataset_selection_list <- renderUI({
    # Use helper function to get organism data safely
    organism_data <- get_organism_data(datasets_config, input$selection_organism)
    
    # Extract data using the helper function
    datasets <- organism_data$datasets
    status <- organism_data$status
    description <- organism_data$description
    message <- organism_data$message
    
    # Check if datasets are available for this organism
    if(length(datasets) == 0) {
      # No datasets available
      content <- div(
        p(style = "color: #f39c12; font-weight: bold;", "⚠️ No datasets currently available"),
        p(style = "color: #666; font-size: 12px;", description),
        if (!is.null(message)) {
          p(style = "color: #7f8c8d; font-size: 11px; font-style: italic;", message)
        },
        disabled(selectInput("selection_dataset", "Select Dataset", choices = c("No datasets available" = "")))
      )
    } else if(status == "Available") {
      # Datasets available - create named choices with description
      choices <- setNames(datasets, paste0(datasets, " - Available"))
      content <- div(
        p(style = "color: #27ae60; font-weight: bold;", paste0("✅ Status: ", status)),
        p(style = "color: #666; font-size: 12px;", description),
        selectInput("selection_dataset", "Select Dataset", choices = choices)
      )
    } else {
      # Datasets not ready (Optional, In Progress, etc.)
      content <- div(
        p(style = "color: #f39c12; font-weight: bold;", paste0("⏳ Status: ", status)),
        p(style = "color: #666; font-size: 12px;", description),
        if (!is.null(message)) {
          p(style = "color: #7f8c8d; font-size: 11px; font-style: italic;", message)
        },
        disabled(selectInput("selection_dataset", "Select Dataset", choices = setNames("", paste0("Status: ", status))))
      )
    }
    
    return(content)
  })
  
  adata <- eventReactive(input$import_dataset, {
    # Validate inputs
    if (is.null(input$selection_organism) || is.null(input$selection_dataset) || input$selection_dataset == "") {
      showNotification("❌ Please select a valid organism and dataset", type = "error")
      return(NULL)
    }

    # Check if dataset name indicates unavailability
    if (grepl("No datasets available|Error|Contact admin|Status:", input$selection_dataset, ignore.case = TRUE)) {
      showNotification("❌ This dataset is not currently available", type = "error", duration = 10)
      return(NULL)
    }

    # ⚡ OPTIMIZATION: Check cache first for faster loading
    cache_key <- paste(input$selection_organism, input$selection_dataset, 
                      input$dataset_size_option %||% "full", sep = "_")
    
    if (exists("get_cached_dataset", mode = "function")) {
      cached_data <- get_cached_dataset(cache_key)
      if (!is.null(cached_data)) {
        showNotification("✅ Dataset loaded from cache (fast loading enabled)", type = "message")
        
        return(cached_data)
      }
    }    # Use validation function to check dataset path
    dataset_info <- validate_dataset_path(input$selection_organism, input$selection_dataset)
    
    # Check if dataset file exists before proceeding
    if (!dataset_info$exists) {
      # Get organism data to provide more informative error message
      organism_data <- get_organism_data(datasets_config, input$selection_organism)
      
      error_message <- paste0("❌ Dataset file not found: ", basename(dataset_info$path))
      if (organism_data$status != "Available") {
        error_message <- paste0(error_message, "\n💡 Status: ", organism_data$status)
        if (!is.null(organism_data$download_command)) {
          error_message <- paste0(error_message, "\n📥 To download: ", organism_data$download_command)
        }
      }
      
      showNotification(error_message, type = "error", duration = 15)
      return(NULL)
    }
    
    dataset_path <- dataset_info$path
    
    # Check if this is the large integrated dataset
    is_large_dataset <- grepl("Fibrotic.*Cross.*Species.*002", input$selection_dataset)
    
    if (is_large_dataset && !is.null(input$dataset_size_option) && input$dataset_size_option != "full") {
      # Use optimized version based on user selection
      size_suffix <- switch(input$dataset_size_option,
                           "sub5k" = "_sub5k",
                           "sub10k" = "_sub10k", 
                           "sub20k" = "_sub20k",
                           "")
      
      optimized_path <- paste0("datasets_optimized/", 
                              tools::file_path_sans_ext(input$selection_dataset), 
                              size_suffix, ".h5ad")
      
      if (file.exists(optimized_path)) {
        dataset_path <- optimized_path
        showNotification(
          paste("✅ Loading", input$dataset_size_option, "version for better performance"), 
          type = "message"
        )
      } else {
        # If optimized version doesn't exist, offer to create it
        showModal(modalDialog(
          title = "Optimized Dataset Not Found",
          p("The selected dataset size is not available. You have two options:"),
          tags$ol(
            tags$li("Use the full dataset (may take 30+ minutes to load)"),
            tags$li("Create optimized versions using the dataset optimization script")
          ),
          footer = tagList(
            actionButton("use_full_dataset", "Use Full Dataset", class = "btn-warning"),
            modalButton("Cancel")
          )
        ))
        return(NULL)
      }
    } else if (is_large_dataset && (is.null(input$dataset_size_option) || input$dataset_size_option == "full")) {
      # Show strong warning for full dataset
      showNotification(
        "⚠️ Loading full 9.2GB dataset. This may take 30+ minutes and use significant memory.",
        type = "warning",
        duration = 15
      )
    }
    
    # Show loading progress
    progress <- Progress$new()
    progress$set(message = "Loading dataset...", value = 0)
    on.exit(progress$close())
    
    progress$set(value = 0.3, detail = "Reading file...")
    
    tryCatch({
      # Standard dataset loading
      adata <- sc$read_h5ad(dataset_path)
      
      progress$set(value = 0.8, detail = "Processing metadata...")
      
      progress$set(value = 1, detail = "Complete!")
      
      showNotification(
        paste("✅ Dataset loaded successfully:", 
              format(adata$n_obs, big.mark = ","), "cells ×", 
              format(adata$n_vars, big.mark = ","), "genes"),
        type = "message"
      )
      
      # Return the loaded dataset
      return(adata)
      
    }, error = function(e) {
      progress$close()
      
      showNotification(
        paste("❌ Error loading dataset:", e$message),
        type = "error",
        duration = NULL
      )
      return(NULL)
    })
  })
  
  
  #ui outputs
    output$textoutput_data_numbers <- renderText({
      req(adata())
      paste0("Number of Cells: ", adata()$n_obs,"   Number of Genes: ", adata()$n_vars)
    })
  
    output$imageoutput_UMAP <- renderImage({
    req(adata())
    sc$pl$umap(adata(),color = c('CellType'), legend_loc = "on data", show=FALSE,save = '.png')
    list(src = "figures/umap.png")
    }, deleteFile = TRUE)
    
    
    
    
    output$umap_col_selection <- renderUI({
      req(adata())
      # Use the column names as choices in selectInput
      selectInput("selection_umap_identity_select", "Select Identity", choices = colnames(adata()$obs))
    })
    
    output$imageoutput_UMAP_selection <- renderImage({
      req(adata(),input$selection_umap_identity_select)
      
      if(input$selection_umap_identity_select != "over_clustering"){
        sc$pl$umap(adata(),color = input$selection_umap_identity_select, show=FALSE, save = 'second.png')
      }else{
        n_clusters <- length(unique(adata()$obs[[input$selection_umap_identity_select]]))
        pal1 <- reticulate::py_to_r(sc$pl$palettes$default_102)
        pal2 <- reticulate::py_to_r(sc$pl$palettes$vega_20_scanpy)
        colors <- c(pal1, pal2)
        if (length(colors) < n_clusters) {
          colors <- rep(colors, length.out = n_clusters)
        }
        sc$pl$umap(adata(),color = input$selection_umap_identity_select, palette = colors[1:n_clusters], show=FALSE, legend_loc = "none", save = 'second.png')
      }
      list(src = "figures/umapsecond.png")
    }, deleteFile = TRUE)
  
    
    
    output$ranks_col_selection <- renderUI({
      req(adata())
      # Use the column names as choices in selectInput
      cluster_names <- unique(adata()$obs$CellType)
      cluster_names <- sort(cluster_names)
      cluster_names <- as.character(cluster_names)
      selectInput("selection_rank_select", "Select Cluster", choices = cluster_names)
    })
    
    
    

    output$imageoutput_CellType_groups <- renderImage({
      req(adata(),input$selection_rank_select)
      group_name <- list(input$selection_rank_select)
      # Convert the group number to a string
      sc$pl$rank_genes_groups(adata(), groups = group_name, sharey = FALSE, show=FALSE, save = '.png')
      list(src = "figures/rank_genes_groups_CellType.png")
    }, deleteFile = TRUE)
    
    

    
    output$table_cluster_markers <- renderDT({
      req(adata(),input$selection_rank_select)
      data_table_dges <- sc$get$rank_genes_groups_df(adata(), group= input$selection_rank_select)
      datatable(
        data_table_dges, 
        options = list(
          pageLength = 5,
          autoWidth = TRUE,
          scrollX = TRUE,
          scrollY = 300,  # Set the height for the scroll
          scroller = TRUE,
          dom = 'l<"toolbar">rtip',  # Removed 'f' from dom
          deferRender = TRUE
        ),
        filter = 'top',
        rownames = FALSE,
        extensions = 'Scroller'
      )%>%
        formatRound(which(sapply(data_table_dges, is.numeric)), digits = 3)
    })
    
    
    ##################################### Cluster Selection
    
    output$cluster_selection_list <- renderUI({
      req(adata())
      
      choices_for_selector <- unique(adata()$obs$CellType)
      choices_for_selector <- sort(choices_for_selector)
      choices_for_selector <- as.character(choices_for_selector) 
      selectizeInput(
        "selected_culusters", 
          "Choose Cluster",
        choices = choices_for_selector, 
        multiple = TRUE, options = list(placeholder = "Select Clusters to Filter"))
    })
    
    
    output$button_cluster_filter <- renderUI({
      req(adata(),input$selected_culusters)
      actionButton("filter_dataset_cluster_selection", "Filter Dataset", class = "btn-primary", icon = icon("filter"), width = '100%')
      
    })

    
    filtered_adata <- eventReactive(input$filter_dataset_cluster_selection, {

      selected_culusters <- list(input$selected_culusters)
      # filtered_adata <- adata()[adata()$obs['CellType'] == selected_culusters]
      
      filtered_adata = adata()[adata()$obs$CellType %in% input$selected_culusters, ]
      
      return(filtered_adata)
    })
    
    
    gene_list_adata <- eventReactive(c(input$import_dataset,input$filter_dataset_cluster_selection), {
      req(adata())
      
      if(is.null(input$filter_dataset_cluster_selection)){
        expression_table <- adata()$var
        gene_list_adata <- rownames(expression_table)
      }else{
        expression_table <- filtered_adata()$var
        gene_list_adata <- rownames(expression_table)
      }
      return(gene_list_adata)
    })
    
    output$gene_list_expression <- renderUI({
      if(input$cluster_selection_visualization_type == 'Visualize Expression of Gene'){
        
      req(gene_list_adata())
        gene_list_adata <- sort(gene_list_adata())
      selectizeInput(
        "gene_selection_cluster_expression", "Select Gene", choices = gene_list_adata
      )}
    })
    
    output$button_visualize_cluster_selection <- renderUI({
      req(gene_list_adata())
      actionButton("visualize_cluster_selection", "Visualize",icon = icon("eye"), width = '100%', class = "btn-primary")
    })
    
    
    output$umap_col_selection_clusters <- renderUI({
      if(input$cluster_selection_visualization_type == 'Visualize Expression of Gene'){
        
      req(input$visualize_cluster_selection)
      # Use the column names as choices in selectInput
      selectInput("selection_umap_identity_select_clusters", "Select Identity", choices = colnames(adata()$obs),
                  selected = 'CellType')}
    })
    
    
    output$info_cluster_selection <- renderText({
      if(input$cluster_selection_visualization_type == 'Visualize Expression of Gene'){
        
      req(input$visualize_cluster_selection)
      
      if(is.null(input$filter_dataset_cluster_selection)){
        paste0("Number of Cells: ", adata()$n_obs,"   Number of Genes: ", adata()$n_vars)
      }else{
        paste0("Number of Cells: ", filtered_adata()$n_obs,"   Number of Genes: ", filtered_adata()$n_vars)
      }
        }
      
      
    })
    
    
    output$imageoutput_UMAP_selection_clusters <- renderImage({
      if(input$cluster_selection_visualization_type == 'Visualize Expression of Gene'){
        
      req(input$selection_umap_identity_select_clusters,input$visualize_cluster_selection,input$gene_selection_cluster_expression)
      if(is.null(input$filter_dataset_cluster_selection)){
        
        if(input$selection_umap_identity_select_clusters != "over_clustering"){
          sc$pl$umap(adata(),color = input$selection_umap_identity_select_clusters, show=FALSE, save = 'clusters.png')
        }else{
          n_clusters <- length(unique(adata()$obs[[input$selection_umap_identity_select_clusters]]))
          pal1 <- reticulate::py_to_r(sc$pl$palettes$default_102)
          pal2 <- reticulate::py_to_r(sc$pl$palettes$vega_20_scanpy)
          colors <- c(pal1, pal2)
          if (length(colors) < n_clusters) {
            colors <- rep(colors, length.out = n_clusters)
          }
          sc$pl$umap(adata(),color = input$selection_umap_identity_select_clusters, palette = colors[1:n_clusters], show=FALSE, legend_loc = "none", save = 'clusters.png')
        }
        list(src = "figures/umapclusters.png")
        
      }else{
        
        if(input$selection_umap_identity_select_clusters != "over_clustering"){
          sc$pl$umap(filtered_adata(),color = input$selection_umap_identity_select_clusters, show=FALSE, save = 'clusters.png')
        }else{
          n_clusters <- length(unique(adata()$obs[[input$selection_umap_identity_select_clusters]]))
          pal1 <- reticulate::py_to_r(sc$pl$palettes$default_102)
          pal2 <- reticulate::py_to_r(sc$pl$palettes$vega_20_scanpy)
          colors <- c(pal1, pal2)
          if (length(colors) < n_clusters) {
            colors <- rep(colors, length.out = n_clusters)
          }
          sc$pl$umap(filtered_adata(),color = input$selection_umap_identity_select_clusters, palette = colors[1:n_clusters], show=FALSE, legend_loc = "none",  save = 'clusters.png')
        }
        
        list(src = "figures/umapclusters.png")
      }
      }
      
    }, deleteFile = TRUE)
    
    output$imageoutput_UMAP_expression_clusters <- renderImage({
      if(input$cluster_selection_visualization_type == 'Visualize Expression of Gene'){
        
      req(input$visualize_cluster_selection,input$gene_selection_cluster_expression)
      if(is.null(input$filter_dataset_cluster_selection)){
        sc$pl$umap(adata(),color = input$gene_selection_cluster_expression, layer = 'scvi_normalized', vmax = 5, show=FALSE, save = 'clusters_exp.png')
        list(src = "figures/umapclusters_exp.png")
      }else{
        sc$pl$umap(filtered_adata(),color = input$gene_selection_cluster_expression, layer = 'scvi_normalized', vmax = 5, show=FALSE, save = 'clusters_exp.png')
        list(src = "figures/umapclusters_exp.png")
      }
      }
    }, deleteFile = TRUE)
    
    output$imageoutput_violin_expression_clusters <- renderImage({
      if(input$cluster_selection_visualization_type == 'Visualize Expression of Gene'){
        
      req(input$visualize_cluster_selection,input$gene_selection_cluster_expression)
      if(is.null(input$filter_dataset_cluster_selection)){
        sc$pl$violin(adata(), keys = input$gene_selection_cluster_expression, groupby = 'CellType', use_raw=F, layer = 'scvi_normalized', show=FALSE, rotation=90, save = "violin_exp.png")
        list(src = "figures/violinviolin_exp.png")
      }else{
        sc$pl$violin(filtered_adata(), keys = input$gene_selection_cluster_expression, groupby = 'CellType', use_raw=F, layer = 'scvi_normalized', show=FALSE, rotation=90, save = 'violin_exp.png')
        list(src = "figures/violinviolin_exp.png")
      }
      }
    }, deleteFile = TRUE)
    
    
    output$imageoutput_violin_expression_clusters_groups <- renderImage({
      if(input$cluster_selection_visualization_type == 'Visualize Expression of Gene'){
        
      req(input$visualize_cluster_selection,input$gene_selection_cluster_expression)
      if(is.null(input$filter_dataset_cluster_selection)){
        sc$pl$violin(adata(), keys = input$gene_selection_cluster_expression, groupby = 'Group', use_raw=F, layer = 'scvi_normalized', show=FALSE, rotation=90, save = "clusters_violin_exp.png")
        list(src = "figures/violinclusters_violin_exp.png")
      }else{
        sc$pl$violin(filtered_adata(), keys = input$gene_selection_cluster_expression, groupby = 'Group', use_raw=F, layer = 'scvi_normalized', show=FALSE, rotation=90, save = 'clusters_violin_exp.png')
        list(src = "figures/violinclusters_violin_exp.png")
      }
      }
    }, deleteFile = TRUE)
    
    ################################################################ Visualize genesets
    

    # observeEvent(input$first_geneset_list, {
    #   genes <- input$first_geneset_list
    # 
    #   # Remove special characters except for "-"
    #   genes <- str_replace_all(genes, "[^A-Za-z0-9\\-\\s,]", "")
    # 
    #   # Split by various separators
    #   gene_list <- str_split(genes, "\\s+|,|\\t")[[1]]
    # 
    #   # Remove empty strings and duplicates, and trim whitespace
    #   gene_list <- unique(trimws(gene_list[gene_list != ""]))
    # 
    #   if(is.null(input$filter_dataset_cluster_selection)){
    #     expression_table <- adata()$var
    #     gene_list_adata <- rownames(expression_table)
    #   }else{
    #     expression_table <- filtered_adata()$var
    #     gene_list_adata <- rownames(expression_table)
    #   }
    # 
    #   gene_list <- intersect(gene_list, gene_list_adata)
    # 
    #   # Create a cleaned string with one gene name per line
    #   cleaned_genes <- paste(gene_list, collapse = "\n")
    # 
    # 
    #   # Update the text area input with the cleaned gene names
    # 
    #   updateTextAreaInput(session, "first_geneset_list", value = cleaned_genes)
    # })
    
    
    first_gene_set_calc <- eventReactive(input$visualize_cluster_selection,{
      req(input$cluster_selection_visualization_type== "Visualize Expression of Geneset", input$first_geneset_list)
      
      genes <- input$first_geneset_list
      
      # Remove special characters except for "-"
      genes <- str_replace_all(genes, "[^A-Za-z0-9\\-\\s,]", "")
      
      # Split by various separators
      gene_list <- str_split(genes, "\\s+|,|\\t")[[1]]
      
      # Remove empty strings and duplicates, and trim whitespace
      gene_list <- unique(trimws(gene_list[gene_list != ""]))
      
      if(is.null(input$filter_dataset_cluster_selection)){
        expression_table <- adata()$var
        gene_list_adata <- rownames(expression_table) 
      }else{
        expression_table <- filtered_adata()$var
        gene_list_adata <- rownames(expression_table)
      }
      
      gene_list <- intersect(gene_list, gene_list_adata)
      
      # Create a cleaned string with one gene name per line
      cleaned_genes <- paste(gene_list, collapse = "\n")
      
      
      # Update the text area input with the cleaned gene names
      
      updateTextAreaInput(session, "first_geneset_list", value = cleaned_genes)
      
      gene_list <- strsplit(input$first_geneset_list, split = "\n")[[1]]
      
      if(!is.null(input$name_first_geneset)){
        net <- data.frame(
          source = rep(input$name_first_geneset, length(gene_list)),
          target = c(gene_list)
        )
      }else{
        net <- data.frame(
          source = rep("First Geneset", length(gene_list)),
          target = c(gene_list)
        )
      }
      
    
      if(is.null(input$filter_dataset_cluster_selection)){
        if(is.null(adata()$raw)){
          dc$run_ulm(adata(), net, weight=NULL, min_n=0, use_raw = FALSE)
          
        }else{
          dc$run_ulm(adata(), net, weight=NULL, min_n=0)
        }
        acts = dc$get_acts(adata(), 'ulm_estimate')
        acts_matrix <- acts$X
        acts_v <- as.vector(acts$X)
        max_e <- max(acts_v[is.finite(acts_v)], na.rm = TRUE)
        acts_matrix[!is.finite(acts_matrix)] <- max_e
        acts$X <- acts_matrix
        acts
        
        return(acts)
        
      }else{
        if(is.null(filtered_adata()$raw)){
          dc$run_ulm(filtered_adata(), net, weight=NULL, min_n=0, use_raw = FALSE)
          
        }else{
          dc$run_ulm(filtered_adata(), net, weight=NULL, min_n=0)
        }
        acts = dc$get_acts(filtered_adata(), 'ulm_estimate')
        acts_matrix <- acts$X
        acts_v <- as.vector(acts$X)
        max_e <- max(acts_v[is.finite(acts_v)], na.rm = TRUE)
        acts_matrix[!is.finite(acts_matrix)] <- max_e
        acts$X <- acts_matrix
        acts
        
        
        return(acts)
        
      }
      
    })
    
    
    output$image_output_enrichment_first_set <- renderImage({
      req(input$visualize_cluster_selection,input$first_geneset_list, input$cluster_selection_visualization_type == 'Visualize Expression of Geneset')
      if(!is.null(input$name_first_geneset)){
        sc$pl$umap(first_gene_set_calc(), color= rownames(first_gene_set_calc()$var), cmap='RdBu_r', show = F, vmax = 5, save = "first_gene_enrichment_umap.png")
      }else{
        sc$pl$umap(first_gene_set_calc(), color= rownames(first_gene_set_calc()$var), cmap='RdBu_r', show = F, vmax = 5, save = "first_gene_enrichment_umap.png")
      }
      list(src = "figures/umapfirst_gene_enrichment_umap.png")
      
    }, deleteFile = TRUE)
    
    
    output$image_output_enrichment_first_violin_clusters <- renderImage({
      req(input$visualize_cluster_selection,input$first_geneset_list, input$cluster_selection_visualization_type == 'Visualize Expression of Geneset')
      if(!is.null(input$name_first_geneset)){
        sc$pl$violin(first_gene_set_calc(), keys = rownames(first_gene_set_calc()$var), groupby='CellType', show = F, rotation=90, save = "first_gene_enrichment_violin.png")
      }else{
        sc$pl$violin(first_gene_set_calc(), keys = rownames(first_gene_set_calc()$var), groupby='CellType', show = F, rotation=90, save = "first_gene_enrichment_violin.png")
      }
      list(src = "figures/violinfirst_gene_enrichment_violin.png")
      
    }, deleteFile = TRUE)
    
    output$image_output_enrichment_first_violin_groups <- renderImage({
      req(input$visualize_cluster_selection,input$first_geneset_list, input$cluster_selection_visualization_type == 'Visualize Expression of Geneset')
      if(!is.null(input$name_first_geneset)){
        sc$pl$violin(first_gene_set_calc(), keys = rownames(first_gene_set_calc()$var), groupby='Group', show = F, rotation=90, save = "first_gene_enrichment_violin_group.png")
      }else{
        sc$pl$violin(first_gene_set_calc(), keys = rownames(first_gene_set_calc()$var), groupby='Group', show = F, rotation=90, save = "first_gene_enrichment_violin_group.png")
      }
      list(src = "figures/violinfirst_gene_enrichment_violin_group.png")
      
    }, deleteFile = TRUE)
    
    
    output$geneset_name_input_second <- renderUI({
      req(input$visualize_cluster_selection)
      textInput("name_second_geneset", "Name second Geneset", placeholder = "Enter Name", width = '100%')
    })
    
    output$geneset_list_input_second <- renderUI({
      req(input$visualize_cluster_selection)
      textAreaInput("second_geneset_list", "Enter Gene Names", placeholder = "Enter Names", width = '100%')
    })
    
    output$geneset_list_go_second <- renderUI({
      req(input$second_geneset_list)
      actionButton("second_geneset_run", "Visualize", width = '100%', icon = icon("eye"), class = "btn-primary")
    })
    
    # observeEvent(input$second_geneset_list, {
    #   genes <- input$second_geneset_list
    #   
    #   # Remove special characters except for "-"
    #   genes <- str_replace_all(genes, "[^A-Za-z0-9\\-\\s,]", "")
    #   
    #   # Split by various separators
    #   gene_list <- str_split(genes, "\\s+|,|\\t")[[1]]
    #   
    #   # Remove empty strings and duplicates, and trim whitespace
    #   gene_list <- unique(trimws(gene_list[gene_list != ""]))
    #   
    #   gene_list <- intersect(gene_list, gene_list_adata())
    #   
    #   # Create a cleaned string with one gene name per line
    #   cleaned_genes <- paste(gene_list, collapse = "\n")
    #   
    #   
    #   # Update the text area input with the cleaned gene names
    #   
    #   updateTextAreaInput(session, "second_geneset_list", value = cleaned_genes)
    # })
    
    
    
    
    
    second_gene_set_calc <- eventReactive(input$visualize_cluster_selection,{
      
      genes <- input$second_geneset_list
      
      # Remove special characters except for "-"
      genes <- str_replace_all(genes, "[^A-Za-z0-9\\-\\s,]", "")
      
      # Split by various separators
      gene_list <- str_split(genes, "\\s+|,|\\t")[[1]]
      
      # Remove empty strings and duplicates, and trim whitespace
      gene_list <- unique(trimws(gene_list[gene_list != ""]))
      
      gene_list <- intersect(gene_list, gene_list_adata())
      
      # Create a cleaned string with one gene name per line
      cleaned_genes <- paste(gene_list, collapse = "\n")
      
      
      # Update the text area input with the cleaned gene names
      
      updateTextAreaInput(session, "second_geneset_list", value = cleaned_genes)
      
      req(input$cluster_selection_visualization_type== "Visualize Expression of Geneset", input$second_geneset_run)
      
      gene_list <- strsplit(input$second_geneset_list, split = "\n")[[1]]
      
      if(!is.null(input$name_second_geneset)){
        net <- data.frame(
          source = rep(input$name_second_geneset, length(gene_list)),
          target = gene_list
        )
      }else{
        net <- data.frame(
          source = rep("Second Geneset", length(gene_list)),
          target = gene_list
        )
      }
      
      
      
      if(is.null(input$filter_dataset_cluster_selection)){
        if(is.null(adata()$raw)){
          dc$run_ulm(adata(), net, weight=NULL, min_n=0, use_raw = FALSE)
          
        }else{
          dc$run_ulm(adata(), net, weight=NULL, min_n=0)
        }
        acts = dc$get_acts(adata(), 'ulm_estimate')
        acts_matrix <- acts$X
        acts_v <- as.vector(acts$X)
        max_e <- max(acts_v[is.finite(acts_v)], na.rm = TRUE)
        acts_matrix[!is.finite(acts_matrix)] <- max_e
        acts$X <- acts_matrix
        acts
        
        return(acts)
        
      }else{
        if(is.null(filtered_adata()$raw)){
          dc$run_ulm(filtered_adata(), net, weight=NULL, min_n=0, use_raw = FALSE)
          
        }else{
          dc$run_ulm(filtered_adata(), net, weight=NULL, min_n=0)
        }
        acts = dc$get_acts(filtered_adata(), 'ulm_estimate')
        acts_matrix <- acts$X
        acts_v <- as.vector(acts$X)
        max_e <- max(acts_v[is.finite(acts_v)], na.rm = TRUE)
        acts_matrix[!is.finite(acts_matrix)] <- max_e
        acts$X <- acts_matrix
        acts
        
        return(acts)
        
      }
      
    })
    
    
    output$image_output_enrichment_second_set <- renderImage({
      req(input$second_geneset_run,input$second_geneset_list, input$cluster_selection_visualization_type == 'Visualize Expression of Geneset')
      if(!is.null(input$name_second_geneset)){
        sc$pl$umap(second_gene_set_calc(), color= rownames(second_gene_set_calc()$var), cmap='RdBu_r', show = F,vmax = 5, save = "second_gene_enrichment_umap.png")
      }else{
        sc$pl$umap(second_gene_set_calc(), color= rownames(second_gene_set_calc()$var), cmap='RdBu_r', show = F,vmax = 5, save = "second_gene_enrichment_umap.png")
      }
      list(src = "figures/umapsecond_gene_enrichment_umap.png")
      
    }, deleteFile = TRUE)
    
    
    output$image_output_enrichment_second_violin_clusters <- renderImage({
      req(input$second_geneset_run,input$second_geneset_list, input$cluster_selection_visualization_type == 'Visualize Expression of Geneset')
      if(!is.null(input$name_second_geneset)){
        sc$pl$violin(second_gene_set_calc(), keys = rownames(second_gene_set_calc()$var), groupby='CellType', show = F, rotation=90, save = "second_gene_enrichment_violin.png")
      }else{
        sc$pl$violin(second_gene_set_calc(), keys = rownames(second_gene_set_calc()$var), groupby='CellType', show = F, rotation=90, save = "second_gene_enrichment_violin.png")
      }
      list(src = "figures/violinsecond_gene_enrichment_violin.png")
      
    }, deleteFile = TRUE)
    
    output$image_output_enrichment_second_violin_groups <- renderImage({
      req(input$second_geneset_run,input$second_geneset_list, input$cluster_selection_visualization_type == 'Visualize Expression of Geneset')
      if(!is.null(input$name_second_geneset)){
        sc$pl$violin(second_gene_set_calc(), keys = rownames(second_gene_set_calc()$var), groupby='Group', show = F, rotation=90, save = "second_gene_enrichment_violin_group.png")
      }else{
        sc$pl$violin(second_gene_set_calc(), keys = rownames(second_gene_set_calc()$var), groupby='Group', show = F, rotation=90, save = "second_gene_enrichment_violin_group.png")
      }
      list(src = "figures/violinsecond_gene_enrichment_violin_group.png")
      
    }, deleteFile = TRUE)
    
    
    
    
    
    
    
    
    
    
    ## co-expression
    
    output$gene_list_coexpresion_first <- renderUI({
      if(input$cluster_selection_visualization_type == 'Calculate Co-Expression'){
        
        req(gene_list_adata())
        gene_list_adata <- sort(gene_list_adata())
        
        selectizeInput(
          "gene_selection_cluster_coexpression_first", "Select First Gene", choices = gene_list_adata
        )}
    })
    
    output$gene_list_coexpresion_second <- renderUI({
      if(input$cluster_selection_visualization_type == 'Calculate Co-Expression'){
        
        req(gene_list_adata())
        gene_list_adata <- sort(gene_list_adata())
        selectizeInput(
          "gene_selection_cluster_coexpression_second", "Select Second Gene", choices = gene_list_adata
        )}
    })
    
    output$imageoutput_UMAP_coexpression_first <- renderImage({
      if(input$cluster_selection_visualization_type == 'Calculate Co-Expression'){
        
        req(input$visualize_cluster_selection, input$gene_selection_cluster_coexpression_first)
        if(is.null(input$filter_dataset_cluster_selection)){
          sc$pl$umap(adata(),color = input$gene_selection_cluster_coexpression_first, layer = 'scvi_normalized', vmax = 5, show=FALSE, save = 'coexp_1.png')
          list(src = "figures/umapcoexp_1.png")
        }else{
          sc$pl$umap(filtered_adata(),color = input$gene_selection_cluster_coexpression_first, layer = 'scvi_normalized', vmax = 5, show=FALSE, save = 'coexp_1.png')
          list(src = "figures/umapcoexp_1.png")
        }
      }
    }, deleteFile = TRUE)
    
    output$imageoutput_UMAP_coexpression_second <- renderImage({
      if(input$cluster_selection_visualization_type == 'Calculate Co-Expression'){
        
        req(input$visualize_cluster_selection,input$gene_selection_cluster_coexpression_second)
        if(is.null(input$filter_dataset_cluster_selection)){
          sc$pl$umap(adata(),color = input$gene_selection_cluster_coexpression_second, layer = 'scvi_normalized', vmax = 5, show=FALSE, save = 'coexp_2.png')
          list(src = "figures/umapcoexp_2.png")
        }else{
          sc$pl$umap(filtered_adata(),color = input$gene_selection_cluster_coexpression_second, layer = 'scvi_normalized', vmax = 5, show=FALSE, save = 'coexp_2.png')
          list(src = "figures/umapcoexp_2.png")
        }
      }
    }, deleteFile = TRUE)
    
    statistics_coexpression <- reactive({
      req(input$gene_selection_cluster_coexpression_first,input$gene_selection_cluster_coexpression_second)
      if(is.null(input$filter_dataset_cluster_selection)){
        counts <- adata()$X
        first_gene_count <- counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_first)]
        second_gene_count <- counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_second)]
      }else{
        counts <- filtered_adata()$X
        first_gene_count <- counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_first)]
        second_gene_count <- counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_second)]
      }
      
      
      
      if(input$gene_selection_cluster_coexpression_first == input$gene_selection_cluster_coexpression_second){
        first_gene <- input$gene_selection_cluster_coexpression_first
        second_gene = paste0(input$gene_selection_cluster_coexpression_second,"_")
      }else{
        first_gene = input$gene_selection_cluster_coexpression_first
        second_gene = input$gene_selection_cluster_coexpression_second
      }
      statistics_coexpression <- cbind(first_gene_count, second_gene_count)
      statistics_coexpression <- as.data.frame(statistics_coexpression)
      colnames(statistics_coexpression) <- c(first_gene, second_gene)
      results <- list(statistics_coexpression,first_gene, second_gene)
      return(results)
    })
    
    # output$spearman_results <- renderText({
    #   req(statistics_coexpression(), input$visualize_cluster_selection, input$gene_selection_cluster_coexpression_first,input$gene_selection_cluster_coexpression_second)
    #   paste(paste0(input$gene_selection_cluster_coexpression_first, "-", input$gene_selection_cluster_coexpression_second),
    #     paste0("Test Statistic: ", statistics_coexpression()[[1]]$statistic),
    #         paste0("Spearman's rho: ", statistics_coexpression()[[1]]$estimate),
    #         paste0("p-value: ", statistics_coexpression()[[1]]$p.value),
    #         collapse = "\n")
    #   
    # })
    
    output$stats_cor_plot <- renderPlot({
      req(statistics_coexpression(), input$visualize_cluster_selection)
      
      if (input$remove_zero_counts == F & input$test_choice == "Spearman") {
        data_for_plot <- statistics_coexpression()[[1]]
        x_col <- colnames(data_for_plot)[1]
        y_col <- colnames(data_for_plot)[2]
        
        data_for_plot <- data_for_plot %>%
          dplyr::rename(x_val = !!x_col, y_val = !!y_col)
        
        ggscatter(data = data_for_plot,
                  x = "x_val",
                  y = "y_val",
                  title = "Spearman Correlation",
                  add = "reg.line",
                  conf.int = TRUE,
                  add.params = list(color = "2c3e50", fill = "lightgray")) +
          stat_cor(method = "spearman", label.x = 0.5, label.sep = "\n", size = 10, color = "red")
        
      } else if (input$remove_zero_counts == T && input$test_choice == "Spearman") {
        data_for_plot <- statistics_coexpression()[[1]] %>%
          filter_all(all_vars(. != 0))
        x_col <- statistics_coexpression()[[2]]
        y_col <- statistics_coexpression()[[3]]
        
        data_for_plot <- data_for_plot %>%
          dplyr::rename(x_val = !!x_col, y_val = !!y_col)
        
        ggscatter(data = data_for_plot,
                  x = "x_val",
                  y = "y_val",
                  title = "Spearman Correlation",
                  add = "reg.line",
                  conf.int = TRUE,
                  add.params = list(color = "2c3e50", fill = "lightgray")) +
          stat_cor(method = "spearman", label.x = 0.5, label.sep = "\n", size = 10, color = "red")
        
      } else if (input$remove_zero_counts == F && input$test_choice == "Pearson") {
        data_for_plot <- statistics_coexpression()[[1]]
        x_col <- statistics_coexpression()[[2]]
        y_col <- statistics_coexpression()[[3]]
        
        data_for_plot <- data_for_plot %>%
          dplyr::rename(x_val = !!x_col, y_val = !!y_col)
        
        ggscatter(data = data_for_plot,
                  x = "x_val",
                  y = "y_val",
                  title = "Pearson Correlation",
                  add = "reg.line",
                  conf.int = TRUE,
                  add.params = list(color = "2c3e50", fill = "lightgray")) +
          stat_cor(method = "pearson", label.x = 0.5, label.sep = "\n", size = 10, color = "red")
        
      } else if (input$remove_zero_counts == T && input$test_choice == "Pearson") {
        data_for_plot <- statistics_coexpression()[[1]] %>%
          filter_all(all_vars(. != 0))
        x_col <- statistics_coexpression()[[2]]
        y_col <- statistics_coexpression()[[3]]
        
        data_for_plot <- data_for_plot %>%
          dplyr::rename(x_val = !!x_col, y_val = !!y_col)
        
        ggscatter(data = data_for_plot,
                  x = "x_val",
                  y = "y_val",
                  title = "Pearson Correlation",
                  add = "reg.line",
                  conf.int = TRUE,
                  add.params = list(color = "2c3e50", fill = "lightgray")) +
          stat_cor(method = "pearson", label.x = 0.5, label.sep = "\n", size = 10, color = "red")
      }
    })
    
    ## calculate top correlated genes
    
    
    output$top_correlated_genes <- renderUI({
      req(input$visualize_cluster_selection, statistics_coexpression())
      fluidRow(
        column(width = 12,
               div(style = "background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 0.25rem; padding: 10px; margin-bottom: 15px;",
                   HTML("<small><i class='fas fa-info-circle'></i> Note: Correlation analysis may take 2-5 minutes for large datasets. The analysis is limited to the 1000 most variable genes for performance.</small>"))),
        column(width = 6,
               actionButton("top_correlated_first_gene",paste0("Find Correlated Genes with ", input$gene_selection_cluster_coexpression_first), icon = icon("chart-line"),class = "btn-primary", width = '100%')),
        column(width = 6,
               actionButton("top_correlated_second_gene",paste0("Find Correlated Genes with ", input$gene_selection_cluster_coexpression_second), icon = icon("chart-line"),class = "btn-primary", width = '100%'))
      )
    })
    
    
    correlation_table_first_gene <- eventReactive(input$top_correlated_first_gene,{
      
      if(is.null(input$filter_dataset_cluster_selection)){
        
        # Standard correlation analysis
        normalized_counts <- as.matrix(adata()$X)
      normalized_counts <- as.data.frame(normalized_counts)
      colnames(normalized_counts) <- gene_list_adata()
      first_gene_count <- normalized_counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_first)]
      
      # Optimization: Limit to most variable genes for faster computation
      if(ncol(normalized_counts) > 1000) {
        gene_vars <- apply(normalized_counts, 2, var, na.rm = TRUE)
        top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:1000])
        normalized_counts <- normalized_counts[, top_genes, drop = FALSE]
      }
      
      # Use vectorized correlation for better performance
      correlation_df <- sapply(names(normalized_counts), function(x) {
        if(input$test_choice == "Spearman") {
          test_result <- cor.test(first_gene_count, normalized_counts[[x]], method = "spearman")
        } else {
          test_result <- cor.test(first_gene_count, normalized_counts[[x]], method = "pearson")
        }
        return(c(correlation = test_result$estimate, p_value = test_result$p.value))
      })
      
      correlation_df <- as.data.frame(t(correlation_df))
      num_tests <- ncol(normalized_counts)
      correlation_df$Bonferroni_p_value <- pmin(1, correlation_df$p_value * num_tests)
      correlation_df$Gene <- rownames(correlation_df)
      colnames(correlation_df)[1:2] <- c("Correlation","p-val")
      correlation_df <- correlation_df[, c("Gene", "Correlation", "p-val", "Bonferroni_p_value")]
      
      # Sort by absolute correlation value
      correlation_df <- correlation_df[order(abs(correlation_df$Correlation), decreasing = TRUE), ]
      
      return(correlation_df)}else{
        normalized_counts <- as.matrix(filtered_adata()$X)
        normalized_counts <- as.data.frame(normalized_counts)
        colnames(normalized_counts) <- gene_list_adata()
        first_gene_count <- normalized_counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_first)]
        
        # Optimization: Limit to most variable genes for faster computation
        if(ncol(normalized_counts) > 1000) {
          gene_vars <- apply(normalized_counts, 2, var, na.rm = TRUE)
          top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:1000])
          normalized_counts <- normalized_counts[, top_genes, drop = FALSE]
        }
        
        correlation_df <- sapply(names(normalized_counts), function(x) {
          if(input$test_choice == "Spearman") {
            test_result <- cor.test(first_gene_count, normalized_counts[[x]], method = "spearman")
          } else {
            test_result <- cor.test(first_gene_count, normalized_counts[[x]], method = "pearson")
          }
          return(c(correlation = test_result$estimate, p_value = test_result$p.value))
        })
        
        correlation_df <- as.data.frame(t(correlation_df))
        num_tests <- ncol(normalized_counts)
        correlation_df$Bonferroni_p_value <- pmin(1, correlation_df$p_value * num_tests)
        correlation_df$Gene <- rownames(correlation_df)
        colnames(correlation_df)[1:2] <- c("Correlation","p-val")
        correlation_df <- correlation_df[, c("Gene", "Correlation", "p-val", "Bonferroni_p_value")]
        
        # Sort by absolute correlation value
        correlation_df <- correlation_df[order(abs(correlation_df$Correlation), decreasing = TRUE), ]
        
        return(correlation_df)
        
      }
    })
    

    output$first_gene_correlation_table <- renderDT({
      req(correlation_table_first_gene())
      datatable(
        correlation_table_first_gene(), 
        options = list(
          pageLength = 5,
          autoWidth = TRUE,
          scrollX = TRUE),
        filter = 'top',
        rownames = FALSE,
        extensions = 'Scroller'
      )%>%
        formatRound(columns=c("Correlation","p-val","Bonferroni_p_value"), digits=3)
      
    })
    
    
    correlation_table_second_gene <- eventReactive(input$top_correlated_second_gene,{
      
      if(is.null(input$filter_dataset_cluster_selection)){
        
        # Standard correlation analysis
        normalized_counts <- as.matrix(adata()$X)
        normalized_counts <- as.data.frame(normalized_counts)
        colnames(normalized_counts) <- gene_list_adata()
        second_gene_count <- normalized_counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_second)]
      
      # Optimization: Limit to most variable genes for faster computation
      if(ncol(normalized_counts) > 1000) {
        gene_vars <- apply(normalized_counts, 2, var, na.rm = TRUE)
        top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:1000])
        normalized_counts <- normalized_counts[, top_genes, drop = FALSE]
      }
      
      correlation_df <- sapply(names(normalized_counts), function(x) {
        if(input$test_choice == "Spearman") {
          test_result <- cor.test(second_gene_count, normalized_counts[[x]], method = "spearman")
        } else {
          test_result <- cor.test(second_gene_count, normalized_counts[[x]], method = "pearson")
        }
        return(c(correlation = test_result$estimate, p_value = test_result$p.value))
      })
      
      correlation_df <- as.data.frame(t(correlation_df))
      num_tests <- ncol(normalized_counts)
      correlation_df$Bonferroni_p_value <- pmin(1, correlation_df$p_value * num_tests)
      correlation_df$Gene <- rownames(correlation_df)
      colnames(correlation_df)[1:2] <- c("Correlation","p-val")
      correlation_df <- correlation_df[, c("Gene", "Correlation", "p-val", "Bonferroni_p_value")]
      
      # Sort by absolute correlation value
      correlation_df <- correlation_df[order(abs(correlation_df$Correlation), decreasing = TRUE), ]
      
      return(correlation_df)
      }else{
        normalized_counts <- as.matrix(filtered_adata()$X)
        normalized_counts <- as.data.frame(normalized_counts)
        colnames(normalized_counts) <- gene_list_adata()
        second_gene_count <- normalized_counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_second)]
        
        # Optimization: Limit to most variable genes for faster computation
        if(ncol(normalized_counts) > 1000) {
          gene_vars <- apply(normalized_counts, 2, var, na.rm = TRUE)
          top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:1000])
          normalized_counts <- normalized_counts[, top_genes, drop = FALSE]
        }
        
        correlation_df <- sapply(names(normalized_counts), function(x) {
          if(input$test_choice == "Spearman") {
            test_result <- cor.test(second_gene_count, normalized_counts[[x]], method = "spearman")
          } else {
            test_result <- cor.test(second_gene_count, normalized_counts[[x]], method = "pearson")
          }
          return(c(correlation = test_result$estimate, p_value = test_result$p.value))
        })
        
        correlation_df <- as.data.frame(t(correlation_df))
        num_tests <- ncol(normalized_counts)
        correlation_df$Bonferroni_p_value <- pmin(1, correlation_df$p_value * num_tests)
        correlation_df$Gene <- rownames(correlation_df)
        colnames(correlation_df)[1:2] <- c("Correlation","p-val")
        correlation_df <- correlation_df[, c("Gene", "Correlation", "p-val", "Bonferroni_p_value")]
        
        # Sort by absolute correlation value
        correlation_df <- correlation_df[order(abs(correlation_df$Correlation), decreasing = TRUE), ]
        
        return(correlation_df)
      }
    })
    
    
    output$second_gene_correlation_table <- renderDT({
      req(correlation_table_second_gene())
      datatable(
        correlation_table_second_gene(), 
        options = list(
          pageLength = 5,
          autoWidth = TRUE,
          scrollX = TRUE),
        filter = 'top',
        rownames = FALSE,
        extensions = 'Scroller'
      ) %>%
        formatRound(columns=c("Correlation","p-val","Bonferroni_p_value"), digits=3)
      
    })
    
    ##########################################################
    
    ######### Differential Expression
    
    output$de_ident_selection_name_first <- renderUI({
      req(adata())
      if(input$de_data == "All Data" && input$de_type == "Clusters"){
        textInput("de_ident_1_name", "Name First Ident", placeholder = "Enter Name")
      }else if(input$de_data == "All Data" && input$de_type == "Groups"){
        textInput("de_ident_1_name", "Name First Ident", placeholder = "Enter Name")
      }else if(input$de_data == "Filtered Data" && input$de_type == "Clusters"){
        filtered_available <- tryCatch({
          !is.null(filtered_adata())
        }, error = function(e) FALSE)
        if(filtered_available){
          textInput("de_ident_1_name", "Name First Ident", placeholder = "Enter Name")
        }else{
          "To continue, first filter the data in the Cluster Selection tab."
        }
      }else if(input$de_data == "Filtered Data" && input$de_type == "Groups"){
        filtered_available <- tryCatch({
          !is.null(filtered_adata())
        }, error = function(e) FALSE)
        if(filtered_available){
          textInput("de_ident_1_name", "Name First Ident", placeholder = "Enter Name")
        }else{
          "To continue, first filter the data in the Cluster Selection tab."
        }
      }
    })
    
    
    output$de_ident_selection_first <- renderUI({
      req(adata(),input$de_ident_1_name)
      if(input$de_data == "All Data" && input$de_type == "Clusters"){
        selectizeInput("de_ident_1", "First Ident", choices = sort(unique(adata()$obs$CellType)), multiple = T)
      }else if(input$de_data == "All Data" && input$de_type == "Groups"){
        selectizeInput("de_ident_1", "First Ident", choices = sort(unique(adata()$obs$Group)), multiple = T)
      }else if(input$de_data == "Filtered Data" && input$de_type == "Clusters"){
        if(!is.null(filtered_adata())){
          selectizeInput("de_ident_1", "First Ident", choices = sort(unique(filtered_adata()$obs$CellType)), multiple = T)
        }else if(is.null(filtered_adata())){
          "You Should First Filter Data"
        }
      }else if(input$de_data == "Filtered Data" && input$de_type == "Groups"){
        if(!is.null(filtered_adata())){
          selectizeInput("de_ident_1", "First Ident", choices = sort(unique(filtered_adata()$obs$Group)), multiple = T)
        }else{
          "You Should First Filter Data"
        }
      }
    })
    
    
    output$de_ident_selection_name_second <- renderUI({
      req(adata(), input$de_ident_1)
      if(input$de_data == "All Data" && input$de_type == "Clusters"){
        textInput("de_ident_2_name", "Name Second Ident", placeholder = "Enter Name")
      }else if(input$de_data == "All Data" && input$de_type == "Groups"){
        textInput("de_ident_2_name", "Name Second Ident", placeholder = "Enter Name")
      }else if(input$de_data == "Filtered Data" && input$de_type == "Clusters"){
        if(!is.null(filtered_adata())){
          textInput("de_ident_2_name", "Name Second Ident", placeholder = "Enter Name")
        }else if(is.null(filtered_adata())){
          "You Should First Filter Data"
        }
      }else if(input$de_data == "Filtered Data" && input$de_type == "Groups"){
        if(!is.null(filtered_adata())){
          textInput("de_ident_2_name", "Name Second Ident", placeholder = "Enter Name")
        }else{
          "You Should First Filter Data"
        }
      }
    })
    
    
    
    output$de_ident_selection_second <- renderUI({
      req(adata(),input$de_ident_2_name)
      if(input$de_data == "All Data" && input$de_type == "Clusters"){
        first_ident <- input$de_ident_1
        idents <- unique(adata()$obs$CellType)
        first_ident <- idents %in% first_ident
        idents <- idents[!first_ident]
        selectizeInput("de_ident_2", "Second Ident", choices = sort(idents),multiple = T)
      }else if(input$de_data == "All Data" && input$de_type == "Groups"){
        first_ident <- input$de_ident_1
        idents <- unique(adata()$obs$Group)
        first_ident <- idents %in% first_ident
        idents <- idents[!first_ident]
        selectizeInput("de_ident_2", "Second Ident", choices = sort(idents),multiple = T)
      }else if(input$de_data == "Filtered Data" && input$de_type == "Clusters"){
        first_ident <- input$de_ident_1
        idents <- unique(filtered_adata()$obs$CellType)
        first_ident <- idents %in% first_ident
        idents <- idents[!first_ident]
        selectizeInput("de_ident_2", "Second Ident", choices = sort(idents),multiple = T)
      }else if(input$de_data == "Filtered Data" && input$de_type == "Groups"){
        first_ident <- input$de_ident_1
        idents <- unique(filtered_adata()$obs$Group)
        first_ident <- idents %in% first_ident
        idents <- idents[!first_ident]
        selectizeInput("de_ident_2", "Second Ident", choices = sort(idents),multiple = T)
      }
    })
    
    output$de_method_selection <- renderUI({
      req(input$de_ident_2)
      selectInput("de_method_selection", "Select Method", choices = c('t-test', 'wilcoxon', 't-test_overestim_var'), selected = "wilcoxon")
    })
    
    output$de_run_analysis <- renderUI({
      req(input$de_method_selection)
      actionButton("de_run_dge", "Run DGE", class = "btn-primary", icon = icon("play"), width = '100%')
    })
    
    
    de_filtering_identifying_clusters <- eventReactive(input$de_run_dge,{
      
      if(input$de_data == "All Data" && input$de_type == "Clusters"){
        adata <- adata()
        obs_data <- adata$obs
        obs_data$ident <- "ident"
        first_identity <- input$de_ident_1
        second_identity <- input$de_ident_2
        first_identity_rows <- obs_data$CellType %in% first_identity
        second_identity_rows <- obs_data$CellType %in% second_identity
        obs_data$ident[first_identity_rows] <- input$de_ident_1_name 
        obs_data$ident[second_identity_rows] <- input$de_ident_2_name       
        adata$obs['ident'] = obs_data$ident
        selection_idents <- list(input$de_ident_1_name, input$de_ident_2_name)
        adata = adata[adata$obs$ident %in% selection_idents, ]

        return(adata)
        

      }else if(input$de_data == "All Data" && input$de_type == "Groups"){
        adata <- adata()
        obs_data <- adata$obs
        obs_data$ident <- "ident"
        first_identity <- input$de_ident_1
        second_identity <- input$de_ident_2
        first_identity_rows <- obs_data$Group %in% first_identity
        second_identity_rows <- obs_data$Group %in% second_identity
        obs_data$ident[first_identity_rows] <- input$de_ident_1_name 
        obs_data$ident[second_identity_rows] <- input$de_ident_2_name       
        adata$obs$ident = obs_data$ident
        selection_idents <- c(input$de_ident_1_name, input$de_ident_2_name)
        adata = adata[adata$obs$ident %in% selection_idents, ]
        return(adata)
        
      }else if(input$de_data == "Filtered Data" && input$de_type == "Clusters"){
        adata <- filtered_adata()
        obs_data <- adata$obs
        obs_data$ident <- "ident"
        first_identity <- input$de_ident_1
        second_identity <- input$de_ident_2
        first_identity_rows <- obs_data$CellType %in% first_identity
        second_identity_rows <- obs_data$CellType %in% second_identity
        obs_data$ident[first_identity_rows] <- input$de_ident_1_name 
        obs_data$ident[second_identity_rows] <- input$de_ident_2_name       
        adata$obs['ident'] = obs_data$ident
        selection_idents <- c(input$de_ident_1_name, input$de_ident_2_name)
        adata = adata[adata$obs$ident %in% selection_idents, ]
        return(adata)
        
      }else if(input$de_data == "Filtered Data" && input$de_type == "Groups"){
        adata <- filtered_adata()
        obs_data <- adata$obs
        obs_data$ident <- "ident"
        first_identity <- input$de_ident_1
        second_identity <- input$de_ident_2
        first_identity_rows <- obs_data$Group %in% first_identity
        second_identity_rows <- obs_data$Group %in% second_identity
        obs_data$ident[first_identity_rows] <- input$de_ident_1_name 
        obs_data$ident[second_identity_rows] <- input$de_ident_2_name       
        adata$obs['ident'] = obs_data$ident
        selection_idents <- c(input$de_ident_1_name, input$de_ident_2_name)
        adata = adata[adata$obs$ident %in% selection_idents, ]
        
        return(adata)
      }
    })
    
    
    de_dge_calculation <- reactive({
      req(de_filtering_identifying_clusters())
        adata <- de_filtering_identifying_clusters()
        sc$tl$rank_genes_groups(adata, 'ident', groups=list(input$de_ident_1_name), reference= input$de_ident_2_name, method= input$de_method_selection,pts = T)
        return(adata)

    })
    
    
    output$imageoutput_dge_ranks <- renderImage({
      req(de_dge_calculation())
      group_name <- list(input$de_ident_1_name)
      sc$pl$rank_genes_groups(de_dge_calculation(), groups = group_name, sharey = FALSE, show=FALSE, save = 'dge.png')
      list(src = "figures/rank_genes_groups_identdge.png")
    }, deleteFile = TRUE)
  
    output$dge_dt <- renderDT({
      req(de_dge_calculation())
      group_name <- list(input$de_ident_1_name)
      result_df = sc$get$rank_genes_groups_df(de_dge_calculation(), group = group_name, key = 'rank_genes_groups')
      
      colnames(result_df) <- c("Gene", "Scores", "LogFC", "p-val", "adj-p", "pct")
      
      if(input$dge_filter_table_logfc > 0){
        result_df <- result_df %>%
          filter(LogFC > 2)
      }
      
      if(input$dge_filter_table_logfc_neg > 0){
        result_df <- result_df %>%
          filter(LogFC < -2)
      }
      
      if(input$dge_filter_table_score_neg > 0){
        result_df <- result_df %>%
          filter(Scores < -2)
      }
      
      if(input$dge_filter_table_score > 0){
        result_df <- result_df %>%
          filter(Scores > 2)
      }
      
      if(input$dge_filter_table_pval > 0){
        result_df <- result_df %>%
          filter(`adj-p` < 0.05)
      }
      
      if(input$dge_filter_table_reset > 0){
          result_df = sc$get$rank_genes_groups_df(de_dge_calculation(), group = group_name, key = 'rank_genes_groups')
          
          colnames(result_df) <- c("Gene", "Scores", "LogFC", "p-val", "adj-p", "pct")
          
      }
      
      datatable(
        result_df, 
        options = list(
          pageLength = 5,
          autoWidth = TRUE,
          scrollX = TRUE,
          dom = 'l<"toolbar">rtip' 
        ),
        filter = 'top',
        rownames = FALSE,
        selection = "single",
        extensions = 'Scroller'
      ) %>%
        formatRound(which(sapply(result_df, is.numeric)), digits = 3)
    })
    
    output$dge_table_filter <- renderUI({
      req(de_dge_calculation())
      fluidRow(
        column(width = 2, actionButton("dge_filter_table_logfc", label = "LogFC > 2", class = "btn-primary", width = '100%')),
        column(width = 2, actionButton("dge_filter_table_logfc_neg", label = "LogFC < -2", class = "btn-primary", width = '100%')),
        column(width = 2, actionButton("dge_filter_table_pval", label = "Adjusted-P Val < 0.05", class = "btn-primary", width = '100%')),
        column(width = 2, actionButton("dge_filter_table_reset", label = "Reset Filters", class = "btn-primary", width = '100%')),
        column(width = 2, actionButton("dge_filter_table_score", label = "Score > 2", class = "btn-primary", width = '100%')),
        column(width = 2, actionButton("dge_filter_table_score_neg", label = "Score < -2", class = "btn-primary", width = '100%'))
      )
    })
    
    
    
    
    output$imageoutput_dge_violin <- renderImage({
      req(de_dge_calculation())
      group_name <- list(input$de_ident_1_name)
      sc$pl$rank_genes_groups_violin(de_dge_calculation(), groups = group_name,show=FALSE, save = '.png')

      list(src = paste0("figures/rank_genes_groups_ident_",group_name,".png"))
      
    }, deleteFile = TRUE)
    
    
    output$radiobutton_dge <- renderUI({
      req(de_dge_calculation())
      fluidRow(column(width = 6, actionButton("dge_violin","Violin Plot", icon = icon("chart-simple"),class = "btn-primary", width = '100%')),
               column(width = 6, actionButton("dge_go","Go Enrichment", icon = icon("chart-bar"),class = "btn-primary", width = '100%')))
    })
    
    output$violin_dge <- renderImage({
      req(de_dge_calculation(), input$dge_dt_rows_selected)

      group_name <- list(input$de_ident_1_name)
      result_df = sc$get$rank_genes_groups_df(de_dge_calculation(), group = group_name, key = 'rank_genes_groups')
      selected_row_index <- input$dge_dt_rows_selected
      selected_row <- result_df[selected_row_index, ]
      selected_value <- selected_row[[1]]
      
      sc$pl$violin(de_dge_calculation(), keys = selected_value, groupby = 'ident', show=FALSE, rotation=90, save = "violin_dge.png")
      list(src = "figures/violinviolin_dge.png")
    }, deleteFile = TRUE)
    
    
    de_enrichment_calc <- reactive({
      # Check if fenr package is available
      if (!requireNamespace("fenr", quietly = TRUE)) {
        showNotification("❌ Package 'fenr' not available. Enrichment analysis is disabled.", 
                        type = "error", duration = 5)
        return(NULL)
      }
      
      req(de_dge_calculation(), input$dge_dt_rows_all)
      group_name <- list(input$de_ident_1_name)
      result_df = sc$get$rank_genes_groups_df(de_dge_calculation(), group = group_name, key = 'rank_genes_groups')

      if(input$selection_organism %in% c("Human", "Integrated")){
        load("enrichment_sets/human.RData")
        all_genes <- de_dge_calculation()$var
        all_genes <- rownames(all_genes)
        selected_genes <- result_df[input$dge_dt_rows_all, ]
        selected_genes <- selected_genes[,1]
        
        
        bp_human <- prepare_for_enrichment(bp_human$terms, bp_human$mapping, all_genes, feature_name = "gene_symbol")
        go_human <- prepare_for_enrichment(go_human$terms, go_human$mapping, all_genes, feature_name = "gene_symbol")
        kegg_human <- prepare_for_enrichment(kegg_human$terms, kegg_human$mapping, all_genes, feature_name = "gene_symbol")
        reactome_human <- prepare_for_enrichment(reactome_human$terms, reactome_human$mapping, all_genes, feature_name = "gene_symbol")
        wiki_pathways_human <- prepare_for_enrichment(wiki_pathways_human$terms, wiki_pathways_human$mapping, all_genes, feature_name = "text_label")
        
        bp_human <- functional_enrichment(all_genes, selected_genes, bp_human)
        go_human <- functional_enrichment(all_genes, selected_genes, go_human)
        kegg_human <- functional_enrichment(all_genes, selected_genes, kegg_human)
        reactome_human <- functional_enrichment(all_genes, selected_genes, reactome_human)
        wiki_pathways_human <- functional_enrichment(all_genes, selected_genes, wiki_pathways_human)
        
        de_enrichment_calc <- list(bp_human, go_human, kegg_human, reactome_human, wiki_pathways_human)
        
        return(de_enrichment_calc)
        
      }else if(input$selection_organism == "Mouse"){
        
        load("enrichment_sets/mouse.RData")
        all_genes <- de_dge_calculation()$var
        all_genes <- rownames(all_genes)
        selected_genes <- result_df[input$dge_dt_rows_all, ]
        selected_genes <- selected_genes[,1]
        
        bp_mouse <- prepare_for_enrichment(bp_mouse$terms, bp_mouse$mapping, all_genes, feature_name = "gene_symbol")
        go_mouse <- prepare_for_enrichment(go_mouse$terms, go_mouse$mapping, all_genes, feature_name = "gene_symbol")
        kegg_mouse <- prepare_for_enrichment(kegg_mouse$terms, kegg_mouse$mapping, all_genes, feature_name = "gene_symbol")
        reactome_mouse <- prepare_for_enrichment(reactome_mouse$terms, reactome_mouse$mapping, all_genes, feature_name = "gene_symbol")
        wiki_pathways_mouse <- prepare_for_enrichment(wiki_pathways_mouse$terms, wiki_pathways_mouse$mapping, all_genes, feature_name = "text_label")
        
        bp_mouse <- functional_enrichment(all_genes, selected_genes, bp_mouse)
        go_mouse <- functional_enrichment(all_genes, selected_genes, go_mouse)
        kegg_mouse <- functional_enrichment(all_genes, selected_genes, kegg_mouse)
        reactome_mouse <- functional_enrichment(all_genes, selected_genes, reactome_mouse)
        wiki_pathways_mouse <- functional_enrichment(all_genes, selected_genes, wiki_pathways_mouse)
        de_enrichment_calc <- list(bp_mouse, go_mouse, kegg_mouse, reactome_mouse, wiki_pathways_mouse)
        return(de_enrichment_calc)
      }else if(input$selection_organism == "Zebrafish"){
        
        load("enrichment_sets/zebrafish.RData")
        all_genes <- de_dge_calculation()$var
        all_genes <- rownames(all_genes)
        selected_genes <- result_df[input$dge_dt_rows_all, ]
        selected_genes <- selected_genes[,1]
        
        bp_zebrafish <- prepare_for_enrichment(bp_zebrafish$terms, bp_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        go_zebrafish <- prepare_for_enrichment(go_zebrafish$terms, go_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        kegg_zebrafish <- prepare_for_enrichment(kegg_zebrafish$terms, kegg_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        reactome_zebrafish <- prepare_for_enrichment(reactome_zebrafish$terms, reactome_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        wiki_pathways_zebrafish <- prepare_for_enrichment(wiki_pathways_zebrafish$terms, wiki_pathways_zebrafish$mapping, all_genes, feature_name = "text_label")
        
        bp_zebrafish <- functional_enrichment(all_genes, selected_genes, bp_zebrafish)
        go_zebrafish <- functional_enrichment(all_genes, selected_genes, go_zebrafish)
        kegg_zebrafish <- functional_enrichment(all_genes, selected_genes, kegg_zebrafish)
        reactome_zebrafish <- functional_enrichment(all_genes, selected_genes, reactome_zebrafish)
        wiki_pathways_zebrafish <- functional_enrichment(all_genes, selected_genes, wiki_pathways_zebrafish)
        de_enrichment_calc <- list(bp_zebrafish, go_zebrafish, kegg_zebrafish, reactome_zebrafish, wiki_pathways_zebrafish)
        return(de_enrichment_calc)
        
      }
      # else if(input$selection_organism == "Integrated"){
      #   load("enrichment_sets/human.RData")
      #   all_genes <- de_dge_calculation()$var
      #   all_genes <- rownames(all_genes)
      #   selected_genes <- result_df[input$dge_dt_rows_all, ]
      #   
      #   bp_human <- prepare_for_enrichment(bp_human$terms, bp_human$mapping, all_genes, feature_name = "gene_symbol")
      #   go_human <- prepare_for_enrichment(go_human$terms, go_human$mapping, all_genes, feature_name = "gene_symbol")
      #   kegg_human <- prepare_for_enrichment(kegg_human$terms, kegg_human$mapping, all_genes, feature_name = "gene_symbol")
      #   reactome_human <- prepare_for_enrichment(reactome_human$terms, reactome_human$mapping, all_genes, feature_name = "gene_symbol")
      #   wiki_pathways_human <- prepare_for_enrichment(wiki_pathways_human$terms, wiki_pathways_human$mapping, all_genes, feature_name = "gene_symbol")
      #   
      #   bp_human <- functional_enrichment(all_genes, selected_genes, bp_human)
      #   go_human <- functional_enrichment(all_genes, selected_genes, go_human)
      #   kegg_human <- functional_enrichment(all_genes, selected_genes, kegg_human)
      #   reactome_human <- functional_enrichment(all_genes, selected_genes, reactome_human)
      #   wiki_pathways_human <- functional_enrichment(all_genes, selected_genes, wiki_pathways_human)
      #   de_enrichment_calc <- list(bp_human, go_human, kegg_human, reactome_human, wiki_pathways_human)
      #   return(de_enrichment_calc)}
    })

    output$de_enrichment_table <- renderDT({
      if (is.null(de_enrichment_calc())) {
        # Return empty datatable with message
        return(datatable(
          data.frame(Message = "Enrichment analysis not available. Please install the 'fenr' package."),
          options = list(pageLength = 5, dom = 't'),
          rownames = FALSE
        ))
      }
      
      req(de_enrichment_calc())
      if(input$de_enrichment_type == "GO"){
        datatable(
          de_enrichment_calc()[[2]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollY = 300,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip',
            columnDefs = list(
                list(
                  targets = which(names(de_enrichment_calc()[[2]]) == "ids") - 1,  # JS indexing starts at 0
                  render = JS(
                    "function(data, type, row, meta) {",
                    "  if (type === 'display') {",
                    "    var displayText = data.length > 50 ? data.substr(0, 50) + '???' : data;",
                    "    return '<span title=\"' + data + '\">' + displayText + '</span>';",
                    "  }",
                    "  return data;",
                    "}"
                  )
                )
              )
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(de_enrichment_calc()[[2]], is.numeric)), digits = 3)
      }else if(input$de_enrichment_type == "BP"){
        datatable(
          de_enrichment_calc()[[1]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollY = 300,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip',
            columnDefs = list(
              list(
                targets = which(names(de_enrichment_calc()[[2]]) == "ids") - 1,  # JS indexing starts at 0
                render = JS(
                  "function(data, type, row, meta) {",
                  "  if (type === 'display') {",
                  "    var displayText = data.length > 50 ? data.substr(0, 50) + '???' : data;",
                  "    return '<span title=\"' + data + '\">' + displayText + '</span>';",
                  "  }",
                  "  return data;",
                  "}"
                )
              )
            )
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(de_enrichment_calc()[[1]], is.numeric)), digits = 3)
      }else if(input$de_enrichment_type == "KEGG"){
        datatable(
          de_enrichment_calc()[[3]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollY = 300,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip',
            columnDefs = list(
              list(
                targets = which(names(de_enrichment_calc()[[2]]) == "ids") - 1,  # JS indexing starts at 0
                render = JS(
                  "function(data, type, row, meta) {",
                  "  if (type === 'display') {",
                  "    var displayText = data.length > 50 ? data.substr(0, 50) + '???' : data;",
                  "    return '<span title=\"' + data + '\">' + displayText + '</span>';",
                  "  }",
                  "  return data;",
                  "}"
                )
              )
            )
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(de_enrichment_calc()[[3]], is.numeric)), digits = 3)
      }else if(input$de_enrichment_type == "Reactome"){
        datatable(
          de_enrichment_calc()[[4]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollY = 300,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip',
            columnDefs = list(
              list(
                targets = which(names(de_enrichment_calc()[[2]]) == "ids") - 1,  # JS indexing starts at 0
                render = JS(
                  "function(data, type, row, meta) {",
                  "  if (type === 'display') {",
                  "    var displayText = data.length > 50 ? data.substr(0, 50) + '???' : data;",
                  "    return '<span title=\"' + data + '\">' + displayText + '</span>';",
                  "  }",
                  "  return data;",
                  "}"
                )
              )
            )
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(de_enrichment_calc()[[4]], is.numeric)), digits = 3)
      }else if(input$de_enrichment_type == "WikiPathways"){
        datatable(
          de_enrichment_calc()[[5]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollY = 300,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip',
            columnDefs = list(
              list(
                targets = which(names(de_enrichment_calc()[[2]]) == "ids") - 1,  # JS indexing starts at 0
                render = JS(
                  "function(data, type, row, meta) {",
                  "  if (type === 'display') {",
                  "    var displayText = data.length > 50 ? data.substr(0, 50) + '???' : data;",
                  "    return '<span title=\"' + data + '\">' + displayText + '</span>';",
                  "  }",
                  "  return data;",
                  "}"
                )
              )
            )
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(de_enrichment_calc()[[5]], is.numeric)), digits = 3)
      }
    })
    
  
    # selectizeInput("de_ident_2", "Second Ident", choices = sort(unique(adata()$obs$CellType)),multiple = T)
    
    
    ########## pseudo_bulk
    
    output$pseudo_ident_selection_first <- renderUI({
      req(adata())
      if(input$pseudo_bulk_data == "All Data"){
        selectizeInput("pseudo_ident_1", "First Ident", choices = sort(unique(adata()$obs$Group)), multiple = T)
      }else if(input$pseudo_bulk_data == "Filtered Data"){
        filtered_available <- tryCatch({
          !is.null(filtered_adata())
        }, error = function(e) FALSE)
        if(filtered_available){
          selectizeInput("pseudo_ident_1", "First Ident", choices = sort(unique(filtered_adata()$obs$Group)), multiple = T)
        }else{
          "To continue, first filter the data in the Cluster Selection tab."
        }
      }
    })
    
    output$pseudo_ident_selection_second <- renderUI({
      req(adata(),input$pseudo_ident_1)
      if(input$pseudo_bulk_data == "All Data"){
        idents <- unique(adata()$obs$Group)
        first_ident <- idents %in% input$pseudo_ident_1
        idents <- idents[!first_ident]
        selectizeInput("pseudo_ident_2", "Second Ident", choices = sort(unique(idents)), multiple = T)
      }else if(input$pseudo_bulk_data == "Filtered Data"){
        idents <- unique(filtered_adata()$obs$Group)
        first_ident <- idents %in% input$pseudo_ident_1
        idents <- idents[!first_ident]
        selectizeInput("pseudo_ident_2", "Second Ident", choices = sort(unique(selectizeInput)), multiple = T)
      }
    })
    
    output$run_pseudo_bulk <- renderUI({
      req(input$pseudo_ident_2)
      actionButton("run_pseudo_bulk_act", "Run Pseudo Bulk", class = "btn-primary",width = '100%')
    })
    
    output$pseudo_pca_selection_ui  <- renderUI({
      req(input$run_pseudo_bulk_act)
      fluidRow(column(width = 6, selectInput("pseudo_pca_selection", label = NULL, choices = c("ident", "CellType"))),
               column(width = 6, actionButton("associations_pseudo_modal", "Show PC Associations", class = "btn-primary", width = '100%')))
    })
    
    
    
    pdata <- eventReactive(input$run_pseudo_bulk_act, {
      if(input$pseudo_bulk_data == "All Data"){
        adata <- adata()
        obs_data <- adata$obs
        obs_data$ident <- "ident"
        first_identity <- input$pseudo_ident_1
        second_identity <- input$pseudo_ident_2
        first_identity_rows <- obs_data$Group %in% first_identity
        second_identity_rows <- obs_data$Group %in% second_identity
        obs_data$ident[first_identity_rows] <- 'Selection' 
        obs_data$ident[second_identity_rows] <- 'Reference'       
        adata$obs['ident'] = obs_data$ident
        selection_idents <- list('Selection', 'Reference' )
        adata = adata[adata$obs$ident %in% selection_idents, ]  
      }else if(input$pseudo_bulk_data == "Filtered Data"){
        adata <- filtered_adata()
        obs_data <- adata$obs
        obs_data$ident <- "ident"
        first_identity <- input$pseudo_ident_1
        second_identity <- input$pseudo_ident_2
        first_identity_rows <- obs_data$Group %in% first_identity
        second_identity_rows <- obs_data$Group %in% second_identity
        obs_data$ident[first_identity_rows] <- 'Selection' 
        obs_data$ident[second_identity_rows] <- 'Reference'       
        adata$obs['ident'] = obs_data$ident
        selection_idents <- list('Selection', 'Reference' )
        adata = adata[adata$obs$ident %in% selection_idents, ]
      }
      
      pdata <- dc$get_pseudobulk(
        adata,
        sample_col = 'Group',
        groups_col = 'CellType',
        layer = 'counts',
        mode = 'sum',
        min_cells = as.integer(10),
        min_counts = as.integer(1000)
      )
      pdata$layers['counts'] <- pdata$X
      sc$pp$normalize_total(pdata, target_sum = 1e4)
      sc$pp$log1p(pdata)
      sc$pp$scale(pdata, max_value = 10)
      sc$tl$pca(pdata)
      dc$swap_layer(pdata, 'counts', X_layer_key = "None", inplace = TRUE)
      dc$get_metadata_associations(
        pdata,
        obs_keys = list('ident', 'CellType', 'psbulk_n_cells', 'psbulk_counts'),
        obsm_key = 'X_pca',
        uns_key = 'pca_anova',
        inplace = TRUE
      )
      return(pdata)
    })
    
    output$pca_pseudo_bulk <- renderImage({
      req(pdata())
      if(input$pseudo_pca_selection == "ident"){
        sc$pl$pca(pdata(), color = c('ident'), show = F, save = "pseudo_pca.png")
        list(src = "figures/pcapseudo_pca.png")
      }else if(input$pseudo_pca_selection == "CellType"){
        sc$pl$pca(pdata(), color = c('CellType'), show = F, save = "pseudo_pca.png")
        list(src = "figures/pcapseudo_pca.png")
      }else if(is.null(input$pseudo_pca_selection)){
        sc$pl$pca(pdata(), color = c('ident'), show = F, save = "pseudo_pca.png")
        list(src = "figures/pcapseudo_pca.png")
      }
    }, deleteFile = TRUE)
    
    output$pca_pseudo_bulk_associations <- renderImage({
      req(pdata())
      dc$plot_associations(
        pdata(),
        uns_key = 'pca_anova',
        obsm_key = 'X_pca',
        stat_col = 'p_adj',
        obs_annotation_cols = list('ident', 'CellType'),
        titles = c('Principle component scores', 'Adjusted p-values from ANOVA'),
        figsize = c(as.integer(7), as.integer(5)),
        n_factors = as.integer(10),
        cmap_stats= as.character('Purples'), 
        cmap_scores=as.character('BrBG'), 
        cmap_cats= as.character('Set2'), 
        save = "adjusted_pca.png",
      )
      list(src = "adjusted_pca.png")
    }, deleteFile = TRUE)
    
    
    observeEvent(input$associations_pseudo_modal, {
      showModal(modalDialog(
        div(
          style = "height: 500px; overflow-y: auto;",
          imageOutput("pca_pseudo_bulk_associations")
        ),
        size = "l", 
        easyClose = TRUE
      ))
    })
    
    output$pdata_clusters <- renderUI({
      req(pdata())
      selectizeInput("pdata_clusters_filter", "Select Clusters to Filter", choices = sort(unique(pdata()$obs$CellType)),
                     multiple = T)
    })
    
    output$run_deseq2_ui <- renderUI({
      req(input$pdata_clusters_filter)
      actionButton("run_deseq2", "Run DESeq2", class = "btn-primary",width = '100%')
    })
    
    results_df <- eventReactive(input$run_deseq2,{
      
      filtered_pseudo <- pdata()[pdata()$obs$CellType %in% input$pdata_clusters_filter, ]
      genes <- dc$filter_by_expr(filtered_pseudo, group = 'ident', min_count = 10, min_total_count = 15)
      filtered_pseudo <- filtered_pseudo[, genes]
      filtered_pseudo <- filtered_pseudo$copy()
      
      group_counts <- table(filtered_pseudo$obs$ident)
      if (any(group_counts < 2)) {
        showNotification("Each group in 'ident' must have at least 2 samples for DESeq2 to work.", type = "error")
        return(NULL)
      }
      
      dds <- pydeseq2_dds$DeseqDataSet(
        adata = filtered_pseudo,
        design_factors = 'ident',
        ref_level=list('ident', 'Reference'),
        refit_cooks = TRUE
      )
      
      dds$deseq2()
      
      stat_res = pydeseq2_ds$DeseqStats(
        dds,
        contrast=list("ident", 'Selection', 'Reference')
      )
      
      
      stat_res$summary()
      
      results_df = stat_res$results_df
      results_df <- as.data.frame(results_df)
      results_df[["Gene Name"]] <- rownames(results_df)
      results_df <- results_df %>% relocate(`Gene Name`, .before = 1)
      mat <- as.data.frame(t(results_df['stat']))
      rownames(mat) <- "Selected_Clusters"
      results_df <- list(results_df, mat)
      return(results_df)
      
    })
    
    output$selection_pseudo_volcano_table_ui <- renderUI({
      req(results_df()[[1]][[1]])
      selectInput("selection_pseudo_volcano_table", label=NULL, choices = c("Result Table", "Volcano"))
    })
    
    
    output$pca_pseudo_bulk_volcano <- renderImage({
      req(results_df()[[1]][[1]])
      
      dataset <- results_df()[[1]]
      dataset[is.na(dataset)] <- 0.000001
      dataset[dataset == ""] <- 0.000001
      
      dc$plot_volcano_df(
        results_df()[[1]],
        x='log2FoldChange',
        y='padj',
        top= as.integer(20),
        save = "pseudo_volcano.png"
      )
      list(src = "pseudo_volcano.png")
    }, deleteFile = TRUE)
    
    output$pca_pseudo_bulk_results_table <- renderDT({
      req(results_df()[[1]])
      dataset <- results_df()[[1]]
      dataset[is.na(dataset)] <- 0
      dataset[dataset == ""] <- 0
      
      datatable(
        dataset, 
        options = list(
          pageLength = 5,
          autoWidth = TRUE,
          scrollX = TRUE,
          dom = 'l<"toolbar">rtip'  # Removed 'f' from dom
        ),
        filter = 'top',
        rownames = FALSE,
        extensions = 'Scroller'
      )
    })
    
    output$enrichment_1_pseudo <- renderUI({
      req(results_df())
      fluidRow(column(width = 6,
                      selectInput("pseudo_enrichment_type", label = NULL, choices = c("GO", "BP", "KEGG", "Reactome", "WikiPathways"))),
               column(width = 6,
                      actionButton("run_pseudo_enrichment", class = "btn-primary", "Calculate", width = '100%')
                      ))
      
      
    })
    
    
    pseudo_enrichment_calc <- eventReactive(input$run_pseudo_enrichment,{
      # Check if fenr package is available
      if (!requireNamespace("fenr", quietly = TRUE)) {
        showNotification("❌ Package 'fenr' not available. Enrichment analysis is disabled.", 
                        type = "error", duration = 5)
        return(NULL)
      }
      
      req(results_df(), input$pca_pseudo_bulk_results_table_rows_all, input$pseudo_enrichment_type)
      
      if(input$selection_organism %in% c("Human", "Integrated")){
        load("enrichment_sets/human.RData")
        all_genes <- rownames(results_df()[[1]])
        
        selected_genes <- rownames(results_df()[[1]][input$pca_pseudo_bulk_results_table_rows_all, ])
        
        
        
        bp_human <- prepare_for_enrichment(bp_human$terms, bp_human$mapping, all_genes, feature_name = "gene_symbol")
        go_human <- prepare_for_enrichment(go_human$terms, go_human$mapping, all_genes, feature_name = "gene_symbol")
        kegg_human <- prepare_for_enrichment(kegg_human$terms, kegg_human$mapping, all_genes, feature_name = "gene_symbol")
        reactome_human <- prepare_for_enrichment(reactome_human$terms, reactome_human$mapping, all_genes, feature_name = "gene_symbol")
        wiki_pathways_human <- prepare_for_enrichment(wiki_pathways_human$terms, wiki_pathways_human$mapping, all_genes, feature_name = "text_label")
        
        bp_human <- functional_enrichment(all_genes, selected_genes, bp_human)
        go_human <- functional_enrichment(all_genes, selected_genes, go_human)
        kegg_human <- functional_enrichment(all_genes, selected_genes, kegg_human)
        reactome_human <- functional_enrichment(all_genes, selected_genes, reactome_human)
        wiki_pathways_human <- functional_enrichment(all_genes, selected_genes, wiki_pathways_human)
        
        de_enrichment_calc <- list(bp_human, go_human, kegg_human, reactome_human, wiki_pathways_human)
        
        return(de_enrichment_calc)
        
      }else if(input$selection_organism == "Mouse"){
        
        load("enrichment_sets/mouse.RData")
        all_genes <- rownames(results_df()[[1]])
        
        selected_genes <- rownames(results_df()[[1]][input$pca_pseudo_bulk_results_table_rows_all, ])
        
        bp_mouse <- prepare_for_enrichment(bp_mouse$terms, bp_mouse$mapping, all_genes, feature_name = "gene_symbol")
        go_mouse <- prepare_for_enrichment(go_mouse$terms, go_mouse$mapping, all_genes, feature_name = "gene_symbol")
        kegg_mouse <- prepare_for_enrichment(kegg_mouse$terms, kegg_mouse$mapping, all_genes, feature_name = "gene_symbol")
        reactome_mouse <- prepare_for_enrichment(reactome_mouse$terms, reactome_mouse$mapping, all_genes, feature_name = "gene_symbol")
        wiki_pathways_mouse <- prepare_for_enrichment(wiki_pathways_mouse$terms, wiki_pathways_mouse$mapping, all_genes, feature_name = "text_label")
        
        bp_mouse <- functional_enrichment(all_genes, selected_genes, bp_mouse)
        go_mouse <- functional_enrichment(all_genes, selected_genes, go_mouse)
        kegg_mouse <- functional_enrichment(all_genes, selected_genes, kegg_mouse)
        reactome_mouse <- functional_enrichment(all_genes, selected_genes, reactome_mouse)
        wiki_pathways_mouse <- functional_enrichment(all_genes, selected_genes, wiki_pathways_mouse)
        de_enrichment_calc <- list(bp_mouse, go_mouse, kegg_mouse, reactome_mouse, wiki_pathways_mouse)
        return(de_enrichment_calc)
      }else if(input$selection_organism == "Zebrafish"){
        
        load("enrichment_sets/zebrafish.RData")
        all_genes <- rownames(results_df()[[1]])
        
        selected_genes <- rownames(results_df()[[1]][input$pca_pseudo_bulk_results_table_rows_all, ])
        
        bp_zebrafish <- prepare_for_enrichment(bp_zebrafish$terms, bp_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        go_zebrafish <- prepare_for_enrichment(go_zebrafish$terms, go_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        kegg_zebrafish <- prepare_for_enrichment(kegg_zebrafish$terms, kegg_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        reactome_zebrafish <- prepare_for_enrichment(reactome_zebrafish$terms, reactome_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        wiki_pathways_zebrafish <- prepare_for_enrichment(wiki_pathways_zebrafish$terms, wiki_pathways_zebrafish$mapping, all_genes, feature_name = "text_label")
        
        bp_zebrafish <- functional_enrichment(all_genes, selected_genes, bp_zebrafish)
        go_zebrafish <- functional_enrichment(all_genes, selected_genes, go_zebrafish)
        kegg_zebrafish <- functional_enrichment(all_genes, selected_genes, kegg_zebrafish)
        reactome_zebrafish <- functional_enrichment(all_genes, selected_genes, reactome_zebrafish)
        wiki_pathways_zebrafish <- functional_enrichment(all_genes, selected_genes, wiki_pathways_zebrafish)
        de_enrichment_calc <- list(bp_zebrafish, go_zebrafish, kegg_zebrafish, reactome_zebrafish, wiki_pathways_zebrafish)
        return(de_enrichment_calc)
        
      }
      # else if(input$selection_organism == "Integrated"){
      #   load("enrichment_sets/human.RData")
      #   all_genes <- rownames(results_df()[[1]])
      #   
      #   selected_genes <- rownames(results_df()[[1]][input$pca_pseudo_bulk_results_table_rows_all, ])
      #   
      #   bp_human <- prepare_for_enrichment(bp_human$terms, bp_human$mapping, all_genes, feature_name = "gene_symbol")
      #   go_human <- prepare_for_enrichment(go_human$terms, go_human$mapping, all_genes, feature_name = "gene_symbol")
      #   kegg_human <- prepare_for_enrichment(kegg_human$terms, kegg_human$mapping, all_genes, feature_name = "gene_symbol")
      #   reactome_human <- prepare_for_enrichment(reactome_human$terms, reactome_human$mapping, all_genes, feature_name = "gene_symbol")
      #   wiki_pathways_human <- prepare_for_enrichment(wiki_pathways_human$terms, wiki_pathways_human$mapping, all_genes, feature_name = "gene_symbol")
      #   
      #   bp_human <- functional_enrichment(all_genes, selected_genes, bp_human)
      #   go_human <- functional_enrichment(all_genes, selected_genes, go_human)
      #   kegg_human <- functional_enrichment(all_genes, selected_genes, kegg_human)
      #   reactome_human <- functional_enrichment(all_genes, selected_genes, reactome_human)
      #   wiki_pathways_human <- functional_enrichment(all_genes, selected_genes, wiki_pathways_human)
      #   de_enrichment_calc <- list(bp_human, go_human, kegg_human, reactome_human, wiki_pathways_human)
      #   return(de_enrichment_calc)}
    })

    output$pseudo_enrichment_table <- renderDT({
      if (is.null(pseudo_enrichment_calc())) {
        # Return empty datatable with message
        return(datatable(
          data.frame(Message = "Enrichment analysis not available. Please install the 'fenr' package."),
          options = list(pageLength = 5, dom = 't'),
          rownames = FALSE
        ))
      }
      
      req(pseudo_enrichment_calc())
      if(input$pseudo_enrichment_type == "GO"){
        datatable(
          pseudo_enrichment_calc()[[2]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollY = 300,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip',
            columnDefs = list(
              list(
                targets = which(names(pseudo_enrichment_calc()[[2]]) == "ids") - 1,  # JS indexing starts at 0
                render = JS(
                  "function(data, type, row, meta) {",
                  "  if (type === 'display') {",
                  "    var displayText = data.length > 50 ? data.substr(0, 50) + '???' : data;",
                  "    return '<span title=\"' + data + '\">' + displayText + '</span>';",
                  "  }",
                  "  return data;",
                  "}"
                )
              )
            )
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(pseudo_enrichment_calc()[[2]], is.numeric)), digits = 3)
      }else if(input$pseudo_enrichment_type == "BP"){
        datatable(
          pseudo_enrichment_calc()[[1]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollY = 300,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip',
            columnDefs = list(
              list(
                targets = which(names(pseudo_enrichment_calc()[[2]]) == "ids") - 1,  # JS indexing starts at 0
                render = JS(
                  "function(data, type, row, meta) {",
                  "  if (type === 'display') {",
                  "    var displayText = data.length > 50 ? data.substr(0, 50) + '???' : data;",
                  "    return '<span title=\"' + data + '\">' + displayText + '</span>';",
                  "  }",
                  "  return data;",
                  "}"
                )
              )
            )
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(pseudo_enrichment_calc()[[1]], is.numeric)), digits = 3)
      }else if(input$pseudo_enrichment_type == "KEGG"){
        datatable(
          pseudo_enrichment_calc()[[3]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollY = 300,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip',
            columnDefs = list(
              list(
                targets = which(names(pseudo_enrichment_calc()[[2]]) == "ids") - 1,  # JS indexing starts at 0
                render = JS(
                  "function(data, type, row, meta) {",
                  "  if (type === 'display') {",
                  "    var displayText = data.length > 50 ? data.substr(0, 50) + '???' : data;",
                  "    return '<span title=\"' + data + '\">' + displayText + '</span>';",
                  "  }",
                  "  return data;",
                  "}"
                )
              )
            )
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(pseudo_enrichment_calc()[[3]], is.numeric)), digits = 3)
      }else if(input$pseudo_enrichment_type == "Reactome"){
        datatable(
          pseudo_enrichment_calc()[[4]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollY = 300,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip',
            columnDefs = list(
              list(
                targets = which(names(pseudo_enrichment_calc()[[2]]) == "ids") - 1,  # JS indexing starts at 0
                render = JS(
                  "function(data, type, row, meta) {",
                  "  if (type === 'display') {",
                  "    var displayText = data.length > 50 ? data.substr(0, 50) + '???' : data;",
                  "    return '<span title=\"' + data + '\">' + displayText + '</span>';",
                  "  }",
                  "  return data;",
                  "}"
                )
              )
            )
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(pseudo_enrichment_calc()[[4]], is.numeric)), digits = 3)
      }else if(input$pseudo_enrichment_type == "WikiPathways"){
        datatable(
          pseudo_enrichment_calc()[[5]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollY = 300,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip',
            columnDefs = list(
              list(
                targets = which(names(pseudo_enrichment_calc()[[2]]) == "ids") - 1,  # JS indexing starts at 0
                render = JS(
                  "function(data, type, row, meta) {",
                  "  if (type === 'display') {",
                  "    var displayText = data.length > 50 ? data.substr(0, 50) + '???' : data;",
                  "    return '<span title=\"' + data + '\">' + displayText + '</span>';",
                  "  }",
                  "  return data;",
                  "}"
                )
              )
            )
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(pseudo_enrichment_calc()[[5]], is.numeric)), digits = 3)
      }
    })
    
    output$pseudo_ora_selection_ui <- renderUI({
      req(results_df())
      fluidRow(column(width = 6,
                      selectInput("pseudo_ora_selection", label = NULL, choices = c("CollecTRI", "PROGENy", "MSigDB Hallmark"))),
               column(width = 6,
                      actionButton("pseudo_run_ora", class = "btn-primary", "Calculate", width = '100%')))

    })
    
    collectri <- eventReactive(input$pseudo_run_ora,{
      req(input$pseudo_ora_selection == "CollecTRI")
      collectri <- read_rds("enrichment_sets/collectri.rds")
      tf_acts <- dc$run_ulm(mat = results_df()[[2]], net = collectri)
      tf_pvals <- tf_acts[[2]]
      tf_acts <- tf_acts[[1]]
      rownames(tf_pvals) <- "Selected_Clusters"
      rownames(tf_acts) <- "Selected_Clusters"
      results <- list(tf_acts, tf_pvals, collectri)
      return(results)
    })
    
    output$collectri_outputs <- renderUI({
      req(collectri())
      fluidRow(column(width = 6,
                      selectInput("collectri_plot_selection", label = NULL, choices = c("Top 25", "Volcano", "Network"))),
               column(width = 6,
                      conditionalPanel(condition = "input.collectri_plot_selection == 'Volcano'",
                                       selectizeInput("collectri_volcano_name", label = NULL, choices = colnames(collectri()[[1]]))),
                      conditionalPanel(condition = "input.collectri_plot_selection == 'Network'",
                                       selectizeInput("collectri_network_names", label = NULL, choices = colnames(collectri()[[1]]), multiple =T)
                                       )))
    })
    
    output$progeny_outputs <- renderUI({
      req(progeny())
      fluidRow(column(width = 6,
                      selectInput("progeny_plot_selection", label =NULL, choices = c("Scores", "Targets"))),
               column(width = 6,
                      conditionalPanel(condition = "input.progeny_plot_selection == 'Targets'",
                                       selectizeInput("progeny_scatter_selection", label = NULL, choices = c("Androgen", "EGFR", "Estrogen", "Hypoxia", "JAK-STAT", "MAPK",
                                                                                                           "NFkB", "PI3K", "TGFb", "TNFa", "Trail", "VEGF", "WNT", "p53"))
                                       )))
    })
    
    output$msigdb_outputs <- renderUI({
      req(msigdb())
      fluidRow(column(width = 6,
                      selectInput("msigdb_plot_selection", label =NULL, choices = c("Top Signatures", "Running Score"))),
               column(width = 6,
                      conditionalPanel(condition = "input.msigdb_plot_selection == 'Running Score'",
                                       selectizeInput("msigdb_running_selection", label = NULL, choices = sort(unique(msigdb()[[2]]$geneset)))
                      )))
    })
    
    
    output$collectri_top <- renderImage({
      req(collectri())
      dc$plot_barplot(
        acts=collectri()[[1]],
        contrast='Selected_Clusters',
        top= as.integer(25),
        vertical=F,
        # figsize=c(6, 10),
        save = "collectri.png"
      )
      list(src = "collectri.png")
    }, deleteFile = TRUE)
    
    output$collectri_volcano <- renderImage({
      req(collectri())
      logFCs = as.data.frame(t(results_df()[[1]]['log2FoldChange']))
      pvals = as.data.frame(t(results_df()[[1]]['padj']))
      rownames(logFCs) <- "Selected_Clusters"
      rownames(pvals) <- "Selected_Clusters"
      
      dc$plot_volcano(
        logFCs=logFCs,
        pvals=pvals,
        contrast='Selected_Clusters',
        name=input$collectri_volcano_name,
        net=collectri()[[3]],
        top=as.integer(10),
        sign_thr=0.05,
        lFCs_thr=0.5,
        save = "collectri_volcano.png"
      )
      list(src = "collectri_volcano.png")
    }, deleteFile = TRUE)
    
    
    output$collectri_network <- renderImage({
      req(collectri(),!is.null(input$collectri_network_names))
      dc$plot_network(
        net=collectri()[[3]],
        obs=results_df()[[2]],
        act=collectri()[[1]],
        n_sources= input$collectri_network_names,
        n_targets= as.integer(15),
        node_size=as.integer(100),
        figsize = c(5,5),
        c_pos_w='darkgreen',
        c_neg_w='darkred',
        vcenter=T,
        save = "collectri_network.png"
      )

      list(src = "collectri_network.png")
    }, deleteFile = TRUE)
    
    progeny <- eventReactive(input$pseudo_run_ora,{
      req(input$pseudo_ora_selection == "PROGENy")
      progeny <- read_rds("enrichment_sets/progeny.rds")
      mat = results_df()[[2]]
      
      pathway_acts = dc$run_mlm(mat=mat, net=progeny)
      pathway_acts <- pathway_acts[[1]]
      rownames(pathway_acts) <- "Selected_Clusters"
      results <- list(pathway_acts, progeny)
      return(results)
    })
    
    output$progeny_top <- renderImage({
      req(progeny())
      dc$plot_barplot(
        acts=progeny()[[1]],
        contrast='Selected_Clusters',
        top=as.integer(25),
        vertical=F,
        save = "progeny.png"
      )
      
      list(src = "progeny.png")
    }, deleteFile = TRUE)
    
    output$progeny_targets <- renderImage({
      req(progeny())
      dc$plot_targets(
        data=results_df()[[1]],
        stat='stat',
        source_name= input$progeny_scatter_selection,
        net=progeny()[[2]],
        top=as.integer(15),
        save = "progeny_targets.png"
      )
      
      list(src = "progeny_targets.png")
    }, deleteFile = TRUE)
    
    msigdb <- eventReactive(input$pseudo_run_ora,{
      req(input$pseudo_ora_selection == "MSigDB Hallmark")
      
      top_genes = results_df()[[1]][results_df()[[1]]['padj'] < 0.05,]
      msigdb <- read_rds("enrichment_sets/msigdb.rds")
      enr_pvals = dc$get_ora_df(
        df=top_genes,
        net=msigdb,
        source='geneset',
        target='genesymbol'
      )
      
      enr_pvals <- enr_pvals[order(enr_pvals$`Combined score`,decreasing = T),]
      
      results <- list(enr_pvals, msigdb)
      saveRDS(enr_pvals,"asdasds.rds")
      return(results)
    })
    
    create_blank_image <- function(path = "blank.png") {
      png(path, width = 50, height = 50)
      par(mar = c(0, 0, 0, 0))  # remove all margins
      plot.new()
      dev.off()
      return(path)
    }
    
    output$msigdb_top <- renderImage({
      req(msigdb())
      data <- head(msigdb()[[1]], 15)
      
      # Fix type and clean NAs
      data <- data[!is.na(data[["Odds ratio"]]), ]
      
      if (nrow(data) == 0) {
        showNotification("No valid Odds ratio values to plot.", type = "error")
        return(list(src = create_blank_image()))
      }else{
        dc$plot_dotplot(
          df = data,
          x = 'Combined score',
          y = 'Term',
          s = 'Odds ratio',
          c = 'FDR p-value',
          scale = 0.1,
          save = "msigdb_top.png"
        )
        
        list(src = "msigdb_top.png")
      }
    }, deleteFile = TRUE)
    
    
    output$msigdb_running <- renderImage({
      req(msigdb(), input$msigdb_running_selection)
      
      dc$plot_running_score(
        df=results_df()[[1]],
        stat='stat',
        net=msigdb()[[2]],
        source='geneset',
        target='genesymbol',
        set_name= input$msigdb_running_selection,
        save = "msigdb_running.png"
      )
      
      list(src = "msigdb_running.png")
    }, deleteFile = TRUE)
    
    
    # ========================================
    # CSV EXPORT HANDLERS
    # Added: 13 octobre 2025
    # Feature: Export analysis results to CSV
    # ========================================
    
    # 1. Export Cell Type Markers
    output$download_markers <- downloadHandler(
      filename = function() {
        paste0("Cell_markers_", Sys.Date(), ".csv")
      },
      content = function(file) {
        tryCatch({
          req(adata(), input$selection_rank_select)
          
          withProgress(message = "Exporting cell markers...", value = 0.5, {
            data_markers <- sc$get$rank_genes_groups_df(adata(), group = input$selection_rank_select)
            
            if (is.null(data_markers) || nrow(data_markers) == 0) {
              stop("No marker data available to export")
            }
            
            write.csv(
              data_markers,
              file,
              row.names = FALSE,
              quote = TRUE
            )
            
            showNotification(
              paste("Exported", nrow(data_markers), "marker genes"),
              type = "message",
              duration = 3
            )
          })
        }, error = function(e) {
          showNotification(
            paste("Error exporting markers:", e$message),
            type = "error",
            duration = 5
          )
          # Write empty CSV to avoid download error
          write.csv(data.frame(Error = e$message), file, row.names = FALSE)
        })
      }
    )
    
    # 2. Export Gene Correlation Results
    output$download_correlation <- downloadHandler(
      filename = function() {
        paste0("Correlation_results_", Sys.Date(), ".csv")
      },
      content = function(file) {
        tryCatch({
          req(statistics_coexpression())
          
          withProgress(message = "Exporting correlation results...", value = 0.5, {
            # statistics_coexpression() returns a list: [[1]] = data, [[2]] = first_gene, [[3]] = second_gene
            corr_data <- statistics_coexpression()[[1]]
            
            if (is.null(corr_data) || nrow(corr_data) == 0) {
              stop("No correlation data available to export")
            }
            
            write.csv(
              as.data.frame(corr_data),
              file,
              row.names = TRUE,
              quote = TRUE
            )
            
            showNotification(
              paste("Exported", nrow(corr_data), "correlation data points"),
              type = "message",
              duration = 3
            )
          })
        }, error = function(e) {
          showNotification(
            paste("Error exporting correlation results:", e$message),
            type = "error",
            duration = 5
          )
          write.csv(data.frame(Error = e$message), file, row.names = FALSE)
        })
      }
    )
    
    # 3. Export Differential Gene Expression (DGE)
    output$download_dge <- downloadHandler(
      filename = function() {
        paste0("DGE_results_", Sys.Date(), ".csv")
      },
      content = function(file) {
        tryCatch({
          req(de_dge_calculation(), input$de_ident_1_name)
          
          withProgress(message = "Exporting DGE results...", value = 0.5, {
            # Extract DGE results the same way as in renderDT
            group_name <- list(input$de_ident_1_name)
            dge_data <- sc$get$rank_genes_groups_df(de_dge_calculation(), group = group_name, key = 'rank_genes_groups')
            
            # Rename columns to match the display
            colnames(dge_data) <- c("Gene", "Scores", "LogFC", "p-val", "adj-p", "pct")
            
            if (is.null(dge_data) || nrow(dge_data) == 0) {
              stop("No DGE data available to export")
            }
            
            write.csv(
              dge_data,
              file,
              row.names = FALSE,
              quote = TRUE
            )
            
            showNotification(
              paste("Exported", nrow(dge_data), "genes from DGE analysis"),
              type = "message",
              duration = 3
            )
          })
        }, error = function(e) {
          showNotification(
            paste("Error exporting DGE results:", e$message),
            type = "error",
            duration = 5
          )
          write.csv(data.frame(Error = e$message), file, row.names = FALSE)
        })
      }
    )
    
    # 4. Export Enrichment Analysis (DE)
    output$download_enrichment <- downloadHandler(
      filename = function() {
        enrichment_type <- input$de_enrichment_type
        if (is.null(enrichment_type)) enrichment_type <- "enrichment"
        paste0("Enrichment_", enrichment_type, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        tryCatch({
          req(de_enrichment_calc())
          
          withProgress(message = "Exporting enrichment results...", value = 0.5, {
            # Get the appropriate enrichment data based on selection
            enrichment_data <- if (!is.null(input$de_enrichment_type)) {
              switch(input$de_enrichment_type,
                     "GO" = de_enrichment_calc()[[2]],
                     "BP" = de_enrichment_calc()[[1]],
                     "KEGG" = de_enrichment_calc()[[3]],
                     "Reactome" = de_enrichment_calc()[[4]],
                     "WikiPathways" = de_enrichment_calc()[[5]],
                     de_enrichment_calc()[[1]])  # default
            } else {
              de_enrichment_calc()[[1]]
            }
            
            if (is.null(enrichment_data) || nrow(enrichment_data) == 0) {
              stop("No enrichment data available to export")
            }
            
            write.csv(
              as.data.frame(enrichment_data),
              file,
              row.names = FALSE,
              quote = TRUE
            )
            
            showNotification(
              "Enrichment results exported successfully",
              type = "message",
              duration = 3
            )
          })
        }, error = function(e) {
          showNotification(
            paste("Error exporting enrichment results:", e$message),
            type = "error",
            duration = 5
          )
          write.csv(data.frame(Error = e$message), file, row.names = FALSE)
        })
      }
    )
    
    # 5. Export Pseudo-bulk Analysis Results
    output$download_pseudobulk <- downloadHandler(
      filename = function() {
        paste0("Pseudobulk_results_", Sys.Date(), ".csv")
      },
      content = function(file) {
        tryCatch({
          req(results_df()[[1]])
          
          withProgress(message = "Exporting pseudo-bulk results...", value = 0.5, {
            dataset <- results_df()[[1]]
            
            if (is.null(dataset) || nrow(dataset) == 0) {
              stop("No pseudo-bulk data available to export")
            }
            
            dataset[is.na(dataset)] <- 0
            dataset[dataset == ""] <- 0
            
            write.csv(
              dataset,
              file,
              row.names = FALSE,
              quote = TRUE
            )
            
            showNotification(
              paste("Exported", nrow(dataset), "genes from pseudo-bulk analysis"),
              type = "message",
              duration = 3
            )
          })
        }, error = function(e) {
          showNotification(
            paste("Error exporting pseudo-bulk results:", e$message),
            type = "error",
            duration = 5
          )
          write.csv(data.frame(Error = e$message), file, row.names = FALSE)
        })
      }
    )
    
    # 6. Export Pseudo-bulk Enrichment
    output$download_pseudo_enrichment <- downloadHandler(
      filename = function() {
        enrichment_type <- input$pseudo_enrichment_type
        if (is.null(enrichment_type)) enrichment_type <- "enrichment"
        paste0("Pseudobulk_enrichment_", enrichment_type, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        tryCatch({
          req(pseudo_enrichment_calc())
          
          withProgress(message = "Exporting pseudo-bulk enrichment...", value = 0.5, {
            # Get the appropriate enrichment data based on selection
            enrichment_data <- if (!is.null(input$pseudo_enrichment_type)) {
              switch(input$pseudo_enrichment_type,
                     "GO" = pseudo_enrichment_calc()[[2]],
                     "BP" = pseudo_enrichment_calc()[[1]],
                     "KEGG" = pseudo_enrichment_calc()[[3]],
                     "Reactome" = pseudo_enrichment_calc()[[4]],
                     "WikiPathways" = pseudo_enrichment_calc()[[5]],
                     pseudo_enrichment_calc()[[1]])  # default
            } else {
              pseudo_enrichment_calc()[[1]]
            }
            
            if (is.null(enrichment_data) || nrow(enrichment_data) == 0) {
              stop("No pseudo-bulk enrichment data available to export")
            }
            
            write.csv(
              as.data.frame(enrichment_data),
              file,
              row.names = FALSE,
              quote = TRUE
            )
            
            showNotification(
              "Pseudo-bulk enrichment exported successfully",
              type = "message",
              duration = 3
            )
          })
        }, error = function(e) {
          showNotification(
            paste("Error exporting pseudo-bulk enrichment:", e$message),
            type = "error",
            duration = 5
          )
          write.csv(data.frame(Error = e$message), file, row.names = FALSE)
        })
      }
    )
    
    
}

options(shiny.sanitize.errors = TRUE)

shinyApp(ui = ui, server = server)

