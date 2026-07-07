library(shiny)
library(ggplot2)

# -----------------------------
# Helper function
# -----------------------------
parse_snps <- function(text_input){
  
  snps <- unlist(strsplit(text_input, "\n"))
  snps <- trimws(snps)
  snps <- snps[snps != ""]
  
  df <- data.frame(
    raw = snps,
    stringsAsFactors = FALSE)
  
  # Extract chromosome and position
  split_coords <- strsplit(df$raw, ":")
  
  df$chr <- sapply(split_coords, `[`, 1)
  df$pos <- as.numeric(sapply(split_coords, `[`, 2))
  
  return(df)
}

# -----------------------------
# UI
# -----------------------------
ui <- fluidPage(
  
  titlePanel("Genetic Coordinate Viewer"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      textAreaInput(
        inputId = "snp_text",
        label = "Enter SNP coordinates",
        placeholder = "chr1:12345\nchr2:56789",
        rows = 10
      ),
      
      actionButton(
        inputId = "run_analysis",
        label = "Run"
      ),
      
      br(),
      br(),
      
      downloadButton(
        outputId = "download_txt",
        label = "Download TXT"
      ),
      
      downloadButton(
        outputId = "download_png",
        label = "Download PNG"
      )
    ),
    
    mainPanel(
      
      h4("Parsed SNPs"),
      
      tableOutput("snp_table"),
      
      br(),
      
      h4("Position Plot"),
      
      plotOutput("snp_plot", height = "500px")
    )
  )
)

# -----------------------------
# Server
# -----------------------------
server <- function(input, output, session){
  
  snp_data <- eventReactive(input$run_analysis, {
    
    req(input$snp_text)
    
    parse_snps(input$snp_text)
  })
  
  # -----------------------------
  # Table output
  # -----------------------------
  output$snp_table <- renderTable({
    
    snp_data()
  })
  
  # -----------------------------
  # Plot output
  # -----------------------------
  output$snp_plot <- renderPlot({
    
    df <- snp_data()
    
    ggplot(df, aes(x = pos, y = chr)) +
      geom_point(size = 3) +
      theme_bw() +
      labs(
        title = "Genetic Coordinates",
        x = "Genomic Position",
        y = "Chromosome"
      )
  })
  
  # -----------------------------
  # Download text file
  # -----------------------------
  output$download_txt <- downloadHandler(
    
    filename = function() {
      paste0("snp_results_", Sys.Date(), ".txt")
    },
    
    content = function(file) {
      
      write.table(
        snp_data(),
        file = file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
  )
  
  # -----------------------------
  # Download PNG
  # -----------------------------
  output$download_png <- downloadHandler(
    
    filename = function() {
      paste0("snp_plot_", Sys.Date(), ".png")
    },
    
    content = function(file) {
      
      df <- snp_data()
      
      p <- ggplot(df, aes(x = pos, y = chr)) +
        geom_point(size = 3) +
        theme_bw() +
        labs(
          title = "Genetic Coordinates",
          x = "Genomic Position",
          y = "Chromosome"
        )
      
      ggsave(
        filename = file,
        plot = p,
        width = 8,
        height = 5,
        dpi = 300
      )
    }
  )
}

# -----------------------------
# Run app
# -----------------------------
shinyApp(ui = ui, server = server)