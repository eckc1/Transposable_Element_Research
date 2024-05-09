library(shiny)
library(ggplot2)
library(plotly)
library(htmlwidgets)

# Define the paths to the CSV files
file_path <- "SingleB73Extract.csv"
new_file_path <- "KCT_data.csv"

# Check if the files exist before proceeding
if (!file.exists(file_path)) {
  stop("File not found: ", file_path)
}
if (!file.exists(new_file_path)) {
  stop("File not found: ", new_file_path)
}

# Load the data
data <- read.csv(file_path)
new_data <- read.csv(new_file_path)

# Pre-process the data
data <- subset(data, !grepl("^scaf", Chromosome) & Chromosome != "FALSE")
data$NumHitsForGenome <- as.numeric(data$NumHitsForGenome)
data$Bin <- cut(data$NumHitsForGenome,
                breaks = c(0,1, 2, 5, 20, 120, Inf),
                labels = c('1','2', '3-5', '6-20', '21-120', '>120'))

# UI definition
ui <- fluidPage(
  tags$head(tags$script(HTML("
    Shiny.addCustomMessageHandler('openPlot', function(message) {
      window.open(message.url, '_blank');
    });
  "))),
  titlePanel("Chromosome Map of Query Locations"),
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("selected_bins",
                         "Select bins to display:",
                         choices = levels(data$Bin),
                         selected = levels(data$Bin)),
      #checkboxInput("show_centromeres", "Centromeres", value = TRUE),
      #checkboxInput("show_knobs", "Knobs", value = TRUE),
      #checkboxInput("show_telomeres", "Telomeres", value = TRUE)
    ),
    mainPanel(
      plotlyOutput("chromosomePlot")
    )
  )
  
)

# Server logic
server <- function(input, output, session) {
  output$chromosomePlot <- renderPlotly({
    filtered_data <- data[data$Bin %in% input$selected_bins, ]
    p <- ggplot(filtered_data, aes(x = sStart, fill = Chromosome)) +
      geom_histogram(binwidth = 1000000, position = "identity", alpha = 0.5) +
      facet_wrap(~ Chromosome) +
      labs(x = "Position", y = "Frequency") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Convert ggplot object to plotly
    ggplotly(p)
  })
  
  observeEvent(input$chromosomePlot_click, {
    click_data <- input$chromosomePlot_click
    plot <- plotly_build(output$chromosomePlot())
    tmp_file <- tempfile(fileext = ".html")
    saveWidget(plot, tmp_file, selfcontained = TRUE)
    session$sendCustomMessage('openPlot', list(url = tmp_file))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)