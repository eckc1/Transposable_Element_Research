library(shiny)
library(ggplot2)

# Define the paths to the CSV files
file_path <- "SingleB73Extract.csv"
#file_path <- "DoubleB73Extract.csv"
#file_path <- "A188_extract.csv"
#file_path <- "W22_Extract.csv"
new_file_path <- "KCT_data.csv"  # Adjust the file name as necessary

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
                breaks = c(0, 1, 2, 5, 20, 120, Inf),
                labels = c('1', '2', '3-5', '6-20', '21-120', '>120'),
                include.lowest = TRUE)

# UI definition
ui <- fluidPage(
  titlePanel("Chromosome Map of Query Locations"),
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("selected_bins",
                         "Select bins to display:",
                         choices = levels(data$Bin),
                         selected = levels(data$Bin)),
      checkboxInput("show_centromeres", "Centromeres", value = TRUE),
      checkboxInput("show_knobs", "Knobs", value = TRUE),
      checkboxInput("show_telomeres", "Telomeres", value = TRUE)
    ),
    mainPanel(
      plotOutput("chromosomePlot")
    )
  )
)

# Server logic
server <- function(input, output) {
  output$chromosomePlot <- renderPlot({
    filtered_data <- data[data$Bin %in% input$selected_bins, ]
    
    # Start with the main plot for filtered data
    p <- ggplot(filtered_data) +
      geom_segment(aes(x = sStart, xend = sEnd, y = Chromosome, yend = Chromosome),
                   color = "grey", size = 1) +
      geom_point(aes(x = sStart, y = Chromosome, color = Bin)) +
      labs(x = "Position", y = "Chromosome", color = "NumHits Bin") +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(angle = 0, hjust = 0))
    
    # Add centromeres, knobs, and telomeres based on checkbox inputs
    if (input$show_centromeres) {
      centromere_data <- new_data[new_data$Classification == 'Centromere', ]
      p <- p + geom_point(data = centromere_data, aes(x = sStart, y = sequence_region, color = "Centromere"), size = 3)
    }
    
    if (input$show_knobs) {
      knob_data <- new_data[new_data$Classification == 'Knob', ]  # Adjust as needed
      p <- p + geom_point(data = knob_data, aes(x = sStart, y = sequence_region, color = "Knob"), size = 3)
    }
    
    if (input$show_telomeres) {
      telomere_data <- new_data[new_data$Classification == 'Telomere', ]  # Adjust as needed
      p <- p + geom_point(data = telomere_data, aes(x = sStart, y = sequence_region, color = "Telomere"), size = 3)
    }
    
    print(p)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
