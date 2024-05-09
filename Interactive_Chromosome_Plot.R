library(shiny)
library(ggplot2)

# Define the path to the CSV file
file_path <- "SingleB73Extract.csv"

# Check if the file exists before proceeding
if (!file.exists(file_path)) {
  stop("File not found: ", file_path)
}

# Load the data outside of the server function so it's only loaded once
data <- read.csv(file_path)
data <- subset(data, !grepl("^scaf", Chromosome) & Chromosome != "FALSE")
data$NumHitsForGenome <- as.numeric(data$NumHitsForGenome)

# Pre-process the data to create bins for NumHitsForGenome
data$Bin <- cut(data$NumHitsForGenome,
                breaks = c(1, 2, 5, 20, 120, Inf),
                labels = c('2', '3-5', '6-20', '21-120', '>120'),
                include.lowest = TRUE)

# UI definition
ui <- fluidPage(
  titlePanel("Chromosome Map of Query Locations"),
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("selected_bins",
                         "Select bins to display:",
                         choices = levels(data$Bin),
                         selected = levels(data$Bin))
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
    
    # Calculate min and max sStart for each chromosome in the filtered data
    chromosome_min = aggregate(sStart ~ Chromosome, filtered_data, min)
    chromosome_max = aggregate(sStart ~ Chromosome, filtered_data, max)
    
    # Merge min and max into a single data frame
    chromosome_ranges <- merge(chromosome_min, chromosome_max, by = "Chromosome")
    names(chromosome_ranges) <- c("Chromosome", "sStart_min", "sStart_max")
    
    ggplot(filtered_data) +
      geom_segment(data = chromosome_ranges, aes(x = sStart_min, xend = sStart_max, y = Chromosome, yend = Chromosome),
                   color = "grey", size = 1) +
      geom_point(aes(x = sStart, y = Chromosome, color = Bin)) +
      theme_minimal() +
      labs(x = "Position", y = "Chromosome", color = "NumHits Bin") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(angle = 0, hjust = 0))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

