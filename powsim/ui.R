#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("GWAS Power Calculator"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            fileInput("input_file",
                      "Choose a gwasim data (recommended) or a VCF file:",
                      multiple = F,
                      accept = c(".vcf", ".rdat", ".vcf.gz"),
                      buttonLabel = "Browse...",
                      placeholder = "chr1.rdat"),
            sliderInput("p_val",
                        "Set p-value threshold:",
                        min = 0.01,
                        max = 0.1,
                        value = 0.05),
            sliderInput("beta",
                        "Set the effect size (Beta):",
                        min = 0.0,
                        max = 5.0,
                        value = 1.0,
                        step = 0.1)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("power_plot"),
            plotOutput("max_q2_plot"),
            tableOutput("input_info"),
            plotOutput("maf_plot")
        )
    )
))
