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
    titlePanel("GWASIM - data-based Genome-Wide Association SIMulations"),
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            fileInput('vcf_file', 'Choose VCF (or rdat) file with variants',
                      accept=c('text/csv',
                               'text/plain',
                               c('.vcf.gz', '.vcf', '.rdat', '.dat'))),
            fileInput('bed_file', 'Choose bed file with region definitions',
                      accept=c('text/csv',
                               'text/plain',
                               '.bed')),
            radioButtons("use_maf",
                         label = "Simulate the effect of:",
                         choiceValues = c('rare', 'common', 'common_rare'),
                         choiceNames = c('rare variants only', 'common variants only', 'common and rare variants')
                         ),
            sliderInput("maf_threshold",
                        "Minor allele frequency threshold:",
                        min = 0.01,
                        max = 0.1,
                        value = 0.05)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tableOutput("vcf_summary"),
            plotOutput("bed_summary")
        )
    )
))
