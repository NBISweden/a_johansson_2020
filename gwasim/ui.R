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
    titlePanel("gwasim example"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("n_mrk",
                        "Max number of loci:",
                        min = 1,
                        max = 10,
                        value = 1),
            sliderInput("rare_thr",
                        "Rare allele maf threshold:",
                        min = 0,
                        max = 0.1,
                        value = 0.01),
            radioButtons("use_rare",
                          "Simulate effect of ",
                         choiceNames = c('rare alleles only', 'common alleles only', 'both common and rare alleles'),
                         choiceValues = c('rare', 'common', 'both')),
            sliderInput("perc_neg",
                        "Percentage of loci with negative effect:",
                        min = 0,
                        max = 100,
                        value = 0),
            textInput("e_mean",
                      "Error mean:",
                      value = "0",
                      width = "40%",
                      placeholder = "Enter mean for the error term, e.g. 0"),
            textInput("e_sd",
                      "Error std. dev.:",
                      value = "1",
                      width = "40%",
                      placeholder = "Enter std. dev. for the error term, e.g. 1"),
            submitButton("Run simulation")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("distPlot")
        )
    )
))
