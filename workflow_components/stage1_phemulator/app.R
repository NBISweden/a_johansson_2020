#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(kableExtra)

# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("Phemulator"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("N",
                  "Number of individuals:",
                  min=1000,
                  max=10000,
                  step = 1000,
                  value=5000),
      sliderInput("m",
                  "Number of markers:",
                  min=1,
                  max=10,
                  step=1,
                  value=5),
      sliderInput("betas_mu",
                  "Mean effect:",
                  min=0,
                  max=5,
                  step=.5,
                  value=1),
      sliderInput("betas_sd",
                  "Effect SD:",
                  min=0,
                  max=3,
                  step=0.1,
                  value=1),
      sliderInput("e_mu",
                  "Mean error:",
                  min=0,
                  max=0.5,
                  step=0.01,
                  value=0.25),
      sliderInput("e_sd",
                  "Error SD:",
                  min=0,
                  max=1,
                  step=0.01,
                  value=1),
      sliderInput("m_neg",
                  "Number of markers with negative effect:",
                  value = 0,
                  step = 1,
                  min = 0,
                  max = 5
      ),
      sliderInput("q",
                  "Minor allele frequency:",
                  min=0,
                  max=0.5,
                  step=0.01,
                  value=.33),
      sliderInput("N_sim",
                "Number of simulations:",
                min = 1,
                max = 100,
                value=30,
      ),
    ),


    # Show a plot of the generated distribution
    mainPanel(
      fluidRow(
        splitLayout(cellWidths = c("70%", "25%"),
                    plotOutput("phenos"),
                    plotOutput("phenos2"))
      ),
      fluidRow(
        splitLayout(cellWidths = c("47.5%", "47.5%"),
                    plotOutput("effects"),
                    plotOutput("errors"))
      ),
      textOutput("varexp"),
      tableOutput("genos"),
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  rvals <- reactiveValues(
    G = NULL,
    betas = NULL,
    errors = c(0),
    phenos = c(0),
    varexp = 0
  )

  observeEvent(c(input$N, input$N_sim, input$m, input$q), {
    p <- 1 - input$q
    q <- input$q
    gt <- sample(x = c(0,1,2),
                 size = (input$N * input$m),
                 replace = T,
                 prob = c(p^2, p*q, q^2))
    rvals$G <- matrix(gt, ncol = input$m, byrow = T)
  })

  observe({
    betas <- abs(rnorm(n = input$m,
                       mean = input$betas_mu,
                       sd = input$betas_sd))
    negative_idx <- sample(c(1:input$m), size=input$m_neg, replace=F)
    signs <- rep(1, times=input$m)
    betas[negative_idx] <- -betas[negative_idx]
    rvals$betas <- betas
  })

  observe({
    N <- dim(rvals$G)[1]
    m <- dim(rvals$G)[2]
    errors <- rnorm(n = N, mean = input$e_mu, sd = input$e_sd)
    rvals$errors <- errors
  })

  observe({
    phenos <- round(rvals$G %*% rvals$betas + rvals$errors, digits = 3)
    rvals$phenos <- phenos
  })

  observe({
    y <- rvals$phenos
    G <- rvals$G
    M <- data.frame(y, G)
    colnames(M) <- c('y', paste0("g", 1:dim(G)[2]))
    mf <- as.formula(paste0("y~", paste0(colnames(M)[-1], collapse="+")))
    model <- lm(formula = mf, data = M)
    rvals$varexp <- (var(M$y) - var(model$residuals)) / var(M$y)
  })

  output$genos <- function(){
    N_ids <- dim(rvals$G)[1]
    maf <- round(colSums(rvals$G)/(2*N_ids), digits=3)
    x <- rvals$G[c(1:5, N_ids, N_ids),]
    x <- apply(x, MARGIN = c(1,2), as.integer)
    x[6,] <- "..."
    x <- rbind(round(rvals$betas, digits = 2), maf, x)
    colnames(x) <- paste0('SNP', 1:input$m)
    rownames(x) <- c("effect", "maf", paste0("ind_", 1:5), "ind_...", paste0("ind_", N_ids))
    x %>% knitr::kable("html") %>%
      kable_styling("striped", full_width = T) %>%
      row_spec(c(1,2), background  = "lightblue")
  }

  output$errors <- renderPlot({
    x <- data.frame(errors = rvals$errors)
    ggplot(x, mapping = aes(x = errors)) +
      geom_histogram(colour='white', fill='olivedrab') +
      theme_bw()
  })

  output$phenos <- renderPlot({
    x <- data.frame(phenos = rvals$phenos)
    ggplot(x, mapping = aes(x = phenos)) +
      geom_histogram(colour = 'white', fill = 'orange') +
      theme_bw() +
      theme(legend.position = 'none')
  })

  output$phenos2 <- renderPlot({
    x <- data.frame(phenos = rvals$phenos)
    ggplot(x, mapping = aes(y = phenos)) +
      geom_boxplot(fill = 'orange') +
      theme_bw() +
      theme(legend.position = 'none')
  })

  output$effects <- renderPlot({
    x <-  rnorm(n = 1000, mean = input$betas_mu, sd = input$betas_sd)
    x <- data.frame(effects = x)
    ggplot(x, mapping = aes(x = effects)) +
      geom_histogram(colour = 'white', fill='slateblue') +
      theme_bw() +
      theme(legend.position = 'none')
  })

  output$varexp <- renderText({
    varexp <- round(rvals$varexp, digits=2) * 100
    paste0("Variance explained by the y ~ 1 + G model: ", varexp, "%")
  })
}

# Run the application
shinyApp(ui = ui, server = server)
