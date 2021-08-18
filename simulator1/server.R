#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(ggplot2)
library(forcats)

mem_limit <- 1024 # 1GB
nrows_limit <- 500
mem_limit_bytes <- mem_limit * 1024^2
options(shiny.maxRequestSize = mem_limit_bytes) # set max upload size

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  output$vcf_summary <-
    renderTable({
      req(input$vcf_file)
      tryCatch(
        {
          if (grepl("\\.vcf$", input$vcf_file$datapath) | grepl("\\.vcf.gz$", input$vcf_file$datapath)) {
            datapath <- input$vcf_file$datapath
            #datapath <- "~/Projects/johansson/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
            vcf_data <- vcfR::read.vcfR(file = datapath,
                                        nrows = nrows_limit,
                                        limit = mem_limit_bytes)
            shp <- dim(vcf_data@gt)
            print(shp)
            df <- data.frame(c("Filename","N markers", "N individuals"), c(datapath, shp))
            colnames(df) <- c("Feature", "Value")
          }
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )
      return(df)
    }) # renderTable

        output$bed_summary <- renderPlot({
            req(input$bed_file)
            tryCatch({
                regions <- read_delim(file = input$bed_file$datapath,
                                      delim = " ",
                                      col_names = c("chr", "start", "end", "boundary"))

                #chr_names <- c(paste0('chr', c('X', 'Y', c(22:1))))
                data <- regions %>% select(1:3) %>%
                    mutate(size = end - start)
                counts <- data %>%
                    group_by(chr) %>%
                    summarise(count = n())

                ggplot(data, mapping = aes(y = chr, x = size)) +
                    geom_boxplot() +
                    theme_bw() +
                  ggtitle("Distribution of region sizes per chromosome")
            })
        })

})
