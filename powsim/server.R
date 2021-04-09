#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(vcfR)
library(stringr)
library(gwasim)
mem_limit <- 1024 # 1GB
mem_limit_bytes <- mem_limit * 1024^2
options(shiny.maxRequestSize = mem_limit_bytes) # set max upload size

# Define server logic
shinyServer(function(input, output) {
       values <- reactiveValues(coords = NULL, maf = NULL, N = NULL, max_q2 = NULL)

        output$input_info <-
            renderTable({
            req(input$input_file)
            tryCatch(
                {
                if (grepl("\\.vcf$", input$input_file$datapath) | grepl("\\.vcf.gz$", input$input_file$datapath)) {
                    vcf_data <- vcfR::read.vcfR(file = input$input_file$datapath, nrows = 500, limit = mem_limit_bytes)
                    shp <- dim(vcf_data@gt)
                    values$N <- shp[2]
                    df <- data.frame(c("Filename","N markers", "N individuals"), c(input$input_file$datapath, shp))
                    colnames(df) <- c("Feature", "Value")
                    markers <- vcf_data@fix[,'ID'] # Extract marker names
                    coords <- vcf_data@fix[, 'POS'] # Extract genomic marker coordinates
                    to_remove <- which(str_detect(markers, ";")) # detect weird markers (multi-marker)
                    markers <- markers[-to_remove]
                    coords <- coords[-to_remove]
                    G <- gwasim::get_genotypes(x = vcf_data, marker_names = markers) %>%
                        gwasim::fix_allele_encoding()
                    rm(vcf_data)
                    values$coords <- coords
                    values$maf <- colSums(G) / (2 * dim(G)[1])
                    rm(G)
                    #gc()
                } else if (grepl("\\.rdat$", input$input_file$datapath)){
                    load(input$input_file$datapath)
                    values$N <- data$N_inds
                    values$maf <- data$maf
                    values$coords <- data[['coords']]
                    df <- data.frame(c("Filename","N markers", "N individuals"), c(input$input_file$datapath, length(maf), N_inds))
                    colnames(df) <- c("Feature", "Value")
                } else {

                }
                },
                error = function(e) {
                    # return a safeError if a parsing error occurs
                    stop(safeError(e))
                }
            )
            return(df)
            }
        )

        output$maf_plot <- renderPlot({
            if (!is.null(values$coords) & !is.null(values$maf)) {
                plot(values$coords, values$maf, ylim=c(0, .5), las=1, cex.axis = .5, type='n', main="Minor Allele Frequency", xlab='bp', ylab='maf', bty='n')
                grid(lty=2, col='grey')
                points(values$coords, values$maf, pch=19, cex=.5)
            }
        })

        output$max_q2_plot <- renderPlot({
            if (!is.null(values$coords) & !is.null(values$maf)) {
             max_q2 <- 2 * values$maf * (1-values$maf) * input$beta^2
             max_q2 <- pmin(max_q2, 1)
             values$max_q2 <- max_q2
             plot(values$coords, max_q2, ylim=c(0, 1), las=1, cex.axis = .5, type='n', main = "Maximal q2", xlab='bp', ylab = 'q2', bty='n')
             grid(lty=2, col='grey')
             points(values$coords, max_q2, pch=19, cex=.5)
            }
        })

        output$power_plot <- renderPlot({
            if (!is.null(values$coords) & !is.null(values$maf) & !is.null(values$N) & !is.null(values$max_q2)) {
                N <- values$N
                max_q2 <- values$max_q2
                NCP <- (N * max_q2^2)/(1 - max_q2^2)
                thr <- qchisq(input$p_val, df = 1, lower.tail = F)
                pow <- pchisq(thr, df  = 1, lower.tail = F, ncp = NCP)
                pow <- pmin(pow, 1)
                plot(values$coords, pow, ylim = c(0, 1), las=1, cex.axis = .5, type='n', main='Statistical Power', xlab='bp', ylab='Power', bty='n')
                grid(lty=3, col='grey')
                points(values$coords, pow, pch=19, cex=.5)
            }
        })


})
