#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application 
ui <- fluidPage(
  
  # Application title
  titlePanel("Intuiting covariance matrices (3x3 example)"),
  
  
  # Inputs begin
  #####
  sidebarLayout(
    sidebarPanel(
      
      sliderInput("n",
                  "No. of simulated datapoints",
                  min = 3,
                  max = 10000,
                  value = 250,
                  step = 1),      
      sliderInput("ab",
                  "A-B correlation",
                  min = -1,
                  max = 1,
                  value = 0.5,
                  step = .01),
      sliderInput("ac",
                  "A-C correlation",
                  min = -1,
                  max = 1,
                  value = 0.3,
                  step = .01),
      sliderInput("bc",
                  "B-C correlation",
                  min = -1,
                  max = 1,
                  value = 0.4,
                  step = .01),
      selectInput("covar", 
                  "Correlation or covariance?", 
                  c(correlation = "correlation", 
                    covariance = "covariance")),

      conditionalPanel(
        condition = "input.covar == 'covariance'",
        sliderInput("avar",
                    "Variance of 'A'",
                    min = .1,
                    max = 20,
                    value = 1,
                    step = .1),
        sliderInput("bvar",
                    "Variance of 'B'",
                    min = .1,
                    max = 20,
                    value = 1,
                    step = .1),
        sliderInput("cvar",
                    "Variance of 'C'",
                    min = .1,
                    max = 20,
                    value = 1,
                    step = .1)
      ),
      textOutput("byline")        # by line
    ),
    #####
    # Inputs end
    
    # Show outputs
    mainPanel( 
      fluidRow(

         tableOutput("covartable"), # Table for covars, a bit misaligned at this time
         plotOutput("EigenPlot"),   # Table for the (3) eigenvalues
         plotOutput("ABPlot")       # Facet plot similar to pairs()
      ))
  )
)

# Define server logic required
server <- function(input, output) {
  library(mvtnorm)
  library(ggplot2)
  library(GGally)
  library(xtable)
  
  get_covar_matrix <- reactive({
    # Constructing the covariance matrix based on correlations (and variances if applicable)
    
    if (input$covar == "correlation") {
      sig <- matrix( c( 1         , input$ab  , input$ac,
                        input$ab  , 1         , input$bc,
                        input$ac  , input$bc  , 1),
                     nrow = 3, byrow = TRUE)
    } else {
      #  putting these into vars for brevity and convert var to sd
      avar <- sqrt(input$avar)
      bvar <- sqrt(input$bvar)
      cvar <- sqrt(input$cvar)
      
      # (Pearson's) correlation = cov(XY)/(sigma(x)*sigma(y)) 
      sig <-  matrix( c( 1 * avar^2                , input$ab * avar * bvar  , input$ac * avar * cvar,
                         input$ab * avar * bvar  , 1 * bvar^2                , input$bc * bvar * cvar,
                         input$ac * avar * cvar  , input$bc * bvar * cvar  , 1 * cvar^2),
                      nrow = 3, byrow = TRUE)
      
    }
    # Works without return() but keeping here for readability (& sanity)
    return(sig)
    })
  
  get_simulated_data <- reactive({
    set.seed(123)
    sig <- get_covar_matrix()
    if ( min( eigen( sig)$values) > 0) {
      # cov. mat is pos. semidefinite
      dat <- rmvnorm( input$n, mean = c(0,0,0), sigma = sig)
      if ( input$covar == "correlation" ) { 
        cn <-  c("A (standardized)","B (standardized)","C (standardized)")
      } else {
        cn <- c("A (centralized - mean of A = 0)","B (centralized - mean of B = 0)","C  (centralized - mean of C = 0)")
      }
      colnames(dat) <- cn
      return(dat)
    } else {
      return(NA)
    }
  })
  
  output$byline <- renderText({"by Marton Kiss (mrkmarton@gmail.com) \n https://github.com/MartynK/Intuiting_Covar"})
  
  output$EigenPlot <- renderPlot({
    # a base R plot with the eigenvalues (to check positive semidefinity)
    
    sig <- get_covar_matrix()
    
    eigs <- data.frame( eigs = eigen(sig)$values)
    
    # eigs$eigs[3] <- -5 # !!!!!! for testing ONLY
    
    eigs$col <- sapply(1:3, function(x) {ifelse( eigs$eigs[x] > 0, "blue", "red")})

    ggplot(eigs, aes( x = 1:3,  y = eigs, fill = col)) +
      geom_bar(stat = "identity", 
               width = 0.7, 
               position = position_dodge(width=0.4)) +
      theme_bw() %+replace%
      theme(axis.line.x = element_blank(), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(angle=90),
            plot.title = element_text(hjust = .5, size = 16, face = "bold"),
            legend.position = "none") +
      xlab(label="") +
      ylab("Eigenvalues") +
      ggtitle("Waterfall plot; (all 3) eigenvalues of the covariance matrix") +
      scale_fill_manual(values=c("dodgerblue","darkorchid"))
  
    })
    
  output$ABPlot <- renderPlot({
    
    dat <- get_simulated_data()

    if (is.na(dat[1]) == FALSE) {

      pl  <- ggpairs(as.data.frame( dat)) +
        ggtitle("Simulated data and their correlations") +
        theme_bw() +
        theme(title = element_text(face="bold", size=16, hjust = .5),
              plot.title = element_text(face="bold", size=16, hjust = .5))
        
      if (input$covar == "covariance") {
        # Adding emphasis to note that only the scales change
        pl <- pl +
                theme(axis.text.x = element_text(face="bold", color="purple", size=16),
                      axis.text.y = element_text(face="bold", color="purple", size=16)
                      )
      }
      # print plot
      (pl)
    }
  })     
  
    output$covartable <- renderTable({ 
      dat <- get_simulated_data()
      if ( is.na(dat[1]) == FALSE) {
        tab <- get_covar_matrix()
      } else {
       tab <- matrix(rep(NA,9),ncol=3)
      }
      colnames(tab) <- c("    A    ","    B    ","    C    ") # Extra spaces were the easiest to control column width with... without success
      rownames(tab) <- c("  A  ","  B  ","  C  ") # doesn't work...
      xtable(tab)
      },
      align = "c",
      caption = "Covariance matrix of ABC",
      caption.placement = "top")
    

  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
  
  