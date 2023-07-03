#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(Ryacas)
library(mathjaxr)
source("helpers.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("SIR-, SIRS-, and SEIR-models"),
    h3("N=1.5 millions, S(0)=N-1, I(0)=1, R(0)=0"),
    br(),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          h4("Parameters"),
            sliderInput("beta",
                        "beta:",
                        min = 0,
                        max = 1,
                        value = 0.5),
            br(),
            sliderInput("alpha",
                        "alpha:",
                        min = 0,
                        max = 1,
                        value = 1/10),
            br(),
            sliderInput("gamma",
                        "gamma/gamma_e:",
                        min = 0,
                        max = 1,
                        value = 1/30),
            br(),
            #checkboxInput("sirs",
            #              "Show SIRS-model",
            #              value = FALSE),
            #br(),
            #checkboxInput("seir",
            #              "Show SEIR-model",
            #              value = FALSE),
            radioButtons("model", label = h4("Model"),
                         choices = list("SIR" = 1, "SIRS" = 2, "SEIR" = 3), 
                         selected = 1),
            br(),
          h4("Show data (for I) for the 4th wave in Munich"),
            checkboxInput("realI",
                          "Yes!",
                          value = FALSE)
          
        ),

        # Show a plot of the generated distribution
        mainPanel(
           h4("Model plot"),
           plotOutput("distPlot"),
           withMathJax(),
           h4("Model equations:"),
           uiOutput("model")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$distPlot <- renderPlot({
    
    plot_model(parms = c(beta=input$beta, alpha=input$alpha, gamma=input$gamma), model=input$model, real=input$realI)
    
  })
  
  output$model <- renderUI({
    if (input$model == 1) {
      withMathJax(
        helpText('$$ \\dot{S} = - \\beta \\frac{I}{N} S $$'),
        helpText('$$\\dot{I} = \\beta \\frac{I}{N} S - \\alpha I $$'),
        helpText('$$ \\dot{R} = \\alpha I $$'),
        sprintf('$$ \\text{Basic reproduction number} \\; \\textit{R}_0 = \\frac{\\beta \\cdot S(0)}{N \\cdot \\alpha} = \\frac{%g \\cdot S(0)}{N \\cdot %g} = %g$$', 
                input$beta, input$alpha, input$beta * 1499999/ (input$alpha * 1500000))
      )
    } else if (input$model == 2) {
      withMathJax(
        helpText('$$ \\dot{S} = - \\beta \\frac{I}{N} S + \\gamma R$$'),
        helpText('$$\\dot{I} = \\beta \\frac{I}{N} S - \\alpha I $$'),
        helpText('$$ \\dot{R} = \\alpha I - \\gamma R $$'),
        sprintf('$$ \\text{Basic reproduction number} \\; \\textit{R}_0 = \\frac{\\beta \\cdot S(0)}{N \\cdot \\alpha} = \\frac{%g \\cdot S(0)}{N \\cdot %g} = %g$$', 
                input$beta, input$alpha, input$beta * 1499999/ (input$alpha * 1500000))
      )
    } else {
      withMathJax(
        helpText('$$ \\dot{S} = - \\beta \\frac{I}{N} S$$'),
        helpText('$$ \\dot{E} = \\beta \\frac{I}{N} S - \\gamma_e E $$'),
        helpText('$$\\dot{I} = \\gamma_e E - \\alpha I $$'),
        helpText('$$ \\dot{R} = \\alpha I$$'),
        sprintf('$$ \\text{Basic reproduction number} \\; \\textit{R}_0 = \\frac{\\beta \\cdot S(0)}{N \\cdot \\alpha} = \\frac{%g \\cdot S(0)}{N \\cdot %g} = %g$$', 
                input$beta, input$alpha, input$beta * 1499999/ (input$alpha * 1500000))
      )
    }
  })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
