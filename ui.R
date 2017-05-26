###################################################################
# R code to accompany "Impact of non-uniform correlation structure
#    on sample size and power in multiple-period cluster randomised
#    trials", by J Kasza et al.
# 
# This file contains the UI for the Shiny app 
###################################################################

library(shiny)
library(plotly)

shinyUI(pageWithSidebar(
  h3("R Shiny App to accompany \"Non-uniform correlation structures\", Kasza et al."),
  sidebarPanel(
    sliderInput("T",
                "Number of periods, T:",
                min = 2,
                max=20,
                step = 1,
                value = 4),
    numericInput("nclust",
                 "Number of clusters assigned to each treatment sequence:",
                 min = 1,
                 max=50,
                 step = 1,
                 value = 1),
    numericInput("m",
                 "Number of subjects in each cluster-period, m:",
                 min = 1,
                 max=1000,
                 step = 1,
                 value = 500),
    sliderInput("rho0",
                label = HTML(paste("Intra-cluster correlation, &rho;",tags$sub(0),":", sep="")),
                #"rho 0:"
                min = 0,
                max = 0.2,
                step = 0.005,
                value = 0.035),
    helpText(HTML(paste("Note: calculations assume &sigma;",tags$sup(2),tags$sub("CP"), "+", "&sigma;",
                        tags$sup(2),tags$sub("e"), "=1",  " so  &rho;",tags$sub(0), 
                        "=&sigma;",tags$sup(2),tags$sub("CP"), "/(&sigma;",tags$sup(2),tags$sub("CP"), "+", "&sigma;",
                        tags$sup(2),tags$sub("e"), ") =&sigma;",tags$sup(2),tags$sub("CP"), sep=""))),
    
   
    
    conditionalPanel(condition="input.conditionedPanels==1",
                     
                     selectInput("select", label = ("What would you like to calculate?"), 
                                 choices = list("Variance" = 1, 
                                                "Design effect" = 3, "Power" = 4), selected = 1),
                     conditionalPanel(
                       condition = "input.select == '4'",
                       sliderInput("effsize",
                                   "Effect size:",
                                   min=0.05,
                                   max=1,
                                   step=0.05, 
                                   value = 0.2))
                     
    ),
    
    conditionalPanel(condition="input.conditionedPanels==1 | input.conditionedPanels == 2",
                     radioButtons("HooperYN", "Include Hooper/Girling model?",
                                  choices = list("No" =1, "Yes"=2), selected = 1),
                     conditionalPanel(
                       condition = "input.HooperYN == '2'",
                       sliderInput("rCons",
                                   "alpha value for Hooper/Girling model:",
                                   min = 0,
                                   max = 1,
                                   step = 0.005,
                                   value = 0.1))
    ),                 
    
    
    conditionalPanel(condition="input.conditionedPanels==3",
                     
                     selectInput("design", label = ("Which cluster randomised design would you like to consider?"), 
                                 choices = list("Stepped wedge" = 1, "Parallel" = 2,
                                                "Parallel with baseline" = 3, 
                                                "CRXO" = 4), selected = 1),
                     
                     sliderInput("maxr",
                                 "Maximum value of decay (1-r):",
                                 min = 0,
                                 max=1,
                                 step = 0.05,
                                 value = 1)
                     
                     
    )
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Impact of exponential decay", value=1, plotlyOutput("varplotlyexp"),
               textOutput("ExpDecayExplan"),
               textOutput("Hooperincltext")), 
      tabPanel("Exponential decay vs. Hussey and Hughes", value=2,  plotlyOutput("HHrelvarplot"),
               textOutput("HHrelvarexplan")),
      
      tabPanel("Exponential decay vs. Hooper/Girling", value=3, plotOutput("Contourplot", height=600, width=600),
               textOutput("HGrelvarexplan")),
      
      
      tabPanel("Design matrices", value=4, textOutput("text1"),
               tableOutput("SWxmat"),
               textOutput("text2"),
               tableOutput("pllelxmat"),
               textOutput("text3"),
               tableOutput("pllelbxmat"),
               textOutput("text4"),
               tableOutput("crxoxmat") )
      , id = "conditionedPanels"
    )
  )
))

