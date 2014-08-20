library(shiny)
library(igraph)

shinyUI(fluidPage(
  sidebarLayout(
    sidebarPanel( h4("Logistic Map"),
                  withMathJax("$$x_{t+1}=r x_t (1-x_t)$$"),
                  h5("Model Parameters"),
                  sliderInput("r", withMathJax("$$r$$"), min=0, max=4, value = 4, ticks=FALSE, step=0.001),
                  sliderInput("x0", withMathJax("$$x_0$$"), min=0, max=1, value = 0.1),
                  sliderInput("iterations", withMathJax("Model Iterations (i.e. edges)"), min=1000, max=10000, value = 1000),
                  h5("Graph Layout Parameters"),
                  sliderInput("bins", "Vertices (i.e. bins)", min = 10, max = 1000, value = 500),
                  checkboxInput("mult_edges", label = "Keep Duplicated Edges", value = TRUE)
    ),
    mainPanel(
      tabsetPanel(type = "tabs", 
                  tabPanel("Graph", 
                           fluidRow(
                           column(5,
                           h5("Logistic Map Graph"),
                           plotOutput('plotGraph')
                           ),
                           column(6,
                                  h5("Metrics"),
                                  dataTableOutput("metrics")
                           )
                           )
                  ),
                  tabPanel("Plot", 
                           h5("Degree Distribution"),
                           plotOutput('plotDegDist')
                  ),
                  tabPanel("About",
                           h5("Help"),
                           "The Logistic Map is a simple non-linear recursive relation that leads to chaotic dynamics, depending on the value of a control parameter r.",
                           "Matematically, it is defined as",
                           withMathJax("$$x_{t+1}=r x_t (1-x_t)$$ where $$x_t\\in [0,1], ~\\forall t$$"),
                           "Here, we use this recursive relation to create a graph.",
                           "First, we bin the domain [0,1], such that each bin will correspond to a vertex in the graph.",
                           "Second, two vertices A and B will be connected by a directed edge if a number Xt in bin A, is mapped to a number Xt+1 in bin B, after one iteration of the Logistic Map.",
                           br(),
                           "The side bar allows the user to chage the value of the Logistic Map parameters r (i.e. control chaotic dynamics), x0 (i.e. starting point) and the number of recursive iterations (i.e. edges).",
                           "The user can also change the number of vertices in the graph (i.e. bins), and whether duplicated edges should be kept or not.",
                           "The results are presented in the form of a graph, a set of quantitative graph metrics and also two plots depicting the degree distribution of the respective graph",
                           h5("Author"),
                           "Marcelo S. Zanetti <ms.zanetti@ufma.br>",
                           br(),
                           "August 20, 2014"
                  )
      )
    )
  )
)
)