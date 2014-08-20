library(shiny)

shinyServer(function(input, output) {
  my_graph_reactive <- reactive({
    
    my_graph         <- graph.empty(0,directed=TRUE)
    my_graph         <- add.vertices(my_graph, input$bins)
    V(my_graph)$name <- c(1:input$bins)
    bins             <- seq(0,1.0001,length=(1+input$bins))
    x0               <- input$x0
    r                <- input$r
    n_iter           <- input$iterations
    for(iteration in c(1:n_iter))
    {
      x1       <- r*x0*(1-x0)
      myBin    <- findInterval(c(x0, x1), bins)
      my_graph <- add.edges(my_graph, myBin, directed=TRUE)
      x0       <- x1
    } 
    if(!input$mult_edges)
    {
      my_graph <- simplify(my_graph, remove.multiple = TRUE, remove.loops=FALSE)
    }
    V(my_graph)$degree <- degree(my_graph, V(my_graph), mode="all")
    
    return (my_graph)
  })
  
  output$metrics = renderDataTable({
    library(ggplot2)
    my_table <- data.frame(metric=character(), value=numeric())
    my_table <- rbind(my_table, data.frame(metric="Mean Degree", value=mean(V(my_graph_reactive())$degree)))
    x        <- centralization.degree (my_graph_reactive(), mode = "all")
    my_table <- rbind(my_table, data.frame(metric="Degree Centralization", value=x$centralization/x$theoretical_max))
    x        <- centralization.closeness (my_graph_reactive(), mode = "all")
    my_table <- rbind(my_table, data.frame(metric="Closeness Centralization", value=x$centralization/x$theoretical_max))
    x        <- centralization.evcent (my_graph_reactive(), directed = TRUE)
    my_table <- rbind(my_table, data.frame(metric="Eigenvector Centralization", value=x$centralization/x$theoretical_max))
    my_table <- rbind(my_table, data.frame(metric="Diameter", value=diameter(my_graph_reactive(), directed = TRUE, unconnected = TRUE, weights = NULL)))
    x        <- V(my_graph_reactive())$degree
    my_table <- rbind(my_table, data.frame(metric="Vertices in Largest Component", value=as.integer(length(x[x>0]))))
    my_table
  })
  
  output$plotDegDist <- renderPlot({
    par(mfrow=c(1,2))
    plot(density(V(my_graph_reactive())$degree), main="")
    plot.ccdf(V(my_graph_reactive())$degree, xlab="Degree", ylab="CCDF(Degree)")
  })

  output$plotGraph <- renderPlot({
    n_iter <- 1000 #input$lout_iterations
    L      <- generateLayout(my_graph_reactive(), niter=n_iter, scalingf=3, smallestN=2)
    e_W    <- 0.7   #input$edge_width
    a_S    <- 0.3   #input$arrow_size
    v_S    <- L[,3] #input$vertex_size
    v_L    <- ""
    plot(my_graph_reactive(), layout=L[,c(1,2)], edge.width=e_W, edge.arrow.size=a_S, vertex.label=v_L, vertex.size=v_S)
  })
})

# AUXILIARY FUNCTIONS
generateLayout <- function(g, seed=1, niter=1000, scalingf=1, smallestN=1)
{
  set.seed(seed)
  
  L     <- layout.kamada.kawai(g, niter=niter)
  L[,1] <- (L[,1]-min(L[,1]))/(max(L[,1])-min(L[,1]))*2-1
  L[,2] <- (L[,2]-min(L[,2]))/(max(L[,2])-min(L[,2]))*2-1
  
  L     <- cbind(L, degree(g, V(g)))
  maxk  <- max(L[,3])
  L[,3] <- smallestN+(scalingf*L[,3]/maxk)
  return(L)
}

plot.ccdf <- function (data, log = "xy", xlab = "x", ylab = "CCDF(x)", xlim
                       = F, ymin = 0, sample.it=F, sample.count = 100, col = 1, cex = 1) 
{
  if ( length(dim(data)) > 0) {
    cat ("Attention input is a matrix! Usign all values in it!")
  }
  x <- sort(data)
  
  if ( (log == "xy") || ( log == "y")) {
    if ( ymin == 0) {
      min = ceiling(log10(length(x)))
    } else {
      min = ceiling(log10(1/ymin))
    }
    ylab.at = (1/(10^(0:min)))
    ylab.lab = ylab.at
    yrange = c(min(ylab.at),1)
  } else {
    ylab.at = 0:10/10
    ylab.lab = ylab.at
    yrange = c(0,1)
  }
  
  if ( length(xlim) != 2) {
    if ( (log == "xy") || (log == "x")) { 
      minx.log = x[ x>0 ][1]
      xlim = c(minx.log, max(x)) 
    } else {
      xlim = c(min(x), max(x))
    }
  }
  
  cat( "yrange: ", yrange, "| xrange: ", xlim, "\n")
  plot(1,1, log=log, xlab = xlab, ylab = ylab, 
       ylim = yrange, xlim = xlim, axes = F, type = 'n')
  add.ccdf(data, col = col, sample.it = sample.it, sample.count = sample.count, cex = cex)
  axis(1)
  axis(2, at = ylab.at, labels = ylab.lab)
  box()
  grid()
}

add.ccdf <- function (data, col = 2, sample.it = F, sample.count = 100, cex = 1) 
{
  if ( length(dim(data)) > 0) {
    cat ("Attention input is a matrix! Usign all values in it!")
  }
  
  x <- sort(data)
  y <- 1-((1:(length(x))-1)/length(x))
  
  xmin <- min(x)
  xmax <- max(x)
  ymin <- min(y)
  ymax <- max(y)
  
  if (sample.it) {
    if ( par("ylog") == T ) {
      ymarks <- 10^(seq(from=log10(ymin), to = log10(ymax), length.out=sample.count))
    } else {
      ymarks <- seq(from=ymin, to = ymax, length.out=sample.count)
    }
      
    if ( par("xlog") == T ) { 
      xsbr <- 10^(seq(from=log10(xmin), to = log10(xmax), length.out=sample.count))
    } else {
      xsbr <- seq(from=xmin, to = xmax, length.out=sample.count)
    }
    xmarks <- cumsum(table(cut(x, breaks = xsbr, labels = F)))
  }
  
  if ( sample.it ){
    points(x= quantile(data, probs = 1-ymarks, type =2), y = ymarks, col = col, pch = col, cex = cex)
    points(x[xmarks], y[xmarks], col=col, pch = col, cex = cex)
  } else {
    points(x,y,col = col, pch = col, cex = cex)
  }
}