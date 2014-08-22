---
title       : Chaotic Map App
subtitle    : From Logistic Map to a Graph
author      : Marcelo Serrano Zanetti
job         : Nerd
framework   : io2012        # {io2012, html5slides, shower, dzslides, ...}
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      # 
widgets     : [mathjax]     # {mathjax, quiz, bootstrap}
mode        : selfcontained # {standalone, draft}
knit        : slidify::knit2slides
---

## The Logistic Map: Definition

$$x_{t+1}=r x_t(1-x_t), ~x_t \in [0,1],~ \forall t$$

1. a simple non-linear recursive relation 
2. can behave chaotically depending on parameter $r$

  a. chaotic dynamic systems are sensitive to small perturbations

  b. two trajectories starting at nearby points diverge exponentially


--- .class #id 

## The Logistic Map: Dynamics depending on $r$

<img src="assets/fig/unnamed-chunk-1.png" title="plot of chunk unnamed-chunk-1" alt="plot of chunk unnamed-chunk-1" style="display: block; margin: auto;" />

--- .class #id 

## From Logistic Map to a Graph

1. we divide the domain $[0,1]$ in $n$ bins, which we draw as vertices
2. we define an initial condition $x_0=l, ~l\in [0,1]$ 
3. for some $t$, we apply the logistic map to $x_t$, to obtain $x_{t+1}$

  a. if $x_t$ belongs to bin $b_k$ and $x_{t+1}$ to bin $b_m$, we draw a directed edge connecting $b_k$ and $b_m$

--- .class #id

## Logistic Graph!

<img src="assets/fig/unnamed-chunk-2.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />
by the way, this is R!
