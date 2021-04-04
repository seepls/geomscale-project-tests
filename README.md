#  Visualize sampling in a polytope 
```R
library(volesti)
library(ggplot2)
library(gridExtra)
library(plotly)

#sampling in a 3D cube
x1<-runif(1000, min = -1, max = 1)
x2<-runif(1000, min = -1, max = 1)
x3<-runif(1000, min = -1, max = 1)

#plot
layout(
  plot_ly(x=x1, y=x2, z=x3, type="scatter3d", mode="markers"),
  title = "\n\nSampling of 3D cube"
)

```
![png](Sampling_cube.png)
#  Sample approximate uniformly points from a random generated polytope using the implmented in volesti random walks, for various walk lengths. For each sample compute the PSRF and report on the results.
# Implement Gibbs sampler for uniform sampling from a polytope in C++ following the code structure of volesti.
