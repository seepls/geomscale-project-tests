# visualize sampling in a polytope.
#visualize sampling in a polytope using Volesti



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


# sampling in a 3D sphere 
points = direct_sampling(n = 500, body = list("type" = "hypersphere", "dimension" = 3))
points_3D= as.data.frame(t(points))

make_3d_circ <- function(center = c(0,0),diameter = 2,non_dim='z'){
  radius = diameter/2
  circle_points <- seq(0,2*pi,length.out = 500)
  d1 <- center[1] + radius*cos(circle_points)
  d2 <- center[2] + radius*sin(circle_points)
  if(non_dim=='z'){
    return(data.frame(V1 = d1, V2 = d2, V3=0))
  }
  else if(non_dim=='y'){
    return(data.frame(V1 = d1, V2 = 0, V3=d2))
  }
  else{
    return(data.frame(V1 = 0, V2 = d1, V3=d2))
  }
  
}

cycle_3D = data.frame(matrix(ncol = 3, nrow = 0))

colnames(cycle_3D) <- c("V1", "V2", "V3")
all_data = rbind(sampled_data_3D,cycle_3D)
#plot
layout(
  plot_ly(x=all_data$V1, y=all_data$V2, z=all_data$V3, type="scatter3d", mode="markers", color = points_3D$color),
  title = "\n\nDirect sampling of 3D Ball"
)
