# Sample approximate uniformly points from a random generated polytope using the implmented in volesti random walks, for various walk lengths.
# For each sample compute the PSRF and report on the results.

library(ggplot2)
library(volesti)
library(devtools)

# p dimension of the estimation problem.
# m number of chains.
# epsilon relative precision level.
# delta desired delta value - the cutoff for potential scale reduction factor. If specified, then the corresponding \code{epsilon} is returned.
# alpha significance level for confidence regions for the Monte Carlo estimators.

target.psrf <- function(p, m, epsilon = .05, delta = NULL, alpha=.05){
  if(is.null(delta)){
    Tee <- as.numeric(minESS(p, epsilon, alpha = alpha)) #min effective sample size
  }
  
  if(is.null(delta) == FALSE){
    Tee <- M <- m/((1+delta)^2 - 1)
    epsilon <- as.numeric(minESS(p, eps = epsilon, ess = M))
  }
  
  del <- sqrt(1 + m/Tee) - 1
  arr <- 1 + del  
  list(psrf = arr, epsilon = epsilon)
}



for (step in c(1,20,100,150,500)){
  for (walk in c("BiW")){
    P <- gen_cube(100, 'H')
    points1 <- sample_points(P,n=1000,random_walk = list("walk" = walk, "walk_length" = step))
    g<-plot(ggplot(data.frame( x=points1[1,], y=points1[2,] )) +
              geom_point( aes(x=x, y=y, color=walk)) + coord_fixed(xlim = c(-1,1),
                                                                   ylim = c(-1,1)) + ggtitle(sprintf("walk length=%s", step, walk)))

    print(target.psrf(p=3,m=step))
  }
  
}

# output : 
'''
c = 1 
$psrf
[1] 1.000062

$epsilon
[1] 0.05

c = 20 
$psrf
[1] 1.00123

$epsilon
[1] 0.05

c= 100
$psrf
[1] 1.006137

$epsilon
[1] 0.05

c= 150
$psrf
[1] 1.009191

$epsilon
[1] 0.05

c=500
$psrf
[1] 1.030317

$epsilon
[1] 0.05
'''
