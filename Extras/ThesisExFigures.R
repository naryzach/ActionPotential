library(SobolSequence)

unifm_2d <- function(n, seed) {
  set.seed(seed)
  return(replicate(2, runif(n)))
}

sobol_2d <- function(n, seed) {
  return(randtoolbox::sobol(n, dim = 2, scrambling = 3, seed = seed))
}

latin_2d <- function(n, seed) {
  set.seed(seed)
  return(lhs::randomLHS(n, k = 2))
}

par(mfrow = c(1, 3))
plot(latin_2d(1000, 2019), main = "Latin Hypercube Method", xlab = '', ylab = '', cex = 2, col = "blue")
plot(sobol_2d(1000, 2019), main = "Sobol Sequence", xlab = '', ylab = '', cex = 2, col = "red")
plot(unifm_2d(1000, 2019), main = "Random Uniform", xlab = '', ylab = '', cex = 2, col = "black")

x_vals = c(-10:10)
ident = data.frame(
  x = x_vals,
  y_2 = rep(5, length(x_vals)),
  y = x_vals^2, 
  value = abs(x_vals)
)
ident <- (ident[order(ident$value),])

par(mfrow = c(1, 1))
plot(ident[,"x"],ident[,"y"],
     main="Identifiable",
     xlab="Estimated Parameter",
     ylab="Objective Function Value",
     col=colorRampPalette(c("blue","red"))(nrow(ident)) )

plot(ident[,"x"],ident[,"y_2"],
     main="Not Identifiable",
     xlab="Estimated Parameter",
     ylab="Objective Function Value",
     col="blue")
