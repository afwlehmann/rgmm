#
# simpleGMM.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#

library(rgmm)

# Compute GMM.
X <- rbind(rmvnorm( 9000, mu=c(30,4,5), Sigma=c(6,7,8)),
           rmvnorm( 8000, mu=c(5,4,5), Sigma=c(17,4,11)),
           rmvnorm(12000, mu=c(0,25,4), Sigma=c(22,5,25)))
colnames(X) <- c("Foo", "Bar", "Baz")
model <- computeGMM(X, K=3)

# Multidimensional contour plot.
dev.new()
plot(model, dataset=X)

# Perspective plot of the first two variables.
dev.new()
marginOneTwo <- marginal(model, c(1,2))
x <- y <- seq(-15, 45, length.out=100)
xy <- as.matrix(expand.grid(x, y))
z <- matrix(dgmm(xy, marginOneTwo), ncol=length(x))
persp(x, y, z, theta=30, phi=30)

# Curve plot of the second variable.
dev.new()
marginTwo <- marginal(model, 2)
plot(function(x) dgmm(x, marginTwo), -15, 45, n=250, ylab="p(x)")
