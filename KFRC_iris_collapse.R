source("simstudy_src.R")
data(iris)
X <- iris[, c(1:4)]

R_good <- compute_dissimilarity(scale(X), bw = c(1, 1, 1, 1), kernel = "sum")
R_bad <- compute_dissimilarity(scale(X), bw = c(0.5, 1, 5, 0.5), kernel = "sum")

sds_good <- c()
sds_bad <- c()
for (i in seq(2, 50, length.out = 49)){
  FRC_good <- FRC(R_good^2, k = 3, m = i)
  FRC_bad <- FRC(R_bad^2, k = 3, m = i)
  U_good <- FRC_good$U
  U_bad <- FRC_bad$U
  sds_good <- c(sds_good, mean(apply(U_good, 1, sd)))
  sds_bad <- c(sds_bad, mean(apply(U_bad, 1, sd)))
}

plot(x = seq(2, 50, length.out = 49), y = sds_good, type = 'l', lwd = 2,
     col = 'blue',
     xlab = "Fuzzifier (m)",
     ylab = "Average row variability")
lines(x = seq(2, 50, length.out = 49), y = sds_bad, type = 'l', lwd = 2,
      col = 'red')
abline(h = 0, lty = 2, col = 'black')

legend("topright", legend = c("No collapse", "Uniform collapse"),
       lwd = 2, col = c("blue", "red"))