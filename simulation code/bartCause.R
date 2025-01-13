library(bartCause)
## fit a simple linear model
n <- 100L
beta.z <- c(.75, -0.5, 0.25)
beta.y <- c(.5, 1.0, -1.5)
sigma <- 2
set.seed(725)
x <- matrix(rnorm(3 * n), n, 3)
tau <- rgamma(1L, 0.25 * 16 * rgamma(1L, 1 * 32, 32), 16)
p.score <- pnorm(x %*% beta.z)
z <- rbinom(n, 1, p.score)
mu.0 <- x %*% beta.y
mu.1 <- x %*% beta.y + tau
y <- mu.0 * (1 - z) + mu.1 * z + rnorm(n, 0, sigma)
# low parameters only for example
fit <- bartc(y, z, x, n.samples = 100L, n.burn = 15L, n.chains = 2L)
summary(fit)
## example to show refitting under the common support rule
fit2 <- refit(fit, commonSup.rule = "sd")
fit3 <- bartc(y, z, x, subset = fit2$commonSup.sub,
              n.samples = 100L, n.burn = 15L, n.chains = 2L)

#### Employ BCF package
library(bcf)