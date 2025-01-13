library(grf)
library(caret)
source('simulation code/simulations.R')
p_main <- 5
p_noise_main <- 4
p <- p_main + p_noise_main

#simulate data
data <- simulate_data(n_samples = 3000, p_main = p_main, p_noise_main = p_noise_main, interaction = TRUE, treatment = TRUE, binary = TRUE, seed = 123)
data1 <- data$simulated_data

Y1 <- as.numeric(data1[,1])
A <- as.numeric(data1[,47])
X1 <- as.matrix(data1[,-c(1, ncol(data1))])
X1 <- X1[, 1:9]
#create test and train data set
train_index <- createDataPartition(Y1, p = 0.8, list = FALSE)

Y_train <- Y1[train_index]
Y_test <- Y1[-train_index]

A_train <- A[train_index]
A_test <- A[-train_index]

X_train <- X1[train_index, ]
X_test <- X1[-train_index, ]

#fit the model
causal_forest_model <- causal_forest(X_train, Y_train, A_train)
 
ate <- average_treatment_effect(causal_forest_model)
print(ate)

tau.hat <- predict(causal_forest_model, X_test, estimate.variance = TRUE)
sigma.hat <- sqrt(tau.hat$variance.estimates)

#Employing Bart (Thisfunction is limited please proceed to use BCF!)
library(bartCause)
X1 <- as.matrix(X1)
confounders_string <- paste(colnames(X1), collapse = " + ")
confounders_formula <- as.formula(paste("~ ", confounders_string))
X1 <- X1[, 1:9]
fit <- bartc(
  response = Y1,
  treatment = A,
  confounders = X1,
  method.rsp = "bart",
  method.trt = "bart",
  estimand = "ate"
)

summary(fit)

### BCF, Carvalho (2020)
library(bcf)
data <- simulate_data(n_samples = 500, p_main = p_main, p_noise_main = p_noise_main, interaction = TRUE, treatment = TRUE, binary = TRUE, seed = 123)
data1 <- data$simulated_data

Y1 <- as.numeric(data1[,1])
A <- as.numeric(data1[,47])
X1 <- as.matrix(data1[,-c(1, ncol(data1))])
X1 <- X1[, 1:9]
bcf_fit = bcf(Y1, A, X1, X1, rep(0.5, 500), nburn=2000, nsim=2000)

tau_post = bcf_fit$tau
tauhat = colMeans(tau_post)
print(tauhat)
mean(tauhat)