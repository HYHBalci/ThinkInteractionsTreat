library(MASS)

simulate_data <- function(n_samples, p_main, p_noise_main, interaction = TRUE, treatment = TRUE, binary = FALSE, heterogeneity = FALSE, seed = 123) {
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Define the number of samples and covariates
  n <- n_samples
  p <- p_main + p_noise_main  # Total number of covariates
  
  # Generate moderately collinear continuous covariates
  mean_vector <- rep(0, p)
  cov_matrix <- matrix(0.5, nrow = p, ncol = p)
  diag(cov_matrix) <- 1  # Set diagonal to 1 for variances
  X <- MASS::mvrnorm(n = n, mu = mean_vector, Sigma = cov_matrix)
  
  # Assign coefficients for main effects
  main_effects <- rep(0, p)
  for(i in 1:p_main){
    main_effects[i] <- 2*log(1+i)
  }
  #Heterogeneity effects of treatment
  main_heterogeneity <- rep(0, p)
  p_heterogeneous <- floor(p/2)
  
  for(j in 1:p_heterogeneous) {
    main_heterogeneity[j] <- log(1+j)
  }
  # Initialize interaction effects
  interaction_effects <- matrix(0, nrow = p, ncol = p)
  non_zero_main_effects <- 1:p_main
  if (interaction) {
    # Define non-zero interactions
    non_zero_interactions <- list(
      "main_main" = 6,  # Interactions of two non-zero main effects
      "main_absent" = 4,  # Interactions between a non-zero and an absent main effect
      "absent_absent" = 3  # Surprising interactions between two absent main effects
    )
    
    # Interactions between two non-zero main effects
    main_main_pairs <- t(combn(non_zero_main_effects, 2))
    chosen_main_main <- 1:non_zero_interactions$main_main  # Select the first set of pairs deterministically
    for (i in chosen_main_main) {
      pair <- main_main_pairs[i, ]
      interaction_effects[pair[1], pair[2]] <- 0.3  # Fixed interaction value
      interaction_effects[pair[2], pair[1]] <- interaction_effects[pair[1], pair[2]]  # Symmetric interaction
    }
    
    # Interactions between a non-zero and an absent main effect
    absent_main_indices <- setdiff(1:p, non_zero_main_effects)
    chosen_main_absent <- as.matrix(expand.grid(non_zero_main_effects, absent_main_indices))
    chosen_main_absent <- chosen_main_absent[1:non_zero_interactions$main_absent, ]  # Select the first set 
    for (i in 1:nrow(chosen_main_absent)) {
      pair <- chosen_main_absent[i, ]
      interaction_effects[pair[1], pair[2]] <- 0.3  # Fixed interaction value
      interaction_effects[pair[2], pair[1]] <- interaction_effects[pair[1], pair[2]]  
    }
    
    # Surprising interactions between two absent main effects
    chosen_absent_absent <- t(combn(absent_main_indices, 2))
    chosen_absent_absent <- chosen_absent_absent[1:non_zero_interactions$absent_absent, ]  
    for (i in 1:nrow(chosen_absent_absent)) {
      pair <- chosen_absent_absent[i, ]
      interaction_effects[pair[1], pair[2]] <- 0.3  
      interaction_effects[pair[2], pair[1]] <- interaction_effects[pair[1], pair[2]]  
    }
  
  if (treatment) {
    if(binary){
      treatment_effect <-  rbinom(n, 1, 0.5)
    } else {
      treatment_effect <-  rnorm(n, 0, 1)
    }
  # Simulate response variable
  if(heterogeneity){
    linear_predictor <- X %*% main_effects + rowSums((X %*% interaction_effects) * X) + treatment_effect + ((diag(treatment_effect) %*% X ) %*% main_heterogeneity) 
  } else {
    linear_predictor <- X %*% main_effects + rowSums((X %*% interaction_effects) * X) + treatment_effect
  }
  
  response <- linear_predictor + rnorm(n, 0, 0.25)  # Add Gaussian noise
  
  # Compile the dataset
  simulated_data <- data.frame(response = response, X)
  # Add all interaction terms (including zero-coefficient interactions)
  interaction_cols <- list()
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      col_name <- paste0("X", i, ":X", j)
      simulated_data[[col_name]] <- X[, i] * X[, j]
      interaction_cols[[col_name]] <- interaction_effects[i, j]
    }
  }
  
  if (treatment) {
    simulated_data$treatment <- treatment_effect
  }
  
  return(list(
    simulated_data = simulated_data, 
    main_effects = main_effects, 
    interaction_effects = interaction_effects,
    interaction_columns = interaction_cols
  ))
}}
}
simulated_data <- simulate_data(n_samples = 200, p_main = 5, p_noise_main = 4, interaction = TRUE, treatment = TRUE, seed = 123)

  