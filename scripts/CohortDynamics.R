# Gates Dupont
library(tidyverse)
library(deSolve)
library(progress)

# ---- Helper functions ----

# See first 5 x 5 of a matrix
peek <- function(mat) return(mat[1:5,1:5])

# Calculate r-star
calculate_r_star <- function(c, gR, g0, mR, mu, Tv, Tw, v, w) {
  
  # Calculate the terms inside the square root
  sqrt_term <- sqrt((c^2) * (mR^2) * gR * g0 * mu * Tv * Tw * (v*Tw + gR * w * Tv))
  
  # Calculate r_star according to the given formula
  r_star <- (c * g0 / mR) + (sqrt_term / ((mR^2) * (v * Tw + gR * w * Tv)))
  
  return(r_star)
}

# Calculate alpha
calculate_alpha <- function(v, mR, r, c, g0, gR, Tv, w, Tw){
  
  numerator1 <- v * (mR * r - c * g0)
  denominator1 <- c * gR * Tv
  
  numerator2 <- w * (mR * r - c * g0)
  denominator2 <- c * Tw
  
  result <- (numerator1/denominator1) + (numerator2/denominator2)
  return(result)
  
}


# ---- Model functions -----

# Population-level model function
within_host_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dV <- r*V - c*V*L
    dL <- g0 + gR*r*V - mR*L
    
    beta <- zeta*V
    
    return(list(c(dV, dL), beta = beta))
  })
}

# Population-level model function
population_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    N <- S + I
    
    dS <- N * (b - q * N) - beta*S*I - mu*S
    dI <- beta*S*I - mu*I - alpha*I
    
    # Return the rate of change
    return(list(c(dS, dI)))
  })
}


# ---- Multiscale setup ----

# Parameters (example values, adjust based on your system)
parameters <- c(b = 0.2,   # Birth rate
                mu = 1 / (20*365), # Death rate
                q = 0.002, # Density dependence on birth rate
                zeta = 0.2, 
                c = 0.5, 
                g0 = 0.3, 
                gR = 0.9, 
                mR = 1/21,
                mu = 1 / (20*365),
                Tv = 0.5,
                Tw = 0.5,
                v = 1,
                w = 1) # Natural leukocyte death rate

# r-star
parameters["r"] <- calculate_r_star(
  c = parameters[["c"]],
  gR = parameters[["gR"]],
  g0 = parameters[["g0"]],
  mR = parameters[["mR"]],
  mu = parameters[["mu"]],
  Tv = parameters[["Tv"]],
  Tw = parameters[["Tw"]],
  v = parameters[["v"]],
  w = parameters[["w"]]
)

parameters[["alpha"]] <- calculate_alpha(
  v = parameters[["v"]],
  mR = parameters[["mR"]],
  r = parameters["r"],
  c = parameters[["c"]],
  gR = parameters[["gR"]],
  g0 = parameters[["g0"]],
  Tv = parameters[["Tv"]],
  Tw = parameters[["Tw"]],
  w = parameters[["w"]]
)

# Time of scenario
time <- seq(0,1000,by = 1)

# Population-level number of susceptible
S_t <- c()
# Population-level number of infecteds
I_t <- c()
# Average transmission rate
beta_t <- c()

# Starting population parameters
I_t[1] <- 1
S_t[1] <- 25 - I_t[1]

# Number of infecteds in each cohort
I_c <- matrix(data = 0, nrow = length(time), ncol = length(time))
I_c[1,1] <- I_t[1]

# Virus and leukocytes for each cohort (time x cohort)
V <- matrix(data = 0, nrow = length(time)+1, ncol = length(time)+1)
L <- matrix(data = 0, nrow = length(time)+1, ncol = length(time)+1)


# ---- Multi-scale simulation ----

# Create a progress bar object
pb <- progress::progress_bar$new(
  format = "  Progress: [:bar] :percent  elapsed: :elapsed",
  total = length(time),
  clear = FALSE
)

# Iterate through time
for(i in 1:length(time)){
  
  # Update progress bar
  pb$tick()
  
  ### WITHIN HOST MODEL ###
  
  # Vector of cohort-time FOIs
  foi_ij <- c()
  
  # Loop over cohorts to calculate cohort-time FOIs
  for(j in 1:i){
    
    # Calculate the number of infecteds
    if(i == 1) I_c[i,j] <- I_t[1]
    if(i != 1) I_c[i,j] <- sum(I_t[1:i]) - sum(I_t[1:(i-1)])
    
    # If this is a new cohort, start V and L
    if(j == 1) V[i,j] <- 1e-07
    if(j == 1) L[i,j] <- 1
    
    # Initial state values for the within-host model
    initial_state_WH <- c(V = V[i,j],   
                          L = L[i,j])
    
    # Solve the within-host model for this cohort
    WH_model <- ode(y = initial_state_WH, times = seq(0, 1, by = 1), 
                    func = within_host_model, parms = parameters) %>%
      as.data.frame()
    
    # Update V and L for this cohort
    V[i+1,j] <- WH_model$V[2]
    L[i+1,j] <- WH_model$L[2]
    
    # Force of infection (FOI) for this cohort at this time
    foi_ij[j] <- WH_model$beta[2] * I_c[i,j]
    
  }
  
  # Average transmission rate across cohorts at this time
  beta_t[i] <- sum(foi_ij) / sum(I_c[i,])
  beta_t[i] <- 0.0005
  
  ### POPULATION MODEL ###
  
  # Assigning calculated beta as a parameter
  parameters["beta"] <- beta_t[i]
  
  # Initial state values
  initial_state_POP <- c(S = S_t[i],   # Susceptible individuals
                         I = I_t[i])     # Infected individuals
  
  # Solve the pop model
  POP_model <- ode(y = initial_state_POP, times = seq(0, 1, by = 1), 
                   func = population_model, parms = parameters) %>%
    as.data.frame()
  
  # Update the vector for susceptibles and infecteds
  S_t[i+1] <- POP_model$S[2]
  I_t[i+1] <- POP_model$I[2]
  
}


# ---- Aggregate and plot simulation ----

# Compile the data
df <- data.frame(
  time = time,
  S = S_t[1:length(time)],
  I = I_t[1:length(time)],
  β = beta_t) %>%
  mutate(N = S+I) %>%
  pivot_longer(
    cols = 2:5,
    names_to = "Variable", 
    values_to = "value") %>%
  mutate(Variable = factor(
    Variable, 
    levels = c("N","S","I","β")))

# Plot
ggplot(df, aes(x = time, y = value, color = Variable)) +
  facet_wrap(~Variable, scales = "free") +
  geom_line(linewidth = 0.5) +
  labs(x = "Time", y = NULL, color = NULL,
       title = "Multi-scale disease dynamics model") +
  theme_minimal() +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        panel.grid.minor = element_blank())







