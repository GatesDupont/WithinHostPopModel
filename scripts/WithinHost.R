library(deSolve)
library(ggplot2)

calculate_r_star <- function(c, gR, g0, mR, mu, Tv, Tw, v, w) {
  
  # Calculate the terms inside the square root
  sqrt_term <- sqrt(c^2 * mR^2 * gR * g0 * mu * Tv * Tw * (v*Tw + gR * w*Tv))
  
  # Calculate r_star according to the given formula
  r_star <- (c * gR / mR) + (sqrt_term / (mR^2 * (Tw + gR * w * Tv)))
  
  return(r_star)
}

# Define the within-host model
within_host_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dV <- r * V - c * V * L
    dL <- g0 + gR * r * V - mR * L
    beta <- zeta * V
    return(list(c(dV, dL), beta = beta))
  })
}

# Initial state
initial_state <- c(V = 0.00001, L = 0.1)  # Example initial values

# Time points
times <- seq(0, 30, by = 1)

# Define parameters
parameters <- c(
  zeta = 0.2, 
  c = 0.5, 
  g0 = 0.3, 
  gR = 0.9, 
  mR = 1/21)

parameters["r"] <- calculate_r_star(
  c = parameters[["c"]],
  gR = parameters[["gR"]],
  g0 = parameters[["g0"]],
  mR = parameters[["mR"]],
  mu = 1 / (20*365),
  Tv = 0.5,
  Tw = 0.5,
  v = 1,
  w = 1
)

# Solve the differential equations
output <- ode(y = initial_state, times = times, func = within_host_model, parms = parameters)

# Convert output to a data frame
output_df <- as.data.frame(output) %>%
  pivot_longer(cols = 2:4) %>%
  mutate(name = factor(name, levels = c("V","L","beta")))

# Plotting
ggplot(data = output_df, aes(x = time, y = value, color = name)) +
  facet_wrap(~ name, scales = "free") +
  geom_line() +
  labs(title = "Dynamics of Virus, Leukocyte, and Transmission Rate",
       color = NULL,
       x = "Time (days)",
       y = "Population/Rate") +
  theme_minimal()

