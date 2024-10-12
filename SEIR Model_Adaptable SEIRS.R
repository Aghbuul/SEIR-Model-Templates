# Adaptable SEIRS Model

library(deSolve)
library(ggplot2)
library(reshape2)

# SEIRS model function
seirs_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I + delta * R
    dE <- beta * S * I - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I - delta * R
    return(list(c(dS, dE, dI, dR)))
  })
}

# Parameters for the SEIRS model
parameters <- c(beta = 0.25, sigma = 0.05, gamma = 0.1, delta = 0.05)
# Initial state for the SEIRS model
initial_state <- c(S = 99, E = 1, I = 0, R = 0)
# Time sequence for the simulation
times <- seq(0, 100, by = 1)

# Solving the SEIRS model
output <- ode(y = initial_state, times = times, func = seirs_model, parms = parameters)

# Visualising the output
output_df <- as.data.frame(output)

output_long <- reshape2::melt(output_df, id.vars = "time")

# Plotting the results
ggplot(output_long, aes(x = time, y = value, color = variable)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(x = "Time", y = "Number of Individuals", color = "Compartment") +
  ggtitle("SEIRS Model Dynamics") +
  theme(plot.title = element_text(hjust = 0.5))

# Modification examples:
# SI Model: Set sigma = 0, gamma = 0, and delta = 0.
# SIS Model: Set sigma = 0, gamma = 0, and delta > 0.
# SIR Model: Set sigma = 0 and delta = 0.
# SIRS Model: Set sigma = 0.
# SEIR Model: Set delta = 0.

# Basically, E corresponds to sigma, R to gamma, and S(at the end) to delta.



