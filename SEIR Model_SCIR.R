# Basic SCIR Model (C = Carrier)
# Note: This is almost identical to the Basic SEIR Model code.

library(deSolve)
library(ggplot2)
library(reshape2)

scir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * (I + C)
    dC <- beta * S * (I + C) - sigma * C
    dI <- sigma * C - gamma * I
    dR <- gamma * I
    return(list(c(dS, dC, dI, dR)))
  })
}

parameters <- c(beta = 0.25, sigma = 0.05, gamma = 0.1)
initial_state <- c(S = 99, C = 1, I = 0, R = 0)
times <- seq(0, 100, by = 1)

output <- ode(y = initial_state, times = times, func = scir_model, parms 
              = parameters)

# Visualising the output

output_df <- as.data.frame(output)

output_long <- reshape2::melt(output_df, id.vars = "time")

# Plot
ggplot(output_long, aes(x = time, y = value, color = variable)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(x = "Time", y = "Number of Individuals", color = "Compartment") +
  ggtitle("SCIR Model Dynamics") +
  theme(plot.title = element_text(hjust = 0.5))
