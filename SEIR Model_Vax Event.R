# Basic SEIR Model with Vaccination Event

library(deSolve)
library(ggplot2)
library(reshape2)

seir_model_vaccination <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dE <- beta * S * I - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I
    return(list(c(dS, dE, dI, dR)))
  })
}

# Event function to handle vaccination
vaccination_event <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    if (time == vacc_time) {
      state["S"] <- state["S"] - vacc_fraction * state["S"]
      state["R"] <- state["R"] + vacc_fraction * state["S"]
    }
    return(state)
  })
}

parameters <- c(beta = 0.01, sigma = 0.05, gamma = 0.1, vacc_time = 25, vacc_fraction = 0.5)
initial_state <- c(S = 99, E = 1, I = 0, R = 0)
times <- seq(0, 100, by = 0.5)

output <- ode(y = initial_state, times = times, func = seir_model_vaccination, parms = parameters,
              events = list(func = vaccination_event, time = times))


# Visualising the output

output_df <- as.data.frame(output)

output_long <- reshape2::melt(output_df, id.vars = "time")

# Plot
ggplot(output_long, aes(x = time, y = value, color = variable)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(x = "Time", y = "Number of Individuals", color = "Compartment") +
  ggtitle("SEIR Model Dynamics with Vaccination") +
  theme(plot.title = element_text(hjust = 0.5))


# Find the peak value of the 'I' compartment
max(output_df$I)

