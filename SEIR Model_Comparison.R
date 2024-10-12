# SEIR Model Comparison

library(deSolve)
library(ggplot2)
library(reshape2)
library(gridExtra)

seir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dE <- beta * S * I - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I
    return(list(c(dS, dE, dI, dR)))
  })
}

# List of parameter sets
parameters_list <- list(
  list(beta = 0.05, sigma = 0.1, gamma = 0.10),
  list(beta = 0.05, sigma = 0.1, gamma = 0.05),
  list(beta = 0.20, sigma = 0.1, gamma = 0.10),
  list(beta = 0.20, sigma = 0.1, gamma = 0.05)
)

initial_state <- c(S = 99, E = 1, I = 0, R = 0)
times <- seq(0, 100, by = 1)

# Create and store plots
plots <- list()

for (i in 1:length(parameters_list)) {
  params <- parameters_list[[i]]
  parameters <- c(beta = params$beta, sigma = params$sigma, gamma = params$gamma)
  
  output <- ode(y = initial_state, times = times, func = seir_model, parms 
                = parameters)
  
  
  # Visualising the output
  
  output_df <- as.data.frame(output)
  
  output_long <- melt(output_df, id.vars = "time")
  
  # Plot
  plot <- ggplot(output_long, aes(x = time, y = value, color = variable)) +
    geom_line(size = 1.5) +
    theme_minimal() +
    labs(x = "Time", y = "Number of Individuals", color = "Compartment") +
    ggtitle(paste("SEIR Model Dynamics\n(beta =", params$beta, ", sigma =", params$sigma, ", gamma =", params$gamma, ")")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  plots[[i]] <- plot
}

# Arrange plots on a single page
grid.arrange(
  grobs = plots,
  ncol = 2,
  nrow = 2
)

# Generate and examine each plot separately
plot1 <- plots[[1]]
plot2 <- plots[[2]]
plot3 <- plots[[3]]
plot4 <- plots[[4]]

print(plot1)
print(plot2)
print(plot3)
print(plot4)
