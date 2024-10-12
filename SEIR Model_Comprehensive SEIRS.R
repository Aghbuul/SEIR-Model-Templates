library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)

# Comprehensive SEIRS model function (SEIRQC)
seirs_comprehensive_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Total population
    total_pop <- S_young + E_young + I_young + R_young + Q_young + C_young +
      S_old + E_old + I_old + R_old + Q_old + C_old
    
    # Vaccination intervention for both age groups
    if (time >= vacc_start & time <= vacc_end) {
      vacc_young <- vacc_rate * S_young
      vacc_old <- vacc_rate * S_old
    } else {
      vacc_young <- 0
      vacc_old <- 0
    }
    
    # Quarantine intervention for both age groups
    if (time >= quarantine_start & time <= quarantine_end) {
      new_quarantined_young <- quarantine_rate * I_young
      new_quarantined_old <- quarantine_rate * I_old
    } else {
      new_quarantined_young <- 0
      new_quarantined_old <- 0
    }
    
    # Differential equations for young age group
    dS_young <- birth_rate * total_pop - beta_young * S_young * (I_young + I_old + kappa * C_young + kappa * C_old) - death_rate * S_young + delta * R_young - vacc_young
    dE_young <- beta_young * S_young * (I_young + I_old + kappa * C_young + kappa * C_old) - sigma * E_young - death_rate * E_young
    dI_young <- (1-p) * sigma * E_young - gamma * I_young - death_rate * I_young - new_quarantined_young
    dR_young <- gamma * I_young - delta * R_young - death_rate * R_young + vacc_young + gamma * Q_young + gamma * C_young
    dQ_young <- new_quarantined_young - gamma * Q_young - death_rate * Q_young
    dC_young <- p * sigma * E_young - gamma * C_young - death_rate * C_young
    
    # Differential equations for old age group
    dS_old <- - beta_old * S_old * (I_young + I_old + kappa * C_young + kappa * C_old) - death_rate * S_old + delta * R_old - vacc_old
    #There should be no birth rate term for the old age group!
    dE_old <- beta_old * S_old * (I_young + I_old + kappa * C_young + kappa * C_old) - sigma * E_old - death_rate * E_old
    dI_old <- (1-p) * sigma * E_old - gamma * I_old - death_rate * I_old - new_quarantined_old
    dR_old <- gamma * I_old - delta * R_old - death_rate * R_old + vacc_old + gamma * Q_old + gamma * C_old
    dQ_old <- new_quarantined_old - gamma * Q_old - death_rate * Q_old
    dC_old <- p * sigma * E_old - gamma * C_old - death_rate * C_old
    
    return(list(c(dS_young, dE_young, dI_young, dR_young, dQ_young, dC_young,
                  dS_old, dE_old, dI_old, dR_old, dQ_old, dC_old)))
  })
}

# Parameters for the comprehensive SEIRS model
parameters <- c(
  beta_young = 0.1,        # Transmission rate for young
  beta_old = 0.1,           # Transmission rate for old
  sigma = 0.05,             # Latency
  gamma = 0.01,              # Recovery rate
  delta = 0.01,             # Loss of immunity rate
  birth_rate = 0.01,        # Birth rate
  death_rate = 0.01,        # Natural death rate
  vacc_rate = 0.05,          # Vaccination rate
  vacc_start = 5,          # Vaccination intervention start time
  vacc_end = 20,            # Vaccination intervention end time
  quarantine_rate = 0.2,   # Quarantine rate
  quarantine_start = 20,    # Quarantine intervention start time
  quarantine_end = 40,      # Quarantine intervention end time
  kappa = 0.5,              # Relative transmissibility of carriers
  p = 0.3                  # Proportion progressing to carrier state
)

# Initial state for the comprehensive SEIRS model
initial_state <- c(S_young = 50, E_young = 1, I_young = 0, R_young = 0, Q_young = 0, C_young = 0,
                   S_old = 49, E_old = 0, I_old = 0, R_old = 0, Q_old = 0, C_old = 0)

# Time sequence for the simulation
times <- seq(0, 100, by = 1)

# Solving the comprehensive SEIRS model
output <- ode(y = initial_state, times = times, func = seirs_comprehensive_model, parms = parameters)



# Visualising the output
output_df <- as.data.frame(output)
output_long <- reshape2::melt(output_df, id.vars = "time")

# Plotting the results
ggplot(output_long, aes(x = time, y = value, color = variable)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(x = "Time", y = "Number of Individuals", color = "Compartment") +
  ggtitle("Comprehensive SEIRS Model Dynamics") +
  theme(plot.title = element_text(hjust = 0.5))



### Clearer Visualisation


library(deSolve)
library(ggplot2)
library(reshape2)

output_df <- as.data.frame(output)
output_long <- reshape2::melt(output_df, id.vars = "time")


## Plotting separate plots by age and compartment

# Simple code (default)
 # output_long$age_group <- ifelse(grepl("_young", output_long$variable), "Young", "Old")
 # output_long$compartment <- sub("_young|_old", "", output_long$variable)

# More complex code setting the order for groups to appear in the plot chart (eg. young before old, etc.)
output_long$age_group <- factor(ifelse(grepl("_young", output_long$variable), "Young", "Old"), levels = c("Young", "Old"))
output_long$compartment <- factor(sub("_young|_old", "", output_long$variable), levels = c("S", "C", "Q", "I", "E", "R"))


ggplot(output_long, aes(x = time, y = value, color = compartment, linetype = compartment)) +
  geom_line(size = 1) +
  facet_wrap(~age_group, scales = "free_y") +
  theme_minimal() +
  labs(x = "Time", y = "Number of Individuals", color = "Compartment", linetype = "Compartment") +
  ggtitle("Comprehensive SEIRS Model Dynamics") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

ggplot(output_long, aes(x = time, y = value, color = age_group, linetype = age_group)) +
  geom_line(size = 1) +
  facet_wrap(~compartment, scales = "free_y") +
  theme_minimal() +
  labs(x = "Time", y = "Number of Individuals", color = "Age Group", linetype = "Age Group") +
  ggtitle("Comprehensive SEIRS Model Dynamics") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")


## Quarantine Intervention Impact Analysis

# Calculate the epidemic dynamics without quarantine
parameters_no_quarantine <- parameters
parameters_no_quarantine["quarantine_rate"] = 0
output_no_quarantine <- ode(y = initial_state, times = times, func = seirs_comprehensive_model, parms = parameters_no_quarantine)

# Calculate the epidemic dynamics with quarantine
output_quarantine <- ode(y = initial_state, times = times, func = seirs_comprehensive_model, parms = parameters)

# Convert to data frames
output_no_quarantine_df <- as.data.frame(output_no_quarantine)
output_quarantine_df <- as.data.frame(output_quarantine)

# Calculate cumulative infections
output_no_quarantine_df$Cumulative_Infected_no_quarantine <- cumsum(rowSums(output_no_quarantine_df[, c("I_young", "I_old")]))
output_quarantine_df$Cumulative_Infected_quarantine <- cumsum(rowSums(output_quarantine_df[, c("I_young", "I_old")]))

# Calculate number of infected individuals at each time point
output_no_quarantine_df$Total_Infected_no_quarantine <- rowSums(output_no_quarantine_df[, c("I_young", "I_old")])
output_quarantine_df$Total_Infected_quarantine <- rowSums(output_quarantine_df[, c("I_young", "I_old")])

# Plot number of infected individuals at each time point for both scenarios
df_infected_quarantine <- data.frame(
  time = rep(times, 2),
  Total_Infected = c(output_no_quarantine_df$Total_Infected_no_quarantine, output_quarantine_df$Total_Infected_quarantine),
  Scenario = rep(c("No Quarantine", "Quarantine"), each = length(times))
)

ggplot(df_infected_quarantine, aes(x = time, y = Total_Infected, color = Scenario)) +
  geom_line(size = 1.5) +
  theme_minimal() +
  labs(x = "Time", y = "Number of Infected Individuals", title = "Number of Infected Individuals Over Time") +
  theme(plot.title = element_text(hjust = 0.5))

# Calculate percentage reduction in cumulative infections
percentage_prevented_quarantine <- (output_no_quarantine_df$Cumulative_Infected_no_quarantine - output_quarantine_df$Cumulative_Infected_quarantine) / output_no_quarantine_df$Cumulative_Infected_no_quarantine * 100

# Plot percentage reduction in cumulative infections by quarantine
df_prevented_quarantine <- data.frame(time = times, percentage_prevented = percentage_prevented_quarantine)
ggplot(df_prevented_quarantine, aes(x = time, y = percentage_prevented)) +
  geom_line(color = "red", size = 1.5) +
  theme_minimal() +
  labs(x = "Time", y = "Percentage Reduction in Infections", title = "Percentage Reduction in Cumulative Infections by Quarantine") +
  theme(plot.title = element_text(hjust = 0.5))


## Vaccine Intervention Impact Analysis

# Calculate the epidemic dynamics without vaccination
parameters_no_vacc <- parameters
parameters_no_vacc["vacc_rate"] = 0
output_no_vacc <- ode(y = initial_state, times = times, func = seirs_comprehensive_model, parms = parameters_no_vacc)

# Calculate the epidemic dynamics with vaccination
output_vacc <- ode(y = initial_state, times = times, func = seirs_comprehensive_model, parms = parameters)

# Convert to data frames
output_no_vacc_df <- as.data.frame(output_no_vacc)
output_vacc_df <- as.data.frame(output_vacc)

# Calculate number of infected individuals at each time point
output_no_vacc_df$Total_Infected_no_vacc <- rowSums(output_no_vacc_df[, c("I_young", "I_old")])
output_vacc_df$Total_Infected_vacc <- rowSums(output_vacc_df[, c("I_young", "I_old")])

# Combine data frames for plotting
df_infected_vacc <- data.frame(
  time = rep(times, 2),
  Total_Infected = c(output_no_vacc_df$Total_Infected_no_vacc, output_vacc_df$Total_Infected_vacc),
  Scenario = rep(c("No Vaccination", "Vaccination"), each = length(times))
)

# Plot number of infected individuals at each time point for both scenarios
ggplot(df_infected_vacc, aes(x = time, y = Total_Infected, color = Scenario)) +
  geom_line(size = 1.5) +
  theme_minimal() +
  labs(x = "Time", y = "Number of Infected Individuals", title = "Number of Infected Individuals Over Time") +
  theme(plot.title = element_text(hjust = 0.5))

# Calculate cumulative infections
output_no_vacc_df$Cumulative_Infected_no_vacc <- cumsum(rowSums(output_no_vacc_df[, c("I_young", "I_old")]))
output_vacc_df$Cumulative_Infected_vacc <- cumsum(rowSums(output_vacc_df[, c("I_young", "I_old")]))

# Calculate percentage reduction in cumulative infections
percentage_prevented_vacc <- (output_no_vacc_df$Cumulative_Infected_no_vacc - output_vacc_df$Cumulative_Infected_vacc) / output_no_vacc_df$Cumulative_Infected_no_vacc * 100

# Plot percentage reduction in cumulative infections by vaccination
df_prevented_vacc <- data.frame(time = times, percentage_prevented = percentage_prevented_vacc)
ggplot(df_prevented_vacc, aes(x = time, y = percentage_prevented)) +
  geom_line(color = "blue", size = 1.5) +
  theme_minimal() +
  labs(x = "Time", y = "Percentage Reduction in Infections", title = "Percentage Reduction in Cumulative Infections by Vaccination") +
  theme(plot.title = element_text(hjust = 0.5))


