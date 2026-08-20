#ST563 FINAL Project
#Jerry, Younha, Silvia

# Import the dataset
data<-read.csv("C:\\Users\\younh\\Documents\\who.csv")

############################
# Data Preprocessing
# Remove the columns Total.expenditure and Alcohol due to large amount of missing values
data <- data[, !colnames(data) %in% c("Total.expenditure", "Alcohol")]
summary(data)

# Remove rows with any missing data
data <- na.omit(data)

#Convert categorical Status variable to numeric
data$Status <- ifelse(data$Status == "Developed", 1, 0)
##############################

########################################################
##### Check Assumptions
# Load necessary libraries
library(ggplot2)
library(car)
library(corrplot)

# 1. Check for normality
# Fit a model to check residual assumptions
lm_full <- lm(Life.expectancy ~ GDP + Schooling + HIV.AIDS + Status + Polio + BMI, data = data)
residuals <- lm_full$residuals

# Q-Q plot
qqnorm(residuals)
qqline(residuals, col = "red")

# 2. Check for homoscedasticity
plot(lm_full$fitted.values, residuals, 
     xlab = "Fitted Values", ylab = "Residuals", 
     main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")

# 3. Check for multicollinearity
# Load the car package
library(car)

# Fit a linear regression model
lm_model <- lm(Life.expectancy ~ GDP + Schooling + HIV.AIDS + Status + Polio + BMI, data = data)

# Calculate VIF values
vif_values <- vif(lm_model)

# Print the VIF values
print("VIF values for predictors:")
print(vif_values)

# 4. Check for influential points 
# Fit the linear regression model
lm_model <- lm(Life.expectancy ~ GDP + Schooling + HIV.AIDS + Status + Polio + BMI, data = data)

# Calculate Cook's Distance
cooksd <- cooks.distance(lm_model)

# Plot Cook's Distance
plot(cooksd, main = "Cook's Distance", ylab = "Cook's Distance", xlab = "Observation Index", type = "h")
abline(h = 0.02, col = "red", lty = 2)  # Threshold line for 0.02

# Identify influential points where Cook's Distance > 0.02
influential_points <- which(cooksd > 0.02)

# Print the indices and Cook's Distance values for these points
cat("Observations with Cook's Distance > 0.02:\n")
influential_data <- data.frame(Index = influential_points, Cooks_Distance = cooksd[influential_points])
print(influential_data)

# Optionally, extract the full rows of data for these influential points
influential_observations <- data[influential_points, ]
print("Data for Influential Observations:")
print(influential_observations)
########################################################

########################################################
##### Best Subset Selection
# Choose Best Subset Selection
library(leaps)
# Best Subset model for each model size
bestsub <- regsubsets(Life.expectancy ~ .-Country,
                      data = data,
                      nvmax = 19)
# summary
mod_summary_bestsub <- summary(bestsub)
mod_summary_bestsub

# choose the best model
metrics_bestsub <- data.frame(aic = mod_summary_bestsub$cp,
                              bic = mod_summary_bestsub$bic,
                              adjR2 = mod_summary_bestsub$adjr2)
metrics_bestsub

# Find best model
round( coef(bestsub, 13), 2)
best_model_bestsub <- lm(Life.expectancy ~Year+Status+Adult.Mortality+infant.deaths
                         +BMI+under.five.deaths+Polio+Diphtheria+HIV.AIDS+GDP
                         +thinness.5.9.years+Income.composition.of.resources+Schooling, data=data)
summary(best_model_bestsub)

# We found 'Status' and 'Polio' are not significant. So we remove them. Then we get 11 variables.
round( coef(bestsub, 11), 2)
finalbest_model_bestsub <- lm(Life.expectancy ~Year+Adult.Mortality+infant.deaths
                              +BMI+under.five.deaths+Diphtheria+HIV.AIDS+GDP
                              +thinness.5.9.years+Income.composition.of.resources+Schooling, data=data)
summary(finalbest_model_bestsub)

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Add a Model index to metrics_bestsub for plotting
metrics_bestsub$Model <- seq_len(nrow(metrics_bestsub))

metrics_bestsub_long <- metrics_bestsub %>%
  pivot_longer(cols = c(aic, bic, adjR2), names_to = "Metric", values_to = "Value")

# Plot the metrics_bestsub separately by Metric
ggplot(metrics_bestsub_long, aes(x = Model, y = Value)) +
  geom_line(size = 1) +                  # Line for each metric
  geom_point(size = 2) +                # Points for emphasis
  labs(
    title = "Model Selection metrics_bestsub (Best Subset Selection)",
    x = "Number of Predictors",
    y = "Metric Value"
  ) +
  theme_minimal() +                     # Clean theme
  scale_x_continuous(breaks = seq_len(nrow(metrics_bestsub))) +  # Show all model indices on x-axis
  facet_wrap(~ Metric, scales = "free_y")  # Separate plots for each metric
########################################################

########################################################
##### Backward Stepwise Selection
# Choose backward stepwise selection
library(leaps)
backward <- regsubsets(Life.expectancy ~ .-Country,
                       data = data,
                       nvmax = 19,
                       method = "backward")
# summary
mod_summary_backward <- summary(backward)
mod_summary_backward

# choose the best model for backward selection
metrics_backward <- data.frame(aic = mod_summary_backward$cp,
                               bic = mod_summary_backward$bic,
                               adjR2 = mod_summary_backward$adjr2)
metrics_backward

#Find best model for backward selection
round( coef(backward, 13), 2)
best_model_backward <- lm(Life.expectancy ~Year+Status+Adult.Mortality+infant.deaths
                          +BMI+under.five.deaths+Polio+Diphtheria+HIV.AIDS+GDP
                          +thinness.5.9.years+Income.composition.of.resources+Schooling, data=data)
summary(best_model_backward)

#We found 'Status' and 'Polio' are not significant. So we remove them. Then we get 11 variables.
round( coef(backward, 11), 2)
finalbest_model_backward <- lm(Life.expectancy ~Year+Adult.Mortality+infant.deaths
                               +BMI+under.five.deaths+Diphtheria+HIV.AIDS+GDP
                               +thinness.5.9.years+Income.composition.of.resources+Schooling, data=data)
summary(finalbest_model_backward)

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Add a Model index to metrics_backward for plotting
metrics_backward$Model <- seq_len(nrow(metrics_backward))

metrics_backward_long <- metrics_backward %>%
  pivot_longer(cols = c(aic, bic, adjR2), names_to = "Metric", values_to = "Value")

# Plot the metrics_backward separately by Metric
ggplot(metrics_backward_long, aes(x = Model, y = Value)) +
  geom_line(size = 1) +                  # Line for each metric
  geom_point(size = 2) +                # Points for emphasis
  labs(
    title = "Model Selection metrics_backward (Backward Selection)",
    x = "Number of Predictors",
    y = "Metric Value"
  ) +
  theme_minimal() +                     # Clean theme
  scale_x_continuous(breaks = seq_len(nrow(metrics_backward))) +  # Show all model indices on x-axis
  facet_wrap(~ Metric, scales = "free_y")  # Separate plots for each metric
########################################################