
library(dplyr)
library(ggplot2)


# Example species parameter table
species_params <- data.frame(
  Species = c("Mouse", "Rat", "Dog"),
  BW = c(0.025, 0.25, 10)
)


model_cl_fit <- function(df_data) {
  validate(
    need(length(unique(df_data$Species[!is.na(df_data$CL)])) >= 2, 
         "Please enter CL data for at least two species.")
  )
  model <- lm(logCL ~ logBW, data = df_data)
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  CI_95 <- confint(model, 2, level = 0.95)
  R2 <- summary(model)$r.squared
  log10_clin_CL <- slope * log10(55) + intercept
  clin_CL <- 10^log10_clin_CL
  
  list(model = model, intercept = intercept, slope = slope, CI_95 = CI_95, R2 = R2, clin_CL = clin_CL)
}


model_v_fit <- function(df_data) {
  validate(
    need(length(unique(df_data$Species[!is.na(df_data$V)])) >= 2, 
         "Please enter V data for at least two species.")
  )
  
  model <- lm(logV ~ logBW, data = df_data)
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  CI_95 <- confint(model, 2, level = 0.95)
  R2 <- summary(model)$r.squared
  log10_clin_V <- slope * log10(55) + intercept
  clin_V <- 10^log10_clin_V
  
  list(model = model, intercept = intercept, slope = slope, CI_95 = CI_95, R2 = R2, clin_V = clin_V)
}

# Plot
plot_log_fit <- function(df_data, model_res, yvar, ylab_txt) {
  ggplot(df_data, aes(x = logBW, y = .data[[yvar]], color = Species)) +
    geom_abline(slope = model_res$slope, intercept = model_res$intercept, color = "darkgrey", size = 0.6) +
    geom_point(size = 3) +
    xlab("log10 Body Weight (kg)") +
    ylab(ylab_txt) +
    scale_color_manual(values = c("red", "blue", "orange")) +
    theme_bw() +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}