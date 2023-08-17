
# load relevant libraries
library(lme4) 
library(nlme)
library(MuMIn)
library(car) 
library(tidyverse) #for all data wrangling
library(dplyr)
library(lubridate)
library(ggplot2)
library(lmerTest)
library(sjPlot)
library(ggeffects)
library(ggpmisc)
library(cowplot) #for manuscript ready figures
library(lme4) #for lmer & glmer models
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(effects)
library(sjstats) #use for r2 functions
library(broom.mixed)
getwd()

inPATH = '../data/training/iapar/model_selection_data.csv'

df = read.csv(inPATH)
str(df)
df$Year <- as.factor(df$Year)
df$Treatment <- as.factor(df$Treatment)
df$Block <- as.factor(df$Block)
df$Ecosite <- as.factor(df$Ecosite)
df$Graze_timing <- as.factor(df$Graze_timing)

response_var <-  'Total_Biomass_sqrt'

#global model graze used in dredge ----
# global.model.graze<-lmer(paste0(response_var, " ~ iAPAR + Graze_timing + Ecosite + Graze_timing:iAPAR + Ecosite:iAPAR + (1|Year) + (1|Year:Block)"),
#                                 na.action = "na.fail", data=df)
global.model.graze<-lmer(paste0(response_var, " ~ iAPAR + Graze_timing:iAPAR + Ecosite:iAPAR + (1|Year) + (1|Year:Block)"),
                         na.action = "na.fail", data=df)

global.model.graze

#dredge model
dr.graze<-dredge(global.model.graze)
dr.graze
#get top model
top_model.graze <- get.models(dr.graze, 1)[[1]]
model_summary.graze <- summary(get.models(dr.graze, 1)[[1]])
plot(top_model.graze)
summary(model_summary.graze)

# Calculate Cook's distance and leverage
cooksd <- cooks.distance(top_model.graze)
leverage <- hatvalues(top_model.graze)

# Identify outliers based on Cook's distance and leverage
outliers <- which(cooksd > 1 | leverage > 2)

# Check for normality of residuals
hist(residuals(top_model.graze))
qqnorm(residuals(top_model.graze))
shapiro.test(residuals(top_model.graze))

# Check for multicollinearity (collinear if VIF > 5-10)
vif(top_model.graze)


# !!!!r2 does not take into account random effects!!!!
r.squaredGLMM(top_model.graze)

r.squaredGLMM(top_model.graze, partial =  "Ecosite")

tidy(top_model.graze)

# Extract the effect sizes from the model summary
effect_sizes.graze <- model_summary.graze$coefficients[, "Estimate"]

# Create a data frame with the variable names and effect sizes
effect_sizes_df.graze <- data.frame(Variable = rownames(model_summary.graze$coefficients),
                                  EffectSize = effect_sizes.graze)

# Create the plot using ggplot2
library(ggplot2)

ggplot(effect_sizes_df.graze, aes(x = Variable, y = EffectSize)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Effect Sizes",
       x = "Variable",
       y = "Effect Size")

# Extract the model coefficients and their standard errors
coefficients.graze <- as.data.frame(coef(summary(model_summary.graze)))
coefficients.graze$Variables <- rownames(coefficients.graze)

# Remove the "(Intercept)" row
coefficients.graze <- coefficients.graze[-1, ]

coefficients.graze$Std.Error <- coefficients.graze$`Std. Error`

str(coefficients.graze)
ggplot(coefficients.graze, aes(x = Variables, y = Estimate, ymin = Estimate - (2 * Std.Error), ymax = Estimate + (2 * Std.Error))) +
  geom_linerange(color = "steelblue") +
  geom_point(color = "steelblue", size = 3) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Model Coefficients",
       x = "Variables",
       y = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Get predicted values for each ecosite
df$predict_sqrt <- predict(top_model.graze)

df$predict_transformed <- df$predict_sqrt ** 2

gg_ecosite_graze_yr <- ggplot() +
  geom_point(data = df, aes(x = iAPAR, y = Total_Biomass, color = Ecosite)) +
  geom_smooth(data = df, aes(x = iAPAR, y = predict_transformed , color = Ecosite), size = 1,method = "lm", se = TRUE) +
  facet_wrap(~ Year) +
  # stat_poly_line() +
  # stat_poly_eq(data = df,aes(x = iAPAR, y = predict, color = Ecosite))+
  theme_bw()
  
# geom_ribbon(data = ggpred.graze, aes(x = x, ymin = (predicted - std.error), ymax = (predicted + std.error), fill = group), alpha = 0.2) 
gg_ecosite_graze_yr
gg_ecosite_graze <- ggplot() +
  geom_point(data = df, aes(x = iAPAR, y = Total_Biomass, color = Ecosite)) +
  geom_smooth(data = df, aes(x = iAPAR, y = predict_transformed , color = Ecosite), size = 1,method = "lm", se = TRUE) +
  # stat_poly_line() +
  # stat_poly_eq(data = df,aes(x = iAPAR, y = predict, color = Ecosite))+
  theme_bw()

# geom_ribbon(data = ggpred.graze, aes(x = x, ymin = (predicted - std.error), ymax = (predicted + std.error), fill = group), alpha = 0.2) 
gg_ecosite_graze


gg_graze <- ggplot() +
  geom_point(data = df, aes(x = iAPAR, y = Total_Biomass, color = Graze_timing)) +
  geom_smooth(data = df, aes(x = iAPAR, y = predict_transformed , color = Graze_timing), size = 1,method = "lm", se = TRUE) +
  # stat_poly_line() +
  # stat_poly_eq(data = df,aes(x = iAPAR, y = predict, color = Ecosite))+
  theme_bw()

# geom_ribbon(data = ggpred.graze, aes(x = x, ymin = (predicted - std.error), ymax = (predicted + std.error), fill = group), alpha = 0.2) 
gg_graze

#global model Treatment used in dredge ----
# global.model.trt<-lmer(paste0(response_var, "~ iAPAR + Treatment + Ecosite + Treatment:iAPAR + Ecosite:iAPAR + (1|Year) + (1|Year:Block)"), na.action = "na.fail", data=df)
global.model.trt<-lmer(paste0(response_var, "~ iAPAR + Treatment:iAPAR + Ecosite:iAPAR + (1|Year) + (1|Year:Block)"), na.action = "na.fail", data=df)

global.model.trt

#dredge model
dr.trt<-dredge(global.model.trt)
dr.trt
#get top model
top_model.trt <- get.models(dr.trt, 1)[[1]]
model_summary.trt <- summary(get.models(dr.trt, 1)[[1]])
model_summary.trt
plot(top_model.trt)

# Extract the effect sizes from the model summary
effect_sizes.trt <- model_summary.trt$coefficients[, "Estimate"]

# Create a data frame with the variable names and effect sizes
effect_sizes_df.trt <- data.frame(Variable = rownames(model_summary.trt$coefficients),
                              EffectSize = effect_sizes.trt)


ggplot(effect_sizes_df.trt, aes(x = Variable, y = EffectSize)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Effect Sizes",
       x = "Variable",
       y = "Effect Size")

# Extract the model coefficients and their standard errors
coefficients.trt <- as.data.frame(coef(summary(model_summary.trt)))
coefficients.trt$Variables <- rownames(coefficients.trt)

# Remove the "(Intercept)" row
coefficients.trt <- coefficients.trt[-1, ]


coefficients.trt$Std.Error <- coefficients.trt$`Std. Error`

str(coefficients.trt)
ggplot(coefficients.trt, aes(x = Variables, y = Estimate, ymin = Estimate - (2 * Std.Error), ymax = Estimate + (2 * Std.Error))) +
  geom_linerange(color = "steelblue") +
  geom_point(color = "steelblue", size = 3) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Model Coefficients",
       x = "Variables",
       y = "Estimate") 


# Functional groups model ----
str(df)

fg.model<-lmer(paste0("sqrt(iAPAR) ~ AG + FORB +BOBU + C3PG + SD + WSPG + (1|Plot)"), na.action = "na.fail", data=df)
fg.model
plot(fg.model)

df$fg_pred_values_sqrt <- predict(fg.model)
standard_errors <- summary(fg.model)$coefficients


gg_fg <- ggplot() +
  geom_point(data = df, aes(x = iAPAR, y = fg_pred_values_sqrt)) +
   # stat_poly_line() +
  # stat_poly_eq(data = df,aes(x = iAPAR, y = predict, color = Ecosite))+
  theme_bw()

# geom_ribbon(data = ggpred.graze, aes(x = x, ymin = (predicted - std.error), ymax = (predicted + std.error), fill = group), alpha = 0.2) 
gg_fg

# Temporal model ----
# Create a data frame to store the results
# Create a data frame to store the results
results <- data.frame()

# Loop over all years
for (year in 2014:2022) {
  # Subset the data by year
  data_year <- df %>% filter(Year == year)
  
  # Run the model for the year
  model_year <- lmer(Total_Biomass_sqrt ~ iAPAR + (1|Block) +
                       iAPAR:Ecosite, data = data_year)
  
  # Get the model summary
  summary_year <- summary(model_year)
  
  # Extract the R2 value, slope, and intercept
  r2_year <- r.squaredGLMM(model_year)[1]
  slope_year <- summary_year$coefficients[2, 1]
  intercept_year <- summary_year$coefficients[1, 1]
  
  # Add the results to the data frame
  results <- bind_rows(results,
                       tibble(Year = year, AIC = AIC(model_year),
                              R2 = r2_year, Slope = slope_year,
                              Intercept = intercept_year))
}

# Name the columns of the data frame
colnames(results) <- c("Year", "AIC", "R2", "Slope", "Intercept")

# Print the data frame
results
