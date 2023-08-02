
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

library(cowplot) #for manuscript ready figures
library(lme4) #for lmer & glmer models
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(effects)
library(sjstats) #use for r2 functions

getwd()

inPATH = '../data/training/iapar/model_selection_data.csv'

df = read.csv(inPATH)
str(df)
df$Year <- as.factor(df$Year)
df$Treatment <- as.factor(df$Treatment)
df$Block <- as.factor(df$Block)
df$Ecosite <- as.factor(df$Ecosite)
df$Graze_timing <- as.factor(df$Graze_timing)

response_var <-  'Total_Biomass'

#global model 1 used in dredge
global.model.graze<-lmer(paste0(response_var, " ~ iAPAR + Graze_timing + Ecosite + Graze_timing:iAPAR + Ecosite:iAPAR + (1|Year) + (1|Year:Block)"),
                                na.action = "na.fail", data=df)

global.model.graze

#dredge model
dr.graze<-dredge(global.model.graze)
dr.graze
#get top model
top_model.graze <- get.models(dr.graze, 1)[[1]]
model_summary.graze <- summary(get.models(dr.graze, 1)[[1]])
plot(top_model.graze)

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


#global model Treatment used in dredge
global.model.trt<-lmer(paste0(response_var, "~ iAPAR + Treatment + Ecosite + Treatment:iAPAR + Ecosite:iAPAR + (1|Year) + (1|Year:Block)"), na.action = "na.fail", data=df)

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
