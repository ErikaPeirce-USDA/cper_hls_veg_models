# load relevant libraries
library(lme4) 
library(nlme)
library(MuMIn)
library(car) 
library(tidyverse) #for all data wrangling
library(dplyr)
library(lubridate)
library(reshape2)
library(ggplot2)
library(lmerTest)
library(sjPlot)
library(emmeans)
library(ggeffects)
library(ggpmisc)
library(cowplot) #for manuscript ready figures
library(lme4) #for lmer & glmer models
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(effects)
library(sjstats) #use for r2 functions
# library(broom.mixed)
# library(mgcv)
library(scales)
library(stargazer)
library(forcats)
library(explore)
library(gridExtra)
library(ggpubr)

mytheme<-function ()
{
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_line(colour = "#f0f0f0"),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.8),# get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.key = element_rect(fill = "transparent", colour = NA), # get rid of key legend fill, and of the surrounding
    axis.line = element_line(colour = "black"), # adding a black line for x and y axis
    legend.position=("right")#where legend is located 
  )}

getwd()
#Load RAP
rap_df <- read.csv('../data/training/iapar/RAP_herb.csv')
str(rap_df)
# rap_df$time <- as.Date(rap_df$time)
rap_df$Plot <- as.factor(rap_df$Plot)
rap_df$Id <- rap_df$Plot
rap_df$Year <- lubridate::year(rap_df$time)
rap_df$Year <- as.factor(rap_df$Year)
#Load model training data
inPATH = '../data/training/iapar/model_selection/cper_biomass_iapar_2014_2022_all.csv'

raw_df = read.csv(inPATH)

# All treatments
df_all <- raw_df

# Only CARM treatments
df <- raw_df %>% 
  filter(Treatment != "NEW_HVY")
  
sum <- df %>%
  group_by(Graze_timing) %>%
  summarise(sum_n = length(Total_Biomass))
sum


str(df)

df$Year <- as.factor(df$Year)
df$Treatment <- as.factor(df$Treatment)
df$Block <- as.factor(df$Block)
df$Ecosite <- as.factor(df$Ecosite)
df$Graze_timing <- as.factor(df$Graze_timing)
df$Graze_timing  <-  factor(df$Graze_timing, levels = c("Season-long","Ungrazed","Pulse"))
df$Ecosite  <-  factor(df$Ecosite, levels = c("Loamy","Sandy","Salt Flats"))

str(df)
response_var <-  'Total_Biomass'
df <- right_join(df, rap_df, by = c("Id","Year"))
df <- filter(df, !is.na(Treatment))
  


#Using a lmer ----
global.model.graze<-lmer(paste0(response_var, " ~ iAPAR + Graze_timing:iAPAR + Ecosite:iAPAR + Graze_timing:Ecosite:iAPAR + (1|Year) + (1|Id)"),
                         na.action = "na.fail", data=df)

Anova(global.model.graze, type = 3)
summary(global.model.graze)
plot(global.model.graze)

emm_graze <- emmeans(global.model.graze, pairwise ~ Graze_timing)
emm_graze

emm_eco<- emmeans(global.model.graze, pairwise ~ Ecosite)
emm_eco

#Emmeans for the interaction between graze and ecosite
emm_graze_eco <- emmeans(global.model.graze, pairwise ~ Graze_timing*Ecosite)
emm_graze_eco

multcomp::cld(emm_graze_eco, alpha = 0.05, Letters = letters)

pwpp_ecograze <- pwpp(emm_graze_eco)
pwpp_ecograze


class(global.model.graze) <- "lmerMod"
stargazer(global.model.graze, type = "html",out="test.html",
          digits = 3,
          star.cutoffs = c(0.05,0.01,0.001),
          digit.separator = "")
# how is this calculating p-values

# Calculate Cook's distance and leverage
cooksd <- cooks.distance(global.model.graze)
leverage <- hatvalues(global.model.graze)
# par(mfrow=c(2,2))
# hist(df$Total_Biomass)
# hist(df$Total_Biomass_sqrt)

# Identify outliers based on Cook's distance and leverage
outliers <- which(cooksd > 1 | leverage > 2)
outliers
# Check for normality of residuals
# hist(residuals(global.model.graze))
# qqnorm(residuals(global.model.graze))
# shapiro.test(residuals(global.model.graze))

# Check for multicollinearity (collinear if VIF > 5-10)
vif(global.model.graze)


# !!!!r2 does not take into account random effects!!!!
r.squaredGLMM(global.model.graze)

model_summary.graze <- summary(global.model.graze)
coefficients(global.model.graze)
# Extract the effect sizes from the model summary
effect_sizes.graze <- model_summary.graze$coefficients[, "Estimate"]

# Create a data frame with the variable names and effect sizes
effect_sizes_df.graze <- data.frame(Variable = rownames(model_summary.graze$coefficients),
                                    EffectSize = effect_sizes.graze)


ggplot(effect_sizes_df.graze, aes(x = Variable, y = EffectSize)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Effect Sizes",
       x = "Variable",
       y = "Effect Size")

# Extract the model coefficients and their standard errors
summary <- summary(model_summary.graze)
coefficients.graze <- as.data.frame(coef(summary(model_summary.graze)))
coefficients.graze$Variables <- rownames(coefficients.graze)

# Remove the "(Intercept)" row
# Calculate the confidence intervals for the regression coefficients
# confint(global.model.graze, level = 0.95)
coefficients.graze <- coefficients.graze[-1, ]

coefficients.graze <- rename(coefficients.graze, Std.Error = `Std. Error`, t.value = `t value`)


newdata1 <- data.frame(Estimate = 0, Std.Error =0, t.value = 0, Variables = "iAPAR:Season-long Grazed")
newdata2 <- data.frame(Estimate = 0, Std.Error =0, t.value = 0, Variables = "iAPAR:Loamy Ecosite")
newdata3 <- data.frame(Estimate = 0, Std.Error =0, t.value = 0,
                       Variables = c("iAPAR:Season-long:Loamy","iAPAR:Season-long:Sandy","iAPAR:Season-long:Salt Flats"))
newdata4 <- data.frame(Estimate = 0, Std.Error =0, t.value = 0, Variables = "iAPAR:Pulse:Loamy")
newdata5 <- data.frame(Estimate = 0, Std.Error =0, t.value = 0, Variables = "iAPAR:Ungrazed:Loamy")

coefficients.graze <- rbind(coefficients.graze,newdata1)
coefficients.graze2 <- rbind(coefficients.graze,newdata2)
coefficients.graze2 <- rbind(coefficients.graze2,newdata3)
coefficients.graze2 <- rbind(coefficients.graze2,newdata4)
coefficients.graze2 <- rbind(coefficients.graze2,newdata5)
coefficients.graze2$Variables<- as.factor(coefficients.graze2$Variables)

coefficients.graze2 <- coefficients.graze2 %>% 
  mutate(Variables=recode(Variables, "iAPAR" = "iAPAR",
  "iAPAR:Graze_timingUngrazed" = "iAPAR:Ungrazed",
  "iAPAR:Graze_timingPulse" = "iAPAR:Pulse Grazed",
  "iAPAR:Graze_timingSeasonlong" = "iAPAR:Season-long Grazed",
  "iAPAR:EcositeLoamy" = "iAPAR:Loamy Ecosite",
  "iAPAR:EcositeSalt Flats" = "iAPAR:Salt Flats Ecosite",
  "iAPAR:EcositeSandy" = "iAPAR:Sandy Ecosite",
  "iAPAR:Graze_timingUngrazed:EcositeSandy" = "iAPAR:Ungrazed:Sandy",
  "iAPAR:Graze_timingPulse:EcositeSandy" = "iAPAR:Pulse:Sandy",
  "iAPAR:Graze_timingPulse:EcositeSalt Flats" = "iAPAR:Pulse:Salt Flats",
  "iAPAR:Graze_timingUngrazed:EcositeSalt Flats" = "iAPAR:Ungrazed:Salt Flats",
  )) #%>%
  # mutate(label = c("iAPAR", "Grazing", "Grazing", "Grazing", "Ecosite", "Ecosite", "Ecosite"))

str(coefficients.graze2)
coefficients.graze2 <- coefficients.graze2 %>%
  mutate(Variables = fct_relevel(Variables, c("iAPAR","iAPAR:Season-long Grazed","iAPAR:Ungrazed","iAPAR:Pulse Grazed",
                                              "iAPAR:Loamy Ecosite","iAPAR:Salt Flats Ecosite","iAPAR:Sandy Ecosite",
                                              "iAPAR:Season-long:Loamy","iAPAR:Season-long:Sandy","iAPAR:Season-long:Salt Flats",
                                              "iAPAR:Ungrazed:Loamy","iAPAR:Ungrazed:Sandy","iAPAR:Ungrazed:Salt Flats",
                                              "iAPAR:Pulse:Loamy","iAPAR:Pulse:Salt Flats","iAPAR:Pulse:Sandy")))%>%

  arrange(Variables)

coefficients.graze2$Variables<- ordered(coefficients.graze2$Variables)
coefficients.graze2 <- coefficients.graze2 %>%
  mutate(label = c("iAPAR", "Grazing", "Grazing", "Grazing", "Ecosite", "Ecosite", "Ecosite",
                   "Grazing:Ecosite","Grazing:Ecosite","Grazing:Ecosite","Grazing:Ecosite",
                   "Grazing:Ecosite","Grazing:Ecosite","Grazing:Ecosite",
                   "Grazing:Ecosite","Grazing:Ecosite"))#%>%
  # 
str(coefficients.graze2)

coefficients.graze2 <- coefficients.graze2 %>%
  mutate(Variables = forcats::fct_rev(Variables))

#coeff graph ----
coeff.gg <- ggplot(coefficients.graze2, aes(x = forcats::fct_inorder(Variables), y = Estimate, ymin = Estimate - (2 * Std.Error), ymax = Estimate + (2 * Std.Error))) +
  geom_linerange(color = "steelblue") +
  geom_point(color = "steelblue", size = 3) +
  mytheme() +
  geom_hline(yintercept = 0,
             colour = "grey60",
             linetype = 2,
             size = 0.8) +
  labs(title = "",
       x = "Variables",
       y = "Coefficient Estimate with 95% CIs") + scale_y_continuous(breaks = breaks_width(2))+
  facet_grid(label~.,scales="free_y")+
  coord_flip()
   # ylim(-2,12)
coeff.gg

# ggsave(filename = "./figures/gg_coeff.pdf",
#        plot = coeff.gg, #this is what you named your plot as, in this case our first plot is g1
#        bg = "transparent",
#        width = 8, height = 8, units = "in",
#        dpi = 300)

# graphing ecosite ----
# Get predicted values for each ecosite
df$predict <- predict(global.model.graze)

# Get predicted values for each ecosite
df$predict_fe <- predict(global.model.graze,re.form=NA)


# gg_ecosite_graze_yr <- ggplot() +
#   geom_point(data = df, aes(x = iAPAR, y = Total_Biomass, color = Ecosite)) +
#   geom_smooth(data = df, aes(x = iAPAR, y = predict , color = Ecosite), size = 1,method = "lm", se = TRUE) +
#   geom_line(data = df, aes(x=iAPAR, y=predict_fe, color=Ecosite))+
#   facet_wrap(~ Year) 
# # stat_poly_line() +
# # stat_poly_eq(data = df,aes(x = iAPAR, y = predict, color = Ecosite))
# 
# # geom_ribbon(data = ggpred.graze, aes(x = x, ymin = (predicted - std.error), ymax = (predicted + std.error), fill = group), alpha = 0.2) 
# gg_ecosite_graze_yr


gg_ecosite_fe <- ggplot() +
  geom_point(data = df, aes(x = iAPAR, y = Total_Biomass, color = Ecosite,shape = Ecosite),size = 1.2,alpha = 0.5) +
  # geom_smooth(data = df, aes(x = iAPAR, y = predict, color = Ecosite), linetype = "solid", size = 1,method = "lm", se = TRUE) +
  geom_smooth(data = df, aes(x = iAPAR, y = predict_fe, color = Ecosite), linetype = "longdash", size = 1,method = "lm", se = TRUE) +
  mytheme() + 
  labs(title = "",x = bquote("iAPAR" (MJ/m^-2)), y = bquote("ANHP" (kg/ha^-1)))+
  scale_shape_manual(values=c(15,17,19))+
  scale_color_manual(values = c("#6460A1", "#1D753D", "#D8731E"))

# geom_ribbon(data = ggpred.graze, aes(x = x, ymin = (predicted - std.error), ymax = (predicted + std.error), fill = group), alpha = 0.2) 
gg_ecosite_fe

gg_ecosite_all <- ggplot() +
  geom_point(data = df, aes(x = iAPAR, y = Total_Biomass, color = Ecosite,shape = Ecosite),size = 1.2,alpha = 0.5) +
  geom_smooth(data = df, aes(x = iAPAR, y = predict, color = Ecosite), linetype = "solid", size = 1,method = "lm", se = FALSE) +
  # geom_smooth(data = df, aes(x = iAPAR, y = predict_fe, color = Ecosite), linetype = "dotdash", size = 1,method = "lm", se = TRUE) +
  mytheme() + 
  labs(title = "",x = bquote("iAPAR" (MJ/m^-2)), y = bquote("ANHP" (kg/ha^-1))) +
  scale_shape_manual(values=c(15,17,19))+
  scale_color_manual(values = c("#6460A1", "#1D753D", "#D8731E"))

# geom_ribbon(data = ggpred.graze, aes(x = x, ymin = (predicted - std.error), ymax = (predicted + std.error), fill = group), alpha = 0.2) 
gg_ecosite_all

grid_ecosite<-plot_grid(gg_ecosite_all, gg_ecosite_fe, labels=c("With random effects", "Without random effects"), ncol = 1, nrow = 2)
grid_ecosite


# graze timing graphs ----
gg_graze_all <- ggplot() +
  geom_point(data = df, aes(x = iAPAR, y = Total_Biomass, color = Graze_timing,shape = Graze_timing),size = 1.2,alpha = 0.5) +
  geom_smooth(data = df, aes(x = iAPAR, y = predict, color = Graze_timing), linetype = "solid", size = 1,method = "lm", se = FALSE) +
  # geom_smooth(data = df, aes(x = iAPAR, y = predict_fe, color = Graze_timing), linetype = "dotdash", size = 1,method = "lm", se = TRUE) +
  mytheme() + 
  labs(title = "",x = bquote("iAPAR" (MJ/m^-2)), y = bquote("ANHP" (kg/ha^-1)))+
  scale_shape_manual(values=c(15,17,19))+
  scale_color_manual(values = c("#E72F85", "#5D879D", "#D8955D"))+
  facet_grid(.~ Ecosite) + geom_ribbon() 
gg_graze_all

gg_graze_fe <- ggplot() +
  geom_point(data = df, aes(x = iAPAR, y = Total_Biomass, color = Graze_timing,shape = Graze_timing),size = 1.2,alpha = 0.5) +
  # geom_smooth(data = df, aes(x = iAPAR, y = predict, color = Graze_timing), linetype = "solid", size = 1,method = "lm", se = FALSE) +
  geom_smooth(data = df, aes(x = iAPAR, y = predict_fe, color = Graze_timing), linetype = "longdash", size = 1,method = "lm", se = FALSE) +
  mytheme() +
  labs(title = "",x = bquote("iAPAR" (MJ/m^-2)), y = bquote("ANHP" (kg/ha^-1)))+
  scale_shape_manual(values=c(15,17,19))+
  scale_color_manual(values = c("#E72F85", "#5D879D", "#D8955D"))+
  facet_grid(.~ Ecosite)
gg_graze_fe 

grid_graze<-plot_grid(gg_graze_all, gg_graze_fe, labels=c("With random effects", "Without random effects"), ncol = 1, nrow = 2)
grid_graze

# ggsave(filename = "./gg_graze_all.pdf",
#        plot = grid_graze, #this is what you named your plot as, in this case our first plot is g1
#        bg = "transparent",
#        width = 8, height = 8, units = "in",
#        dpi = 300)

saltflats <- df %>%
  filter(Ecosite == "Salt Flats")

gg_graze_fe <- ggplot() +
  geom_point(data = saltflats, aes(x = iAPAR, y = Total_Biomass, color = Graze_timing,shape = Graze_timing),size = 1.2,alpha = 0.5) +
  # geom_smooth(data = df, aes(x = iAPAR, y = predict, color = Graze_timing), linetype = "solid", size = 1,method = "lm", se = FALSE) +
  geom_smooth(data = df, aes(x = iAPAR, y = predict, color = Graze_timing), linetype = "longdash", size = 1,method = "lm", se = FALSE) +
  mytheme() #+
  # labs(title = "",x = bquote("iAPAR" (MJ/m^-2)), y = bquote("ANHP" (kg/ha^-1)))+
  # scale_shape_manual(values=c(15,17,19))+
  # scale_color_manual(values = c("#E72F85", "#5D879D", "#D8955D"))
gg_graze_fe

# Model with HEAVY Treatment ----
gg_hvy1 <- ggplot() +
  geom_point(data = df_all, aes(x = iAPAR, y = Total_Biomass, color = Graze_management,shape = Graze_management),size = 1.2,alpha = 0.5) +
  # geom_smooth(data = df, aes(x = iAPAR, y = predict, color = Graze_timing), linetype = "solid", size = 1,method = "lm", se = FALSE) +
  geom_smooth(data = df_all, aes(x = iAPAR, y = Total_Biomass, color = Graze_management), linetype = "longdash", size = 1,method = "lm", se = FALSE) +
  mytheme() +
  labs(title = "",x = bquote("iAPAR" (MJ/m^-2)), y = bquote("ANHP" (kg/ha^-1)))+
  scale_shape_manual(values=c(15,17,19,16))+
  scale_color_manual(values = c("#E72F85", "#5D879D", "#D8955D","purple"))+
  facet_grid(.~ Ecosite)
gg_hvy1

df_all2 <- df_all %>%
  filter(Year %in% c('2019','2020','2021','2022'))%>%
  filter(Ecosite == "Loamy")

gg_hvy2 <- ggplot() +
  geom_point(data = df_all2, aes(x = iAPAR, y = Total_Biomass, color = Graze_management,shape = Graze_management),size = 1.2,alpha = 0.5) +
  # geom_smooth(data = df, aes(x = iAPAR, y = predict, color = Graze_timing), linetype = "solid", size = 1,method = "lm", se = FALSE) +
  geom_smooth(data = df_all2, aes(x = iAPAR, y = Total_Biomass, color = Graze_management), linetype = "longdash", size = 1,method = "lm", se = FALSE) +
  mytheme() +
  labs(title = "",x = bquote("iAPAR" (MJ/m^-2)), y = bquote("ANHP" (kg/ha^-1)))+
  scale_shape_manual(values=c(15,17,19,16))+
  scale_color_manual(values = c("#E72F85", "#5D879D", "#D8955D","#228B22"))+
  facet_grid(.~ Ecosite)
gg_hvy2

all_trt_model<-lmer(paste0(response_var, " ~ iAPAR + Graze_management:iAPAR  + (1 + iAPAR|Year) + (1|Id) "),
                         na.action = "na.fail", data=df_all2)

anova(all_trt_model)
summary(all_trt_model)
plot(all_trt_model)

emm_graze <- emmeans(all_trt_model, pairwise ~ Graze_management)
emm_graze


multcomp::cld(emm_graze, alpha = 0.05, Letters = letters)

pwpp_ecograze <- pwpp(emm_graze_eco)
pwpp_ecograze


# graphing RAP ----
#summary graphs to check infestation rates
# Calculate the slope based on the axis scales
x_range <- range(df2$iAPAR)
y_range <- range(df2$RAP_herb_kg_ha)
slope <- diff(y_range) / diff(x_range)

# Add a 1:1 reference line with the calculated slope
p + geom_abline(intercept = y_range[1] - slope * x_range[1], slope = slope, color = "red")

gg_rap_pred <- ggplot() +
  geom_point(data = df, aes(x = Total_Biomass, y = RAP_herb_kg_ha, color = Graze_timing,shape = Graze_timing),size = 1.2,alpha = 0.5) +
  # geom_smooth(data = df, aes(x = iAPAR, y = predict, color = Graze_timing), linetype = "solid", size = 1,method = "lm", se = FALSE) +
  geom_smooth(data = df, aes(x = Total_Biomass, y = RAP_herb_kg_ha, color = Graze_timing), linetype = "solid", size = 1,method = "lm", se = FALSE) +
  mytheme() +
  ylim(0,3500)+
  labs(title = "",y = bquote("RAP Predicted ANHP" (MJ/m^-2)), x = bquote("ANHP" (kg/ha^-1)))+
  scale_shape_manual(values=c(15,17,19))+
  scale_color_manual(values = c("#E72F85", "#5D879D", "#D8955D"))+
  geom_abline(slope = 1, intercept = 0)+
  facet_grid(.~ Ecosite)
gg_rap_pred

gg_iapar_pred <- ggplot() +
  geom_point(data = df, aes(x = Total_Biomass, y = predict_fe, color = Graze_timing, shape = Graze_timing),size = 1.2,alpha = 0.5) +
  # geom_smooth(data = df, aes(x = iAPAR, y = predict, color = Graze_timing), linetype = "solid", size = 1,method = "lm", se = FALSE) +
  geom_smooth(data = df, aes(x = Total_Biomass, y = predict_fe, color = Graze_timing), linetype = "solid", size = 1,method = "lm", se = FALSE) +
  mytheme() +
  ylim(0,3500)+
  labs(title = "",y = bquote("iAPAR Predicted ANHP" (MJ/m^-2)), x = bquote("ANHP" (kg/ha^-1)))+
  scale_shape_manual(values=c(15,17,19))+
  scale_color_manual(values = c("#E72F85", "#5D879D", "#D8955D"))+
  geom_abline(slope = 1, intercept = 0)+
  facet_grid(.~ Ecosite)
gg_iapar_pred

rap_grid1 <- plot_grid(gg_rap_pred, gg_iapar_pred, labels=c("RAP predictions", "iAPAR predictions"), ncol = 1, nrow = 2)
rap_grid1

# ggsave(filename = "./figures/gg_rap_iapar_pred.pdf",
#        plot = rap_grid1, #this is what you named your plot as, in this case our first plot is g1
#        bg = "transparent",
#        width = 8, height = 8, units = "in",
#        dpi = 300)

# calculate prediction errors ----

### MAPE ----
df$rap_mape <- abs(df$RAP_herb_kg_ha - df$Total_Biomass)/df$Total_Biomass
rap_mape <- mean(df$rap_mape)

df$iapar_mape <- abs(df$predict_fe - df$Total_Biomass)/df$Total_Biomass
iapar_mape <- mean(df$iapar_mape)


# all metrics ----
library(Metrics)
mae(predicted = df$predict_fe, actual = df$Total_Biomass)

error_iapar <- df %>%
  group_by(Ecosite, Graze_timing)%>%
  summarise(#Total_Biomass_mean = mean(Total_Biomass),
            mae = mean(abs(Total_Biomass - predict_fe)),
            mae_pct = (mae / mean(Total_Biomass)) *100,
            mpe = (sum(Total_Biomass - predict_fe) / sum(Total_Biomass)) *100,
            mape = mean(iapar_mape) * 100,
            R2pearson = (cor.test(Total_Biomass,predict_fe, method = "pearson")$estimate)^2
            )
error_iapar
error_iapar$Source <- "iAPAR"

error_iapar_all <- df %>%
  summarise(#Total_Biomass_mean = mean(Total_Biomass),
            mae = mean(abs(Total_Biomass - predict_fe)),
            mae_pct = (mae / mean(Total_Biomass)) *100,
            mpe = (sum(Total_Biomass - predict_fe) / sum(Total_Biomass)) *100,
            mape = mean(iapar_mape) * 100,
            R2pearson = (cor.test(Total_Biomass,predict_fe, method = "pearson")$estimate)^2
  )
error_iapar_all
error_iapar_all$Source <- "iAPAR_all"
error_iapar_all$Ecosite <- "All"
error_iapar_all$Graze_timing <- "All"

error_rap <- df %>%
  group_by(Ecosite,Graze_timing)%>%
  summarise(mae = mean(abs(Total_Biomass - RAP_herb_kg_ha)),
            mae_pct = (mae / mean(Total_Biomass)) *100,
            mpe = (sum(Total_Biomass - RAP_herb_kg_ha) / sum(Total_Biomass)) *100,
            mape = mean(rap_mape) *100,
            R2pearson = (cor.test(Total_Biomass,RAP_herb_kg_ha, method = "pearson")$estimate)^2
            ) 
error_rap
error_rap$Source <- "RAP"

error_rap_all <- df %>%
  summarise(mae = mean(abs(Total_Biomass - RAP_herb_kg_ha)),
            mae_pct = (mae / mean(Total_Biomass)) *100,
            mpe = (sum(Total_Biomass - RAP_herb_kg_ha) / sum(Total_Biomass)) *100,
            mape = mean(rap_mape) *100,
            R2pearson = (cor.test(Total_Biomass,RAP_herb_kg_ha, method = "pearson")$estimate)^2
  ) 
error_rap_all
error_rap_all$Source <- "RAP_all"
error_rap_all$Ecosite <- "All"
error_rap_all$Graze_timing <- "All"

# error_combined$Graze_Ecosite <- paste(error_combined$Graze_timing, "_", error_combined$Ecosite)

error_combined <- data.frame()
error_combined <- rbind(error_iapar, error_rap)
error_combined <- rbind(error_combined, error_rap_all)
error_combined <- rbind(error_combined, error_iapar_all)
error_combined <- error_combined %>%
  select(Source, everything())
error_combined

library(kableExtra)
library(flextable)
library(dplyr)
library(forcats)
library(tidyverse)
library(flextable)
library(officer)
library(scales)

my_theme_tbl <- function(x, ...) {
  x <- colformat_double(x, big.mark = "'", decimal.mark = ",", digits = 1)
  x <- set_table_properties(x, layout = "fixed")
  x <- border_remove(x)
  std_border <- fp_border(width = 1, color = "grey")
  x <- border_outer(x, part="all", border = std_border )
  x <- border_inner_h(x, border = std_border, part="all")
  x <- hline_top(x, part = "body", 
                 border = fp_border(color = "black", 
                                    width = 2, 
                                    style = "solid"))
  x <- colformat_num(x, decimal.mark = ".")
  x <- colformat_double(x, digits = 4)
  autofit(x)
}

set_flextable_defaults(
  font.color = "black",
  border.color = "black",
  theme_fun = "my_theme_tbl")

colourer <- col_numeric(
  palette = c("transparent", "red"),
  domain = c(0, 1))

ft1 <-flextable(
  error_combined,
  col_keys = names(error_combined),
  cwidth = 0.75,
  cheight = 0.20,
  defaults = list()
)

ft1

ft1 <- ft1 %>% 
  merge_v(j="Ecosite")%>%
  merge_v(j="Source")%>%
  # valign(j=c("Ecosite"),valign = "top")%>%
  # align(j=c("Ecosite"),align = "center",part = "body")%>%
  align(j=1:8,align = "center", part = "all")%>%
  hline(part = "body", i = 2:8)#%>%
  # add_header_row(
  #   top = TRUE,                # New header goes on top of existing header row
  #   values = c("",     # Header values for each column below
  #              "",
  #              "",    # This will be the top-level header for this and two next columns
  #              "",
  #              "",
  #              "Error Metrics",         # This will be the top-level header for this and two next columns
  #              ""))%>%
  # merge_at(i = 1, j = 4:8, part = "header")
  # bg(bg = colourer,
  #    j = "Ecosite",
  #    part = "body")
ft1

# save_as_docx(
#   "RAP vs iAPAR production predictions" = ft1,
#   path = "./tables/rap_vs_iAPAR.docx")

# Graphing error metrics ----
gg_metrics <- ggplot(error_combined, aes(x = Graze_timing, y = mae, fill= Source))+
  geom_col(position = "dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.3)) +
  facet_grid(~ Ecosite) + mytheme()
gg_metrics

# Graphing error ----

# Create a scatterplot to compare model errors
rap_error.gg <- ggplot(df, aes(x = Total_Biomass)) +
  geom_point(aes(y = rap_error, color = Graze_timing)) +
  # geom_point(aes(y = iapar_error, color = ifelse(rap_mape > iapar_mape, "RAP Overpredicts", "IAPAR Overpredicts"))) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(. ~ Ecosite) +
  labs(
    x = "RAP predicted",
    y = "Residuals",
    color = "Comparison"
  ) +
    ylim(-1000,2500)+ xlim(0,3500)+
  mytheme()

rap_error.gg

iapar_error.gg <- ggplot(df, aes(x = Total_Biomass)) +
  geom_point(aes(y = iapar_error, color = Graze_timing)) +
  # geom_point(aes(y = iapar_error, color = ifelse(rap_mape > iapar_mape, "RAP Overpredicts", "IAPAR Overpredicts"))) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(. ~ Ecosite) +
  labs(
    x = "iAPAR predicted",
    y = "Residuals",
    color = "Comparison"
  ) +
    ylim(-1000,2500)+xlim(0,3500)+
  mytheme()
iapar_error.gg

error_grid <- plot_grid(rap_error.gg,iapar_error.gg,ncol = 1, nrow = 2)
error_grid

# graphing NDVI and grazing date ----
grazing_dates <- read.csv("../data/ground/grazing/CARM_ActualGrazingInfov3_2013-2023.csv")
apar<- read.csv("../data/training/iapar/cper_all_year_iapar_2014_2022.csv")
sos <- read.csv("../data/training/iapar/cper_sos_2014_2022.csv")
# cleaning sos
sos$SOS_date <- ymd(sos$SOS_date)
sos <- sos %>%
  select(Year, Id, SOS_date)%>%
  separate(Id, into = c("Pasture","Plot"), sep = "_")%>%
  group_by(Pasture)%>%
  summarize(SOS_date = mean(SOS_date))
  

#cleaning iapar
apar$Date <- ymd(apar$Date)
apar <- apar %>%
  select(Id,Year,Date, APAR_adjusted,NDVI_smooth_avg)%>%
  separate(Id, into = c("Pasture","Plot"), sep = "_")%>%
  group_by(Pasture,Year,Date)%>%
  summarize(APAR_mean = mean(APAR_adjusted),NDVI_mean = mean(NDVI_smooth_avg))
#   filter(!is.na(Date))
  
str(apar)

# cleaning grazing dates
grazing_dates$Pasture <- grazing_dates$PastureCode
grazing_dates$DateInPasture <- as.Date(grazing_dates$DateInPasture, format = "%m/%d/%Y")
grazing_dates$DateOutPasture <- as.Date(grazing_dates$DateOutPasture, format = "%m/%d/%Y")
grazing_dates$Year <- as.factor(grazing_dates$Year)
grazing_dates <- grazing_dates %>%
  select(Year,DateInPasture,DateOutPasture,Pasture)


#merging and cleaning data
df_select <- df %>%
  distinct(Year,Pasture,Graze_timing,Treatment,Date)%>%
  rename(SamplingDate = Date)

grazing_timing <- right_join(df_select, grazing_dates, by = c("Pasture","Year"))

#adding graze timing dates and factors to filtered data set
filtered_df2 <- merge(apar,grazing_timing, by = c("Pasture","Year"),all.x = TRUE)

#adding sos date to filtered data set
filtered_df2 <- merge(filtered_df2,sos, by = c("Pasture","Year"),all.x = TRUE)

# str(filtered_ndvi)
filtered_df2$Date <- ymd(filtered_df2$Date)
filtered_df3 <- filtered_df2 %>%
  filter(Year == "2021") %>%
  filter(!is.na(Treatment))%>%
  # filter(!is.na(DateInPasture))%>%
  ungroup()
filtered_df3$Date <- ymd(filtered_df3$Date)
filtered_df3$SOS_date <- ymd(filtered_df3$SOS_date)
filtered_df3$SamplingDate<- ymd(filtered_df3$SamplingDate)
str(filtered_df3)


# Split your data frame into a list of data frames for each Year/Pasture combination
df_list <- sample(split(filtered_df3, list(filtered_df3$Id, filtered_df3$Year)))


# Obtain the overall min and max values for NDVI_smooth_avg and APAR_modified
overall_min_ndvi <- min(sapply(df_list, function(filtered_df3) min(filtered_df3$NDVI_smooth_avg, na.rm = TRUE)))
overall_max_ndvi <- max(sapply(df_list, function(filtered_df3) max(filtered_df3$NDVI_smooth_avg, na.rm = TRUE)))

overall_min_apar <- min(sapply(df_list, function(filtered_df3) min(filtered_df3$APAR_adjusted, na.rm = TRUE)))
overall_max_apar <- max(sapply(df_list, function(filtered_df3) max(filtered_df3$APAR_adjusted, na.rm = TRUE)))

generate_plots <- function(df) {
  # Create NDVI_smooth_avg plot

  plot_ndvi <- ggplot(df, aes(x = Date, y = NDVI_smooth_avg)) +
    geom_line() +
    scale_color_brewer(palette = "Set1") +
    scale_y_continuous(limits = c(overall_min_ndvi, overall_max_ndvi)) +  # standardize y-axis
    theme_minimal() +
    labs(title = paste(unique(df$Year), unique(df$Id), "smoothed NDVI", sep = " / "),
         x = "Date",
         y = "smoothed NDVI") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
    geom_vline(xintercept = df$DateInPasture, color = "green") +
    geom_vline(xintercept = df$DateOutPasture, color = "red")+ 
    mytheme() 
  
  # Create APAR plot
  # Create a subset for the shaded region
  shaded_region <- df[df$Date >= df$SOS_date & df$Date <= df$SamplingDate,]
  # 
  # Create APAR plot
  plot_apar <- ggplot(df, aes(x = Date, y = APAR_adjusted)) +
    geom_line() +
    geom_ribbon(data = shaded_region, aes(x = Date, ymin = overall_min_apar, ymax = APAR_adjusted), fill = "grey") +
    geom_vline(xintercept = as.numeric(df$SOS_date), linetype = "dashed", color = "red") + # Add SOS_date line
    geom_vline(xintercept = as.numeric(df$SamplingDate), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = df$DateInPasture, color = "green") +
    geom_vline(xintercept = df$DateOutPasture, color = "red")+ 
        scale_color_brewer(palette = "Set1") +
    scale_y_continuous(limits = c(overall_min_apar, overall_max_apar)) +  # standardize y-axis # standardize y-axis
    mytheme() +
    labs(title = paste(unique(df$Year), unique(df$Id), "APAR", sep = " / "),
         x = "Date",
         y = "Adjusted APAR") +
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
  # Return list of plots
  list(plot_ndvi, plot_apar)
}

invisible(lapply(df_list, function(filtered_df3) dim(filtered_df3)[1]))


# Generate plots for each Year/Pasture combination
plots_list <- lapply(df_list, generate_plots)


# Create a list to save grid of plots for each page
plot_pages <- list()

# Number of pages
num_pages <- ceiling(length(plots_list) / 4)
p1 <- plots_list["21N_P1.2021"][[1]]

p2 <- plots_list["24W_P1.2021"][[1]]
p2
p1 <- plots_list["21N_P1.2021"][[1]][2]
p1
p2 <- plots_list["17S_P1.2021"][[1]][2]
p2
p3 <- plots_list["31E_P1.2021"][[1]][2]
p3
p4 <- plots_list["18S_P1.2021"][[1]][2]
p4

grid_apar <- plot_grid(p1[[1]],p2[[1]],p3[[1]],p4[[1]])
grid_apar

# ggsave(filename = "./figures/gg_apar_examples.pdf",
#        plot = grid_apar, #this is what you named your plot as, in this case our first plot is g1
#        bg = "transparent",
#        width = 8, height = 8, units = "in",
#        dpi = 300)

# Iterate over each page
for (i in 1:num_pages) {
  # Select plots for this page
  plots <- plots_list[((i - 1) * 4 + 1):(i * 4)]
  
  # Flatten the list of lists
  plots <- unlist(plots, recursive = FALSE)
  
  # Arrange plots into a grid
  grid <- ggarrange(plotlist = plots, ncol = 2, nrow = 2)
  
  # Add grid to list of pages
  plot_pages[[i]] <- grid
}

# Create a new PDF file
pdf("APAR_timeseries.pdf")

# Iterate over plot_pages and print each one to the PDF
for(page in plot_pages) {
  print(page)
}

# Close the PDF file
invisible(dev.off())


