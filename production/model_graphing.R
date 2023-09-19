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

inPATH = '../data/training/iapar/model_selection_data.csv'

df = read.csv(inPATH)
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

#Using a lmer ----
global.model.graze<-lmer(paste0(response_var, " ~ iAPAR + Graze_timing:iAPAR + Ecosite:iAPAR + Graze_timing:Ecosite:iAPAR  + (1|Year) + (1|Id)"),
                         na.action = "na.fail", data=df)

anova(global.model.graze)
summary(global.model.graze)
plot(global.model.graze)

emm_graze <- emmeans(global.model.graze, pairwise ~ Graze_timing)
emm_graze

emm_eco<- emmeans(global.model.graze, pairwise ~ Ecosite)
emm_eco

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
par(mfrow=c(2,2))
hist(df$Total_Biomass)
hist(df$Total_Biomass_sqrt)

# Identify outliers based on Cook's distance and leverage
outliers <- which(cooksd > 1 | leverage > 2)
outliers
# Check for normality of residuals
hist(residuals(global.model.graze))
qqnorm(residuals(global.model.graze))
shapiro.test(residuals(global.model.graze))

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

# Create the plot using ggplot2
library(ggplot2)

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


# making fake data for ecosite no intercept----
fake_iapar <- as.numeric(seq(0, 250)) #adjust to range of iAPAR in data

cf_iapar <- 9.4174

fake_loamy <- fake_iapar * cf_iapar

cf_iapar_sandy <- 2.6775
fake_sandy <- (fake_iapar * cf_iapar) +  (fake_iapar * cf_iapar_sandy)


cf_iapar_salt <- 5.2734
fake_salt <- (fake_iapar * cf_iapar) +  (fake_iapar * cf_iapar_salt)



df_fake <- data.frame(fake_iapar, fake_loamy, fake_sandy, fake_salt)
names(df_fake) <- c('iAPAR', 'Loamy', 'Sandy', 'Salt Flats')


df_fake <- gather(df_fake,  key='Ecosite', value='Biomass_pred', Loamy:"Salt Flats", factor_key=TRUE)

# making fake data for ecosite WITH intercept----
fake_iapar2 <- as.numeric(seq(0, 250)) #adjust to range of iAPAR in data

cf_iapar <- 9.4174

fake_loamy2 <- 152.2067 + fake_iapar * cf_iapar

cf_iapar_sandy <- 2.6775
fake_sandy2 <- 152.2067 + (fake_iapar * cf_iapar) +  (fake_iapar * cf_iapar_sandy)


cf_iapar_salt <- 5.2734
fake_salt2 <- 152.2067 + (fake_iapar * cf_iapar) +  (fake_iapar * cf_iapar_salt)


df_fake2 <- data.frame(fake_iapar2, fake_loamy2, fake_sandy2, fake_salt2)
names(df_fake2) <- c('iAPAR', 'Loamy', 'Sandy', 'Salt Flats')


df_fake2 <- gather(df_fake2,  key='Ecosite', value='Biomass_pred', Loamy:"Salt Flats", factor_key=TRUE)


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
# ggsave(filename = "./03_plots/grid_weights.pdf",
#        plot = grid_weights, #this is what you named your plot as, in this case our first plot is g1
#        bg = "transparent",
#        width = 8, height = 8, units = "in",
#        dpi = 300)

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
