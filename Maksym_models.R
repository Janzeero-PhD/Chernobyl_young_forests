library(tidyverse)
library(lme4)
library(DHARMa)
library(effects)
library(vegan)

df.all.heights <- read_csv('https://raw.githubusercontent.com/Janzeero-PhD/Chernobyl_young_forests/master/df_with_heights.csv')
df.all.added <- read_csv('https://raw.githubusercontent.com/Janzeero-PhD/Chernobyl_young_forests/master/df_trees.csv')
df.plots <- read_csv('https://raw.githubusercontent.com/Janzeero-PhD/Chernobyl_young_forests/master/df_plots.csv')

df.all.heights$IDPlots <- as.factor(df.all.heights$IDPlots)
df.all.added$IDPlots <- as.factor(df.all.added$IDPlots)
df.plots$id <- as.factor(df.plots$id)

# N large trees

df.plots <- df.plots %>%
  mutate(response = log(n_large)) %>%
  mutate(dist_stand_scaled = scale(dist),
         prop_stand_scaled = scale(prop),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

ggplot(df.plots, aes(x = log(n_large))) +
  geom_histogram()

ggplot(df.plots, aes(x = cont_code, y = log(n_large))) +
  geom_point(position = position_jitter())

fit_null <- lm(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, 
               data = df.plots %>% filter(., response > 0)
)
fit_n_large <- lm(response ~ cont_code + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, 
                  data = df.plots %>% filter(., response > 0)
)


simulationOutput <- simulateResiduals(fittedModel = fit_n_large, n = 250)

plot(simulationOutput)

hist(residuals(fit_n_large))

anova(fit_null, fit_n_large)

summary(fit_n_large)

effects <- effect("cont_code", fit_n_large) %>%
  as.data.frame(.)

ggplot(effects, aes(x = factor(cont_code, levels = c("0", "1", "2")),
                    y = exp(fit))) +
  geom_point(data = df.plots, 
             aes(x = factor(cont_code, levels = c("0", "1", "2")), y = n_large), 
             col = "grey", position = position_jitter(width = 0.1)) +
  geom_point() +
  geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper)), width = 0.2)

# N understory

df.plots <- df.plots %>%
  mutate(response = log(n_micro)) %>%
  mutate(dist_stand_scaled = scale(dist),
         prop_stand_scaled = scale(prop),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

ggplot(df.plots, aes(x = (n_micro))) +
  geom_histogram()

ggplot(df.plots, aes(x = cont_code, y = log(n_large))) +
  geom_point(position = position_jitter())

fit_null <- glm(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, 
               data = df.plots %>% filter(., response > 0), family = poisson(link = "log")
)
fit_n_micro <- glm(response ~ cont_code + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, 
                  data = df.plots %>% filter(., response > 0), family = poisson(link = "log")
)

simulationOutput <- simulateResiduals(fittedModel = fit_n_micro, n = 250)

plot(simulationOutput) # totally ugly

hist(residuals(fit_n_micro))

anova(fit_null, fit_n_micro)

summary(fit_n_micro)

effects <- effect("cont_code", fit_n_micro) %>%
  as.data.frame(.)

ggplot(effects, aes(x = factor(cont_code, levels = c("0", "1", "2")),
                    y = exp(fit))) +
  geom_point(data = df.plots, 
             aes(x = factor(cont_code, levels = c("0", "1", "2")), y = n_micro), 
             col = "grey", position = position_jitter(width = 0.1)) +
  geom_point() +
  geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper)), width = 0.2)

## BA_pine_prop

df.plots <- df.plots %>%
  mutate(response = (prop_BA_pine)) %>%
  mutate(dist_stand_scaled = scale(dist),
         prop_stand_scaled = scale(prop),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

ggplot(df.plots, aes(x = (prop_BA_pine))) +
  geom_histogram()

ggplot(df.plots, aes(x = cont_code, y = (prop_BA_pine))) +
  geom_point(position = position_jitter())

library(betareg)
y.transf.betareg <- function(y){ # function to transform (0..1) distribution to the one suitable for betareg to use
  n.obs <- sum(!is.na(y))
  (y * (n.obs - 1) + 0.5) / n.obs
}

fit_null <- betareg(y.transf.betareg(response) ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, 
               data = df.plots #%>% filter(., response > 0)
)
fit_BA_pine <- betareg(y.transf.betareg(response) ~ cont_code + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, 
                  data = df.plots #%>% filter(., response > 0)
)


simulationOutput <- simulateResiduals(fittedModel = fit_BA_pine, n = 250) # cannot run DHARMA for betareg

plot(residual(fit_BA_pine, type = 'response')) # fit raw residuals of the betareg model

anova(fit_null, fit_BA_pine) # cannot run anova for betareg

library(lmtest)
lrtest(fit_null, fit_BA_pine) # likelihood test... 0.15 = too high :(

summary(fit_BA_pine)

effects <- effect("cont_code", fit_BA_pine) %>%
  as.data.frame(.)

ggplot(effects, aes(x = factor(cont_code, levels = c("0", "1", "2")),
                    y = (fit))) +
  geom_point(data = df.plots, 
             aes(x = factor(cont_code, levels = c("0", "1", "2")), y = prop_BA_pine), 
             col = "grey", position = position_jitter(width = 0.1)) +
  geom_point() +
  geom_errorbar(aes(ymin = (lower), ymax = (upper)), width = 0.2)
ggsave('pine_BA_effects.jpg')

## BA_birch_prop

df.plots <- df.plots %>%
  mutate(response = (prop_BA_birch)) %>%
  mutate(dist_stand_scaled = scale(dist),
         prop_stand_scaled = scale(prop),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

ggplot(df.plots, aes(x = (prop_BA_birch))) +
  geom_histogram()

ggplot(df.plots, aes(x = cont_code, y = (prop_BA_birch))) +
  geom_point(position = position_jitter())

fit_null <- lm(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, 
               data = df.plots #%>% filter(., response > 0)
)
fit_BA_birch <- lm(response ~ cont_code + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, 
                   data = df.plots #%>% filter(., response > 0)
)


simulationOutput <- simulateResiduals(fittedModel = fit_BA_birch, n = 250)

plot(simulationOutput)

hist(residuals(simulationOutput))

anova(fit_null, fit_BA_birch)

summary(fit_BA_birch)

effects <- effect("cont_code", fit_BA_birch) %>%
  as.data.frame(.)

ggplot(effects, aes(x = factor(cont_code, levels = c("0", "1", "2")),
                    y = (fit))) +
  geom_point(data = df.plots, 
             aes(x = factor(cont_code, levels = c("0", "1", "2")), y = prop_BA_birch), 
             col = "grey", position = position_jitter(width = 0.1)) +
  geom_point() +
  geom_errorbar(aes(ymin = (lower), ymax = (upper)), width = 0.2)

# height diversity
ggplot(df.plots, aes(x = shannon_height)) +
  geom_histogram()

df.plots <- df.plots %>%
  mutate(response = shannon_height) %>%
  mutate(dist_stand_scaled = scale(dist),
         prop_stand_scaled = scale(prop),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

fit_null <- lm(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, data = df.plots)
fit_shannon_height <- lm(response ~ cont_code + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, data = df.plots)

simulationOutput <- simulateResiduals(fittedModel = fit_shannon_height, n = 250)

plot(simulationOutput)

anova(fit_null, fit_shannon_height)

effects <- effect("cont_code", fit_shannon_height) %>%
  as.data.frame(.)

ggplot(effects, aes(x = factor(cont_code, levels = c("0", "1", "2")),
                    y = fit )) +
  geom_point(data = df.plots, 
             aes(x = factor(cont_code, levels = c("0", "1", "2")), y = shannon_height), 
             col = "grey", position = position_jitter(width = 0.1)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower , ymax = upper ), width = 0.2)

# heights:
ggplot(df.all.heights, aes(x = log(Height_m))) +
  geom_histogram()

df.all.heights <- df.all.heights %>%
  mutate(response = log(Height_m)) %>%
  mutate(dist_stand_scaled = scale(dist_stand),
         prop_stand_scaled = scale(prop_stand),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))
df.all.heights$IDPlots <- as.factor(df.all.heights$IDPlots)

fit_null <- lmer(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                 data = df.all.added %>% filter(., DBH > 0), REML = FALSE)

fit_height <- lmer(response ~ Cont + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                   data = df.all.added %>% filter(., DBH > 0), REML = FALSE)

simulationOutput <- simulateResiduals(fittedModel = fit_height, n = 250)

plot(simulationOutput)

hist(residuals(fit_height))

anova(fit_null, fit_height)

effects <- effect("Cont", fit_height) %>%
  as.data.frame(.)

ggplot(effects, aes(x = factor(Cont, levels = c("Low", "Medium", "High")),
                    y = exp(fit))) +
  geom_point() +
  geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper)), width = 0.2)
ggsave('height_effect.jpg')
