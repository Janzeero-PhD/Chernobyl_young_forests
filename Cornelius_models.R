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

### Models:

# Mortality

df.all.added <- df.all.added %>%
  mutate(response = alive == "no") %>%
  mutate(dist_stand_scaled = scale(dist_stand),
         prop_stand_scaled = scale(prop_stand),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

fit <- glmer(response ~ Cont + (1 | IDPlots), data = df.all.added, family = binomial(link = "logit"))

summary(fit)

predict(fit, newdata = data.frame(Cont = unique(df.all.added$Cont)), re.form = NA, type = "response")

fit_null <- glmer(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                  data = df.all.added, family = binomial(link = "logit"))

fit_mortal <- glmer(response ~ Cont + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), data = df.all.added, family = binomial(link = "logit"))

simulationOutput <- simulateResiduals(fittedModel = fit_mortal, n = 250)

plot(simulationOutput)

anova(fit_null, fit_mortal)

summary(fit_more)

predict(fit_more, newdata = expand.grid(Cont = unique(df.all.added$Cont),
                                        dist_stand_scaled = 0,
                                        prop_stand_scaled = 0,
                                        elevation_scaled = 0,
                                        dist_water_scaled = 0), 
        re.form = NA, type = "response")


effects <- effect("Cont", fit_more) %>%
  as.data.frame(.)

ggplot(effects, aes(x = factor(Cont, levels = c("Low", "Medium", "High")),
                    y = fit)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)

# DBH

df.all.added <- df.all.added %>%
  mutate(response = log(DBH)) %>%
  mutate(dist_stand_scaled = scale(dist_stand),
         prop_stand_scaled = scale(prop_stand),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

ggplot(df.all.added, aes(x = log(DBH))) +
  geom_histogram()


ggplot(df.all.added, aes(x = Cont, y = (DBH))) +
  geom_point(position = position_jitter())

fit_null <- lmer(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                 data = df.all.added %>% filter(., DBH > 0), REML = FALSE)

fit_dbh <- lmer(response ~ Cont + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                data = df.all.added %>% filter(., DBH > 0), REML = FALSE)

simulationOutput <- simulateResiduals(fittedModel = fit_dbh, n = 250)

plot(simulationOutput)

hist(residuals(fit_dbh))

anova(fit_null, fit_dbh)

summary(fit_dbh)

effects <- effect("Cont", fit_dbh) %>%
  as.data.frame(.)

ggplot(effects, aes(x = factor(Cont, levels = c("Low", "Medium", "High")),
                    y = exp(fit))) +
  geom_point() +
  geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper)), width = 0.2)

# Deformation

df.all.added <- df.all.added %>%
  mutate(response = deformated == "no") %>%
  mutate(dist_stand_scaled = scale(dist_stand),
         prop_stand_scaled = scale(prop_stand),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

fit_null <- glmer(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                  data = df.all.added, family = binomial(link = "logit"))
fit_deform <- glmer(response ~ Cont + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), data = df.all.added, 
                    family = binomial(link = "logit"))


simulationOutput <- simulateResiduals(fittedModel = fit_deform, n = 250)

plot(simulationOutput)

anova(fit_null, fit_deform)

summary(fit_more)

effects <- effect("Cont", fit_deform) %>%
  as.data.frame(.)

ggplot(effects, aes(x = factor(Cont, levels = c("Low", "Medium", "High")),
                    y = fit)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)

# Shannon diam

ggplot(df.plots, aes(x = shannon_diam)) +
  geom_histogram()

df.plots <- df.plots %>%
  mutate(response = shannon_diam) %>%
  mutate(dist_stand_scaled = scale(dist),
         prop_stand_scaled = scale(prop),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

fit_null <- lm(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, data = df.plots)
fit_shannon_diam <- lm(response ~ cont_code + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, data = df.plots)

simulationOutput <- simulateResiduals(fittedModel = fit_shannon_diam, n = 250)

plot(simulationOutput)

anova(fit_null, fit_shannon_diam)

effects <- effect("cont_code", fit_shannon_diam) %>%
  as.data.frame(.)

ggplot(effects, aes(x = factor(cont_code, levels = c("0", "1", "2")),
                    y = fit )) +
  geom_point(data = df.plots, 
             aes(x = factor(cont_code, levels = c("0", "1", "2")), y = shannon_diam), 
             col = "grey", position = position_jitter(width = 0.1)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower , ymax = upper ), width = 0.2)

# Shannon species

ggplot(df.plots, aes(x = shannon_exp)) +
  geom_histogram()

df.plots <- df.plots %>%
  mutate(response = round(shannon_exp * 100, 0)) %>%
  mutate(dist_stand_scaled = scale(dist),
         prop_stand_scaled = scale(prop),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

fit_null <- MASS::glm.nb(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, 
                         data = df.plots %>% filter(., response > 0))
fit_shannon_species <- MASS::glm.nb(response ~ cont_code + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, 
                                    data = df.plots %>% filter(., response > 0))

simulationOutput <- simulateResiduals(fittedModel = fit_shannon_species, n = 250)

plot(simulationOutput)

anova(fit_null, fit_shannon_species)

effects <- effect("cont_code", fit_shannon_species) %>%
  as.data.frame(.)

ggplot(effects, aes(x = factor(cont_code, levels = c("0", "1", "2")),
                    y = fit / 100)) +
  geom_point(data = df.plots, 
             aes(x = factor(cont_code, levels = c("0", "1", "2")), y = shannon_exp), 
             col = "grey", position = position_jitter(width = 0.1)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower / 100, ymax = upper / 100), width = 0.2)