# Library Load ----

library(tidyverse)
library(car)
library(AICcmodavg)
library(lme4)
library(lubridate)
library(lmerTest)
library(merTools) 

# Data ----
# library(readxl)
# lake <- read_excel("Montana Lakes Data FINAL.xlsx")
lakes <- read.csv("lake_data.csv")


# . Pulling apart the year, month, day ----
# Puts this column into a date format R can recognize
lakes$DATE <- as.Date(as.character(lakes$DATE), format = "%m/%d/%Y")  
lakes$Year <- year(lakes$DATE)  # Makes a column for Year only
lakes$year.s <- scale(lakes$Year)
lakes$Month <- month(lakes$DATE) # Makes a column for Month only
lakes$Month <- as.character(lakes$Month)  # Changes Month to categorical


# Puts all Flathead lake.names together ----
lakes$lake.name <- lakes$LAKE.NAME
lakes$site =  lakes$LAKE.NAME                                      # Making variable name easier to deal with
lakes$lake.name[grep("Flathead", lakes$lake.name)] <- "Flathead"

# Remove trailing whitespace from lake names ----
lakes$lake.name <- trimws(lakes$lake.name, "right")

# Data Exploration ----
# Convert variables to numbers
# lakes$SECCHI <- as.numeric(lake$SECCHI)
# lakes$temp <- as.numeric(lake$`LAKE TEMP`)

# How many measurements by month each year?
with(lakes, table(Year, Month))

# . Secchi data cleanup ----
# Filter the data
lake <- lakes %>% 
  # This 156 foot reading was an extreme outlier
  filter(!is.na(SECCHI) & 
           SECCHI < 40 & 
           Month %in% c("6", "7", "8"))

# . Temperature data cleanup ----
# Set aside data for temperature analysis
# so we don't have to do all the date manipulation
# again
lake_temp <- lakes
lake_temp$temp <- lake_temp$LAKE.TEMP

# Have some NA values in response to remove
any(is.na(lake_temp$temp))

lake_temp <- lake_temp %>% 
  filter(!is.na(temp) & temp > 0 & temp < 40 &
           Month %in% c("6", "7", "8"))

# This code allows us to select minimum number of years
# for lakes that we want to include
check <- with(lake_temp, table(lake.name, Year))
pres_check <- check
pres_check[pres_check > 0] <- 1
rowSums(check)
samples <- rowSums(pres_check)

years_sampled <- data.frame(names(samples), samples)
names(years_sampled) <- c("lake.name", "n_years")

enough_data <- years_sampled[years_sampled$n_years >=5, ]

lake_temp <- lake_temp[lake_temp$lake.name %in%
                         enough_data$lake.name, ]



# Secchi Models ----
# Have a look at secchi
hist(lake$SECCHI, breaks = 60)
hist(log(lake$SECCHI))

ggplot(lake, aes(x = SECCHI))+
  geom_histogram()+
  facet_wrap(~lake.name)

# Null Model
null <- lm(log(SECCHI) ~ 1, data = lake)

# Effect of year.s plus Month on Secchi depth.
year.s.month <- lm(log(SECCHI) ~ year.s + Month, data = lake)

# Effect of year.s and Month together on Secchi depth.
year.s.x.month <- lm(log(SECCHI) ~ year.s * Month, data = lake)

# . Model selection on fixed effects ----
# AIC Rankings (choosing the best model)
model.list <- list(null, year.s.month, year.s.x.month)
model.names <- c("null", "year.s.month", "year.s.x.month")

aictab(cand.set = model.list, modnames = model.names)


# . Model selection on random effects ----
# Interaction models only to select best random effects structure
# Random effect on intercept
int.1 <- lmer(log(SECCHI) ~ year.s * Month + (1|lake.name), data = lake)
# Random effect on intercept and slope
int.2 <- lmer(log(SECCHI) ~ year.s * Month + (year.s|lake.name), data = lake)


model.list3 <- list(int.1, int.2)
model.names3 <- c("int.1", "int.2")

aictab(cand.set = model.list3, modnames = model.names3)

# . Residual check ----
# Note that a sqrt() transformation cleans up normality by group
# and reduces number of outliers, but it does not seem to make a 
# whole lot of difference in the inference or in the predictions...
# could go either way
resids = residuals(int.2)
resid.test = data.frame(lake, resids)
boxplot(resids~lake.name, data = resid.test); abline(h = 0)


# . Statistical significance ----
# Year not statistically significant, but month is
Anova(int.2, type = "III")
summary(int.2)


# . Statistical significance ----
Anova(int.2, type = "III")
summary(int.2)


# . Random effects coeffs ----
# Here is how we get mean and 95% CIs for the random
# intercept and slope adjustments out of this monster.

# .. Get the fixed and random effects ----
cf <- data.frame(summary(int.2)$coefficients) # Fixed
coeffs <- data.frame(ranef(int.2))            # Random

# Pivot the table to unstack parameters and rename
coeffs <- pivot_wider(coeffs, 
                      names_from = term,
                      values_from = c(condval, condsd)
                      )
names(coeffs)[3:6] <- c("Intercept", "June", 
                        "Intercept_sd", "June_sd")

# Get slope adjustments and sds for July and August
# Mean slope adjustment
coeffs$July <- cf[5, 1]
coeffs$July_sd <- cf[5, 2]
coeffs$August <- cf[6, 1]
coeffs$August_sd <- cf[6, 2]

# Add fixed effect coefficients to coeffs table
coeffs$fixed_int <- cf[1, 1]
coeffs$fixed_int_sd <- cf[1, 2]
coeffs$fixed_slope <- cf[2, 1]
coeffs$fixed_slope_sd <- cf[2, 2]

# June is the baseline slope in this model

# # July slope needs to be adjusted for (fixed) main effects
# # and interaction
# coeffs$july_adj <- coeffs$fixed_slope + coeffs$July + coeffs$June 
# coeffs$july_sd_adj <- coeffs$fixed_slope_sd + coeffs$July_sd + coeffs$June_sd 
# 
# # So does August
# coeffs$august_adj <- coeffs$fixed_slope + coeffs$August + coeffs$June 
# coeffs$august_adj_sd <- coeffs$fixed_slope_sd + coeffs$August_sd + coeffs$June_sd 

# .. Calculate derived parameters in a new df ----
derived <- data.frame(coeffs)

# ... Intercepts and CIs ----
derived$Intercept <- coeffs$Intercept + coeffs$fixed_int
derived$Intercept_lwr <- coeffs$Intercept + coeffs$fixed_int - 
  1.96 * (coeffs$Intercept_sd + coeffs$fixed_int_sd)
derived$Intercept_upr <- coeffs$Intercept + coeffs$fixed_int + 
  1.96 * (coeffs$Intercept_sd + coeffs$fixed_int_sd)

# ... Slopes and CIs ----
# June
derived$June <- coeffs$June + coeffs$fixed_slope 
derived$June_lwr <- coeffs$June + coeffs$fixed_slope - 
  1.96 * (coeffs$June_sd + coeffs$fixed_slope_sd)
derived$June_upr <- coeffs$June + coeffs$fixed_slope + 
  1.96 * (coeffs$June_sd + coeffs$fixed_slope_sd)

# July
derived$July <- coeffs$June + coeffs$July + coeffs$fixed_slope 
derived$July_lwr <- coeffs$June + coeffs$July + coeffs$fixed_slope - 
  1.96 * (coeffs$July_sd + coeffs$fixed_slope_sd + coeffs$June_sd )
derived$July_upr <- coeffs$June + coeffs$July + coeffs$fixed_slope + 
  1.96 * (coeffs$July_sd + coeffs$fixed_slope_sd + coeffs$June_sd )

# August
derived$August <- coeffs$June + coeffs$August + coeffs$fixed_slope 
derived$August_lwr <- coeffs$June + coeffs$August + coeffs$fixed_slope - 
  1.96 * (coeffs$August_sd + coeffs$fixed_slope_sd + coeffs$June_sd )
derived$August_upr <- coeffs$June + coeffs$August + coeffs$fixed_slope + 
  1.96 * (coeffs$August_sd + coeffs$fixed_slope_sd  + coeffs$June_sd )

# ... Re-combine parameters to look "tidy" ----
# Can actually make it tidy at some point in future if
# can figure out multiple queries in cols for pivot_longer()
# Pivot the table back into long format so we can plot by month
# Means
coeff_ests <- derived %>% 
  dplyr::select(grp, Intercept, June, July, August) %>% 
  pivot_longer(
    cols = c(Intercept, June, July, August),
    values_to = "fit",
    names_to = "Parameter"
  )
# Lowers
coeff_lwr <- derived %>% 
  dplyr::select(grp, Intercept_lwr, June_lwr, July_lwr, August_lwr) %>% 
  pivot_longer(
    cols = c(Intercept_lwr, June_lwr, July_lwr, August_lwr),
    values_to = "lwr",
    names_to = "Parameter"
  )
# Uppers
coeff_upr <- derived %>% 
  dplyr::select(grp, Intercept_upr, June_upr, July_upr, August_upr) %>% 
  pivot_longer(
    cols = c(Intercept_upr, June_upr, July_upr, August_upr),
    values_to = "upr",
    names_to = "Parameter"
  )

# Smash them back together
plotter <- data.frame(coeff_ests, 
           lwr = coeff_lwr$lwr, 
           upr = coeff_upr$upr)

plotter$Parameter <- factor(plotter$Parameter, 
                            levels = c("Intercept", 
                                       "June", 
                                       "July", 
                                       "August"
                                       )
                            )


# Make a data set that we can use for plotting vertical 
# lines for references to mean (Intercept) and zero (slope)
refs <- data.frame(
  Parameter = factor(c("Intercept", "June", "July", "August")),
  xintercept = c(cf[1,1] , 0, 0, 0)
)

# .. Coefficient Plots ----
# Now, we can plot them to show which lakes have 
# higher or lower Secchi than average based on the intercept
# parameters, and which lakes changed based on the slope
# parameters. The lakes with slope parameters that do not
# overlap zero changed across years. Lakes with increases
# are to the right of the dashed lines and those displaying
# decreases are to the left of the dashed line
jpeg(filename = "secchi_coefficients.jpg",
     height = 1800,
     width = 2000,
     res = 300
     )
ggplot(plotter, aes(x = fit, y = grp)) +
  geom_point() +
  geom_segment(aes(x = lwr, xend = upr, yend = grp)) +
  geom_vline(data = refs, aes(xintercept = xintercept), lty = 2) +
  facet_wrap(~ Parameter, scales = "free_x", nrow=1) +
  xlab("Parameter value") +
  ylab("Lake") +
  theme_bw() +
  theme(
    axis.title.x = element_text(vjust = -1),
    axis.title.y = element_text(vjust = 3),
    strip.background = element_blank()
  )    
dev.off()

# . Model predictions ----
# .. All lakes that showed significant increase ----
# Need to get which are significant first
coeffs2 <- plotter %>% 
  filter(Parameter != "Intercept" & (lwr > 0 | upr < 0))

sig_lakes <- lake[lake$lake.name %in% coeffs2$grp, ]


# Summarize raw data with means and 95% CIs
new_d <- sig_lakes %>% 
  group_by(Year, year.s, Month, lake.name) %>% 
  summarize(SECCHI.AVG = mean(SECCHI), SECCHI.SD = sd(SECCHI)) %>% 
  mutate(SECCHI.LWR = SECCHI.AVG - 1.96*SECCHI.SD,
         SECCHI.UPR = SECCHI.AVG + 1.96*SECCHI.SD,
         month = factor(month.name[as.numeric(Month)], month.name)
         ) %>% 
  data.frame()

# Make predictions and get 95% CIs
log_preds <- predictInterval(int.2,
                         newdata = new_d,
                         level = 0.95,
                         type = "linear.prediction",
                         n.sims = 1000,
                         include.resid.var = TRUE)

real_preds <- apply(log_preds, 2, exp)
lake.preds <- data.frame(new_d, real_preds)
  
# Graph predictions on raw data
jpeg(filename = "secchi_updated.jpg",
     height = 2000,
     width = 2000,
     res = 300
     )
ggplot(lake.preds, aes(x = Year, y = SECCHI.AVG))+
  geom_point() +
  geom_segment(aes(xend = Year, y = SECCHI.LWR, yend=SECCHI.UPR)) +
  geom_line(aes(y = fit)) +
  geom_ribbon(aes(xmax = Year, ymin=lwr, ymax = upr,
                  color = NULL), alpha = 0.15) +
  scale_x_continuous(breaks = seq(2005, 2020, 5),
                     labels = seq(2005, 2020, 5),
                     limits = c(2005, 2020)
                     ) +
  xlab("Year") +
  ylab("Secchi depth (ft)") +
  facet_grid(lake.name ~ month) +
  theme_bw() +
  theme(
    axis.title.x = element_text(vjust = -1),
    axis.title.y = element_text(vjust = 3),
    strip.background = element_blank(),
    panel.spacing.x = unit(1, "line")
  )
dev.off()


# Temperature Models ----
# . Data exploration ----
hist(lake_temp$temp)
boxplot(temp ~ Month, data = lake_temp)
boxplot(temp ~ Year, 
        data = lake_temp[lake_temp$Month == "7", ])

boxplot(temp ~ Year + Month, data = lake_temp)
# View(lake_temp)

# . Models ----
# Null Model
null <- lm(log(temp) ~ 1, data = lake_temp)

#Year and Month Additive Model
year.month <- lm(log(temp) ~ year.s + Month, data = lake_temp)

# year.s and Month Interactive Model
year.x.month <- lm(log(temp) ~ year.s * Month, data = lake_temp)


# . Model selection for fixed effects ----
# AIC Ranking (choosing best fixed effects model)
model.list <- list(null, year.month, year.x.month)
model.names <- c("null", "year.month", "year.x.month")

aictab(cand.set = model.list, modnames = model.names)


# . Model selection for random effects ----
# Random effect of lake on intercept only 
int.1 <- lmer(log(temp) ~ year.s * Month + (1 | lake.name), data = lake_temp)
int.1.residuals <- residuals(int.1)
hist(int.1.residuals)

# Random effect of lake on intercept and year slope
int.2 <- lmer(log(temp) ~ year.s * Month + (year.s|lake.name ), data = lake_temp)
lake_temp$resids <- residuals(int.2)
hist(lake_temp$resids)
boxplot(resids ~ Month, data = lake_temp); abline(h = 0)


# AIC Rankings (Choosing best model)

model.list3 <- list(int.1, int.2)
model.names3 <- c("int.1", "int.2")

aictab(cand.set = model.list3, modnames = model.names3)


# . Statistical significance ----
Anova(int.2, type = "III")
summary(int.2)


# . Random effects coeffs ----
# Here is how we get mean and 95% CIs for the random
# intercept and slope adjustments out of this monster.

# .. Get the fixed and random effects ----
cf <- data.frame(summary(int.2)$coefficients) # Fixed
coeffs <- data.frame(ranef(int.2))            # Random

# Pivot the table to unstack parameters and rename
coeffs <- pivot_wider(coeffs, 
                      names_from = term,
                      values_from = c(condval, condsd)
                      )
names(coeffs)[3:6] <- c("Intercept", "June", 
                        "Intercept_sd", "June_sd")

# Get slope adjustments and sds for July and August
# Mean slope adjustment
coeffs$July <- cf[5, 1]
coeffs$July_sd <- cf[5, 2]
coeffs$August <- cf[6, 1]
coeffs$August_sd <- cf[6, 2]

# Add fixed effect coefficients to coeffs table
coeffs$fixed_int <- cf[1, 1]
coeffs$fixed_int_sd <- cf[1, 2]
coeffs$fixed_slope <- cf[2, 1]
coeffs$fixed_slope_sd <- cf[2, 2]

# June is the baseline slope in this model

# July slope needs to be adjusted for (fixed) main effects
# and interaction
# coeffs$july_adj <- coeffs$fixed_slope + coeffs$July + coeffs$June 
# coeffs$july_sd_adj <- coeffs$fixed_slope_sd + coeffs$July_sd + coeffs$June_sd 
# 
# # So does August
# coeffs$august_adj <- coeffs$fixed_slope + coeffs$August + coeffs$June 
# coeffs$august_adj_sd <- coeffs$fixed_slope_sd + coeffs$August_sd + coeffs$June_sd 

# .. Calculate derived parameters in a new df ----
derived <- data.frame(coeffs)

# ... Intercepts and CIs ----
derived$Intercept <- coeffs$Intercept + coeffs$fixed_int
derived$Intercept_lwr <- coeffs$Intercept + coeffs$fixed_int - 
  1.96 * (coeffs$Intercept_sd + coeffs$fixed_int_sd)
derived$Intercept_upr <- coeffs$Intercept + coeffs$fixed_int + 
  1.96 * (coeffs$Intercept_sd + coeffs$fixed_int_sd)

# ... Slopes and CIs ----
# June
derived$June <- coeffs$June + coeffs$fixed_slope 
derived$June_lwr <- coeffs$June + coeffs$fixed_slope - 
  1.96 * (coeffs$June_sd + coeffs$fixed_slope_sd)
derived$June_upr <- coeffs$June + coeffs$fixed_slope + 
  1.96 * (coeffs$June_sd + coeffs$fixed_slope_sd)

# July
derived$July <- coeffs$June + coeffs$July + coeffs$fixed_slope 
derived$July_lwr <- coeffs$June + coeffs$July + coeffs$fixed_slope - 
  1.96 * (coeffs$July_sd + coeffs$fixed_slope_sd + coeffs$June_sd )
derived$July_upr <- coeffs$June + coeffs$July + coeffs$fixed_slope + 
  1.96 * (coeffs$July_sd + coeffs$fixed_slope_sd + coeffs$June_sd )

# August
derived$August <- coeffs$June + coeffs$August + coeffs$fixed_slope 
derived$August_lwr <- coeffs$June + coeffs$August + coeffs$fixed_slope - 
  1.96 * (coeffs$August_sd + coeffs$fixed_slope_sd + coeffs$June_sd )
derived$August_upr <- coeffs$June + coeffs$August + coeffs$fixed_slope + 
  1.96 * (coeffs$August_sd + coeffs$fixed_slope_sd  + coeffs$June_sd )

# ... Re-combine parameters to look "tidy" ----
# Can actually make it tidy at some point in future if
# can figure out multiple queries in cols for pivot_longer()
# Pivot the table back into long format so we can plot by month
# Means
coeff_ests <- derived %>% 
  dplyr::select(grp, Intercept, June, July, August) %>% 
  pivot_longer(
    cols = c(Intercept, June, July, August),
    values_to = "fit",
    names_to = "Parameter"
  )
# Lowers
coeff_lwr <- derived %>% 
  dplyr::select(grp, Intercept_lwr, June_lwr, July_lwr, August_lwr) %>% 
  pivot_longer(
    cols = c(Intercept_lwr, June_lwr, July_lwr, August_lwr),
    values_to = "lwr",
    names_to = "Parameter"
  )
# Uppers
coeff_upr <- derived %>% 
  dplyr::select(grp, Intercept_upr, June_upr, July_upr, August_upr) %>% 
  pivot_longer(
    cols = c(Intercept_upr, June_upr, July_upr, August_upr),
    values_to = "upr",
    names_to = "Parameter"
  )

# Smash them back together
plotter <- data.frame(coeff_ests, 
           lwr = coeff_lwr$lwr, 
           upr = coeff_upr$upr)

plotter$Parameter <- factor(plotter$Parameter, 
                            levels = c("Intercept", 
                                       "June", 
                                       "July", 
                                       "August"
                                       )
                            )


# Make a data set that we can use for plotting vertical 
# lines for references to mean (Intercept) and zero (slope)
refs <- data.frame(
  Parameter = factor(c("Intercept", "June", "July", "August")),
  xintercept = c(cf[1,1] , 0, 0, 0)
)

# .. Coefficient plots ----
# Now, we can plot them to show which lakes have 
# higher or lower Secchi than average based on the intercept
# parameters, and which lakes changed based on the slope
# parameters. The lakes with slope parameters that do not
# overlap zero changed across years. Lakes with increases
# are to the right of the dashed lines and those displaying
# decreases are to the left of the dashed line
jpeg(filename = "temp_coefficients.jpg",
     height = 1800,
     width = 2000,
     res = 300
     )
ggplot(plotter, aes(x = fit, y = grp)) +
  geom_point() +
  geom_segment(aes(x = lwr, xend = upr, yend = grp)) +
  geom_vline(data = refs, aes(xintercept = xintercept), lty = 2) +
  facet_wrap(~ Parameter, scales = "free_x", nrow=1) +
  xlab("Parameter value") +
  ylab("Lake") +
  theme_bw() +
  theme(
    axis.title.x = element_text(vjust = -1),
    axis.title.y = element_text(vjust = 3),
    strip.background = element_blank()
  )    
dev.off()


# . Model predictions ----
# .. All lakes that showed significant increase ----
# Need to get which are significant first
coeffs2 <- plotter %>% 
  filter(Parameter != "Intercept" & (lwr > 0 | upr < 0))

sig_lakes <- lake_temp[lake_temp$lake.name %in% coeffs2$grp, ]


# Summarize raw data with means and 95% CIs
new_d <- sig_lakes %>% 
  group_by(Year, year.s, Month, lake.name) %>% 
  summarize(temp.AVG = mean(temp), temp.SD = sd(temp)) %>% 
  mutate(temp.LWR = temp.AVG - 1.96*temp.SD,
         temp.UPR = temp.AVG + 1.96*temp.SD,
         month = factor(month.name[as.numeric(Month)], month.name)
         ) %>% 
  data.frame()

# Make predictions and get 95% CIs
log_preds <- predictInterval(int.2,
                         newdata = new_d,
                         level = 0.95,
                         type = "linear.prediction",
                         n.sims = 1000,
                         include.resid.var = TRUE)

real_preds <- apply(log_preds, 2, exp)
lake.preds <- data.frame(new_d, real_preds)
  
# Graph predictions on raw data
jpeg(filename = "temps_updated.jpg",
     height = 2500,
     width = 2000,
     res = 300
     )
ggplot(lake.preds, aes(x = Year, y = temp.AVG))+
  geom_point() +
  geom_segment(aes(xend = Year, y = temp.LWR, yend=temp.UPR)) +
  geom_line(aes(y = fit)) +
  geom_ribbon(aes(xmax = Year, ymin=lwr, ymax = upr,
                  color = NULL), alpha = 0.15) +
  scale_x_continuous(breaks = seq(1990, 2020, 5),
                     labels = seq(1990, 2020, 5),
                     limits = c(1990, 2020)
                     ) +
  xlab("Year") +
  ylab(expression(paste("Temperature (", degree, "F)"))) +
  facet_grid(lake.name ~ month) +
  theme_bw() +
  theme(
    axis.title.x = element_text(vjust = -1),
    axis.title.y = element_text(vjust = 3),
    strip.background = element_blank(),
    panel.spacing.x = unit(1, "line")
  )
dev.off()
