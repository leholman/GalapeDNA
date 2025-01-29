## DEMO of asymetric matrix problem
## Luke E Holamn 22Jan2025

#libraries
library("corMLPE")
library("nlme")

#
data_corMLPE <- read.csv("data.csv")

# Ensure 'Start' and 'End' are factors
data_corMLPE$Start <- as.factor(data_corMLPE$Start)
data_corMLPE$End <- as.factor(data_corMLPE$End)

# Scale predictors
data_corMLPE$GeographicDistance_scaled <- scale(data_corMLPE$GeographicDistance)
data_corMLPE$TempDistance_scaled <- scale(data_corMLPE$TempDistance)
data_corMLPE$OceanResistance_scaled <- scale(data_corMLPE$OceanResistance)

# Define the corMLPE correlation structure without a grouping factor
cor_mlpe <- corMLPE(form = ~ Start + End)


# Fit the model using gls
model_corMLPE <- gls(
  eDNA_Distance ~ GeographicDistance + TempDistance + OceanResistance,
  data = data_corMLPE,
  correlation = corMLPE(form = ~ Start + End),
  method = "REML"
)

# Fit the model using gls (scaled)
model_corMLPE_scaled <- gls(
  eDNA_Distance ~ GeographicDistance_scaled + TempDistance_scaled + OceanResistance_scaled,
  data = data_corMLPE,
  correlation = corMLPE(form = ~ Start + End),
  method = "REML"
)

## output

summary(model_corMLPE)

summary(model_corMLPE_scaled)

library(nlme)


testfit <-lmer(eDNA_Distance ~ GeographicDistance_scaled + TempDistance_scaled + OceanResistance_scaled + (1|Start) + (1|End) , data=data_corMLPE)

testfit <-lmer(eDNA_Distance ~ GeographicDistance + TempDistance + OceanResistance + + (1|Start) + (1|End)  , data=data_corMLPE)


summary(testfit)


library(nlme)

testfit <- lme(
  eDNA_Distance ~ GeographicDistance + TempDistance + OceanResistance, 
  random = ~ 1 | Start + End, 
  correlation = corSymm(form = ~ 1 | Pair),
  data = data_corMLPE
)


# Convert 'Start' and 'End' to factors
data_corMLPE$Start <- as.factor(data_corMLPE$Start)
data_corMLPE$End <- as.factor(data_corMLPE$End)

# Fit the lme model
testfit <- lme(
  eDNA_Distance ~ GeographicDistance + TempDistance + OceanResistance, 
  random = ~ 1 | Start/End, 
  data = data_corMLPE
)

testfit <- lme(
  eDNA_Distance ~ GeographicDistance + TempDistance + OceanResistance, 
  random = ~ 1 | Start + End, 
  data = data_corMLPE
)


summary()

install.packages("effects")  # If not already installed
library(effects)

# Compute partial effects
partial_effects <- allEffects(testfit)

# Plot partial effects
plot(partial_effects)


## informamtion from Nate


```
# fake data
set.seed(1)
intercept <- 5; slope <- 1
row_sd <- 1.5; col_sd <- 0.5; obs_sd <- 1.0
n_pop <- 100

covariate <- matrix(rexp(n_pop^2), n_pop, n_pop)
covariate <- (covariate + t(covariate))/2 #symmetric covariate
row_eff <- rnorm(n_pop, 0, row_sd)
col_eff <- rnorm(n_pop, 0, col_sd)
obs_eff <- matrix(rnorm(n_pop^2, 0, obs_sd), n_pop, n_pop)
response <- intercept + slope*covariate + outer(row_eff, col_eff, "+") + obs_eff
rownames(response) <- rownames(covariate) <- paste0("pop", 1:n_pop)
colnames(response) <- colnames(covariate) <- paste0("pop", 1:n_pop)

df <- merge(
  as.data.frame.table(response, responseName="response"),
  as.data.frame.table(covariate, responseName="covariate"),
)
colnames(df)[1:2] <- c("row", "col")
df <- df[df$row != df$col, ] # no "self-comparisons"
head(df)

# fit dyadic model
library(lme4)
fit <- lmer(response ~ covariate + (1|row) + (1|col), data=df)
print(summary(fit))

summary(fit)
aov(fit)


```
