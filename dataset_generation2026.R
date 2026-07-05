## Dataset Generation for Simulation 2026 - Mixed Data

library(MASS)     
library(mvtnorm)
library(copula)
library(triangle) 
library(rebmix)
library(GGally)
library(bindata)
library(compositions)

setwd("C:/Users/01459189/OneDrive/phd/hpc - mixed")

## Parameters: 
set.seed(20092025)
n <- 5000
p <- 4  # number of variables

############### SIMULATING CONTINUOUS VARIABLES ########################

# Standard deviations
sd_vec <- c(5, 15, 10, 15)
# Correlation matrix
rho <- 0.5
R <- matrix(rho, nrow=p, ncol=p)
diag(R) <- 1
# Covariance matrix
Sigma <- diag(sd_vec) %*% R %*% diag(sd_vec)

# 1. Multivariate Normal Highly Correlated
mu <- c(15, 100, 75, 175)
#5 → e.g. years of work experience
#20 → e.g. IQ
#75 → e.g. index score or percentage-like variable
#300 → e.g. height

d_norm <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
d_norm[d_norm < 0] <- 0
d_norm <- as.data.frame(d_norm)

# 2 Multivariate lognormal Highly Correlated

# Target means and SDs
mean_vec <- c(50000, 4, 75, 500)
sd_vec   <- c(20000, 1.5, 15, 200)

# Convert to log-normal parameters
sigma_log <- sqrt(log(1 + (sd_vec / mean_vec)^2))
mu_log    <- log(mean_vec) - (sigma_log^2 / 2)

# correlation matrix 
rho <- 0.5
R <- matrix(rho, nrow=p, ncol=p)
diag(R) <- 1

# Covariance matrix on log scale
Sigma_log <- diag(sigma_log) %*% R %*% diag(sigma_log)

# Simulate multivariate log-normal
log_norm_data <- rlnorm.rplus(n,meanlog=mu_log,varlog=Sigma_log)
#log_norm_data <- mvrnorm(n=n, mu=mu_log, Sigma=Sigma_log)
d_lognorm <- as.data.frame(log_norm_data)  # all positive

# 3 Multivariate Normal low correlation
mu <- c(15, 100, 75, 175)

# Standard deviations
sd_vec <- c(5, 15, 10, 15)

# Correlation matrix
rho <- 0.01
R <- matrix(rho, nrow=p, ncol=p)
diag(R) <- 1

# Covariance matrix
Sigma <- diag(sd_vec) %*% R %*% diag(sd_vec)

d_norm_low <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
d_norm_low[d_norm_low < 0] <- 0
d_norm_low <- as.data.frame(d_norm_low)

# 4 Lognormal low correlation

# Target means and SDs
mean_vec <- c(50000, 4, 75, 500)
sd_vec   <- c(20000, 1.5, 15, 200)

# Convert to log-normal parameters
sigma_log <- sqrt(log(1 + (sd_vec / mean_vec)^2))
mu_log    <- log(mean_vec) - (sigma_log^2 / 2)

# correlation matrix 
rho <- 0.01
R <- matrix(rho, nrow=p, ncol=p)
diag(R) <- 1

# Covariance matrix on log scale
Sigma_log <- diag(sigma_log) %*% R %*% diag(sigma_log)

# Simulate multivariate log-normal
log_norm_data <- rlnorm.rplus(n,meanlog=mu_log,varlog=Sigma_log)
d_lognorm_low <- as.data.frame(log_norm_data)  # all positive

# 5. Multivariate Chi-Square
df <- c(2,5,10,4)  # degrees of freedom for each variable
d_chisq <- sapply(df, function(d) rchisq(n, df=d))
colnames(d_chisq) <- paste0("Chi2_", df)

d_chisq <- as.data.frame(d_chisq)


####################################### MIXED DATA #############################

## Binary: 
set.seed(20092025)
p <- 0.33 # probability of 1 ("success")

d_bin <- rbinom(n, size = 1, prob = p)


## Categorical 

d_cat <- factor(
  sample(c("Rural", "Urban", "Other"), 5000, replace = TRUE, prob = c(0.4, 0.66, 0.02)),
  ordered = FALSE,
  levels = c("Rural", "Urban", "Other")
)


## Mixture 1: 1 Normal Continuous and 1 binary
mix1 <- as.data.frame(cbind(d_norm[,1], d_bin))
plot <- ggpairs(mix1[, 1:2], columnLabels = c("Normal", "Binary"))
ggsave("mix1.png", plot = plot, width = 8, height = 6, units = "in", dpi = 300)
save(mix1, file="d_mix1.RData")

## Mixture 2: 1 norm continuous, 1 binary, 1 cat
mix2 <- as.data.frame(cbind(d_norm[,2], d_bin, d_cat))
plot <- ggpairs(mix2[, 1:3], columnLabels = c("Normal", "Binary", "Categorical"))
ggsave("mix2.png", plot = plot, width = 8, height = 6, units = "in", dpi = 300)
save(mix2, file="d_mix2.RData")

## Mixture 3: 1 lognormal continuous, 1 binary 
mix3 <- as.data.frame(cbind(d_lognorm[,1], d_bin))
plot <- ggpairs(mix3[, 1:2], columnLabels = c("Lognormal", "Binary"))
ggsave("mix3.png", plot = plot, width = 8, height = 6, units = "in", dpi = 300)
save(mix3, file="d_mix3.RData")

## Correlated categorical data

Sigma <- matrix(c(1, 0.6, 0.6, 1), 2, 2)  # correlation 0.6
latent <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
x1 <- as.numeric(latent[,1] > qnorm(0.7))
x2 <- as.numeric(latent[,2] > qnorm(0.5))
cor(x1, x2)  # appr cor

# Mixture 5: 2 high correlated continuous and 2 high correlated binary
mix5 <- as.data.frame(cbind(d_norm[,1:2], x1, x2))
plot <- ggpairs(mix5[, 1:4], columnLabels = c("Normal 1", "Normal 2", "Binary 1", "Binary 2"))
ggsave("mix5.png", plot = plot, width = 8, height = 6, units = "in", dpi = 300)
save(mix5, file="d_mix5.RData")

# Mixture 4: 1 continuous, 1 binary, 1 cat
mix4 <- as.data.frame(cbind(d_lognorm[,1], x1, d_cat))
plot <- ggpairs(mix4[, 1:3], columnLabels = c("Lognormal", "Binary", "Categorical"))
ggsave("mix4.png", plot = plot, width = 8, height = 6, units = "in", dpi = 300)
save(mix4, file="d_mix4.RData")

set.seed(20092025)
Sigma <- matrix(c(1, 0.05, 0.05, 1), 2, 2)  # correlation 0.05
latent <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
x1 <- as.numeric(latent[,1] > qnorm(0.7))
x2 <- as.numeric(latent[,2] > qnorm(0.5))
cor(x1, x2)  # appr cor

# Mixture 6: 2 low correlated continuous and 2 low correlated binary

mix6 <- as.data.frame(cbind(d_norm_low[,1:2], x1, x2))
plot <- ggpairs(mix6[, 1:4], columnLabels = c("Normal 1", "Normal 2", "Binary 1", "Binary 2"))
ggsave("mix6.png", plot = plot, width = 8, height = 6, units = "in", dpi = 300)
save(mix6, file="d_mix6.RData")


