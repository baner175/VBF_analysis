rm(list = ls())
library(truncdist)
library(latex2exp)
library(knitr)
library(nloptr)
library(SigDetect)

# source('signal_search_functions.R')

### GENERIC PARAMETERS #########################################################
mean_sig <- 125; sd_sig <- 3
u <- 160; l <- 110
k <- 500
bin_ends <- seq(l, u, length.out = k+1)
xi <- (bin_ends[1:k] + bin_ends[2:(k+1)])/2

# Defining signal density and calculating signal region
fs <- function(x) dtrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)
Fs <- function(x) ptrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)

### LOADING DATA AND CALCULATING S_exp ##########################################
cat <- 0
scenario <- 'HLLHC'
mu <- 1
mu_part_1 <- floor(mu); mu_part_2 <- (10*mu)%%10
B <- 1e4

asimov_phys_file <- paste0('../truthBkgModels/h_Cat', cat, '_',
                           scenario, '_mu', mu_part_1, 'p', mu_part_2, '_mass.csv')
phys_asimov <- read.csv(asimov_phys_file, header = TRUE)[,2]
asimov_bkg_file <- paste0('../truthBkgModels/h_Cat', cat, '_',
                          scenario, '_mu0p0', '_mass.csv')

bkg_asimov <- read.csv(asimov_bkg_file, header = TRUE)[,2]
phys_mu1_asimov <- read.csv(paste0('../truthBkgModels/h_Cat', cat, '_', scenario,
                                   '_mu1p0_mass.csv'),
                            header = TRUE)[,2]
S_exp <- sum(phys_mu1_asimov) - sum(bkg_asimov)

### PROPOSAL BACKGROUND ########################################################
# Defining BW parameters
bkg_loc <- 91.2
bkg_scale <- 2.49/2
h_beta <- function(x, beta){
  exp(beta*x)/(((x-bkg_loc)^2) + bkg_scale^2) 
}
h_alpha <- function(x, alpha){
  exp(alpha*x)/(((x-bkg_loc)^2) + bkg_scale^2) 
}
d_log_h <- function(x) x
d2_log_h <- 0

### SIGNAL SEARCH ##############################################################

S_hat_vec <- mu_hat_vec <- mu_hat_etaf_vec <- c()


start_time <- Sys.time()
set.seed(12345)
for(i in 1:B){
  ni <- sapply(1:k, function(i){
    rpois(1, phys_asimov[i])
  })
  ni <- as.numeric(ni)
  N <- sum(ni)
  
  mi <- sapply(1:k, function(i){
    rpois(1, bkg_asimov[i])
  })
  mi <- as.numeric(mi)
  M <- sum(mi)
  res_wbkg <- sig_detect_binned_wbkg(phys.counts = ni,
                                     bkg.counts = mi,
                                     search_region = c(l,u),
                                     fs = fs,
                                     prop.bkg = h_beta,
                                     d_log_bkg = d_log_h,
                                     d2_log_bkg = d2_log_h,
                                     eta0 = 0,
                                     is.density = FALSE,
                                     beta_lower = -10,
                                     beta_upper = 10,
                                     beta_start = 0.01)
  eta_hat <- res_wbkg$eta_hat
  S_hat <- N*eta_hat
  mu_hat <- S_hat/S_exp
  mu_hat_etaf <- N*mu/sum(phys_asimov)
  S_hat_vec <- c(S_hat_vec, S_hat)
  mu_hat_vec <- c(mu_hat_vec, mu_hat)
  mu_hat_etaf_vec <- c(mu_hat_etaf_vec, mu_hat_etaf)
  cat(sprintf('\rIteration: %d/%d', i,B))
}
end_time <- Sys.time()
print(end_time - start_time)

hist(mu_hat_vec, breaks = 50,
     col = rgb(1, 0, 0, 0.4),
     main = paste0('Distribution of µ_hat: Cat', cat, '-', scenario),
     xlab = "Value")

hist(mu_hat_etaf_vec, breaks = 50,
     col = rgb(0, 0, 1, 0.4),
     add = TRUE)

legend("topright",
       legend = c("µ-hat (eta-estimated)", "µ-hat (eta-fixed)"),
       fill = c(rgb(1, 0, 0, 0.4),
                rgb(0, 0, 1, 0.4)))

abline(v = mu, col = 'red', lty = 2, lwd = 4)
