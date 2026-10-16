rm(list = ls())
library(truncdist)
library(latex2exp)
library(knitr)
library(nloptr)
library(SigDetect)

### GENERIC PARAMETERS #########################################################
mean_sig <- 125; sd_sig <- 3
u <- 160; l <- 110
k <- 500
bin_ends <- seq(l, u, length.out = k+1)
xi <- (bin_ends[1:k] + bin_ends[2:(k+1)])/2

# Defining signal density and calculating signal region
fs <- function(x) dtrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)

### LOADING DATA ###############################################################
cats <- 1:3
scenarios <- c('HLLHC', 'Run2', 'Run3')
mu_seq <- c(1, 1.2, 1.5, 1.8, 2)
### PROPOSAL BACKGROUND ########################################################
# Defining BW parameters
bkg_loc <- 91.2
bkg_scale <- 2.49/2
h_beta <- function(x, beta){
  exp(beta*x)/(((x-bkg_loc)^2) + bkg_scale^2) 
}
d_log_h <- function(x) x
d2_log_h <- 0

### SIGNAL SEARCH ##############################################################
eta_Z_vec <- eta_true_vec <- mu_vec <- eta_hat_vec <- mu_hat_vec <- scen_vec <- cat_vec <- S_hat_toy_vec <- S_toy_vec <- c()

for(scenario in scenarios){
  for(cat in cats){
    for(mu in mu_seq){
      cat(sprintf('\rWorking on Scenario: %s, Cat: %d, mu: %.1f', scenario, cat, mu))
      
      mu_part_1 <- floor(mu); mu_part_2 <- (10*mu)%%10
      # reading bin counts from Asimov data with mu = 1
      asimov_mu1_file <- paste0('../truthBkgModels/h_Cat', cat, '_',
                                scenario, '_mu1p0_mass.csv')
      ni_asm_mu1 <- read.csv(asimov_mu1_file, header = TRUE)[,2]
      
      # reading bin counts from Asimov data with user defined mu:
      asimov_file <- paste0('../truthBkgModels/h_Cat', cat, '_',
                            scenario,'_mu',mu_part_1,'p',mu_part_2,'_mass.csv')
      ni_asm <- read.csv(asimov_file, header = TRUE)[,2]
      
      # reading bin counts from Asimov data with mu = 0
      asimov_mu0_file <- paste0('../truthBkgModels/h_Cat', cat, '_',
                                scenario, '_mu0p0_mass.csv')
      mi_asm <- read.csv(asimov_mu0_file, header = TRUE)[,2]
      
      # Calculating True injected signal counts, i.e., S_expected
      S_expected <- sum(ni_asm_mu1) - sum(mi_asm)
      
      # Reading bin counts from toy data:
      phys_toy_file <- paste0('../unbinned_data/Cat',cat,'_',
                              scenario,'_mu',mu_part_1,'p',mu_part_2,'.csv')
      phys_data <- read.csv(phys_toy_file, header = FALSE)[,1]
      ni <- hist(phys_data, breaks = bin_ends, plot = FALSE)$counts
      bkg_toy_file <- paste0('../unbinned_data/Cat',cat,'_',
                             scenario,'_mu',0,'p',0,'.csv')
      bkg_data <- read.csv(bkg_toy_file, header = FALSE)[,1]
      mi <- hist(bkg_data, breaks = bin_ends, plot = FALSE)$counts
      
      true_eta <- S_expected*mu/sum(ni_asm)
      
      # Signal search on the toy data
      res_toy_wbkg <- sig_detect_binned_wbkg(phys.counts = ni,
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
      eta_hat_toy <- res_toy_wbkg$eta_hat
      mu_hat_toy <- sum(ni)*eta_hat_toy/S_expected
      S_hat_toy <- sum(ni)*eta_hat_toy
      S_toy <- sum(ni) - sum(mi)
      eta_hat_SE <- res_toy_wbkg$SE
      eta_Z <- (eta_hat_toy - true_eta)/eta_hat_SE
      eta_Z_vec <- c(eta_Z_vec, eta_Z)
      eta_hat_vec <- c(eta_hat_vec, eta_hat_toy)
      eta_true_vec <- c(eta_true_vec, true_eta)
      mu_vec <- c(mu_vec, mu)
      mu_hat_vec <- c(mu_hat_vec, mu_hat_toy)
      scen_vec <- c(scen_vec, scenario)
      cat_vec <- c(cat_vec, cat)
      S_hat_toy_vec <- c(S_hat_toy_vec, S_hat_toy)
      S_toy_vec <- c(S_toy_vec, S_toy)
    }
  }
}
res_df <- data.frame(
  'Scenario' = scen_vec,
  'Cat' = cat_vec,
  'mu' = mu_vec,
  'true_eta' = eta_true_vec,
  'mu_hat' = mu_hat_vec,
  'eta_hat' = eta_hat_vec,
  'eta_Z' = eta_Z_vec,
  'S_hat' = S_hat_toy_vec*(S_hat_toy_vec>0),
  'S_toy' = S_toy_vec
)

hist(res_df$eta_Z, probability = TRUE, breaks = 15, main = 'Histogram of standardized eta hat')
shapiro.test(res_df$eta_Z)
qqnorm(res_df$eta_Z)
qqline(res_df$eta_Z, col = 'red', lwd = 2)

write.csv(res_df, file = 'sig_detect_res_all.csv', row.names = FALSE)
res_df <- read.csv('sig_detect_res_all.csv', header = TRUE)
res_df[(res_df$S_hat <= res_df$S_toy + sqrt(res_df$S_toy)) & (res_df$S_hat >= res_df$S_toy - sqrt(res_df$S_toy)),
       c('Scenario', 'Cat', 'mu', 'mu_hat', 'S_toy', 'S_hat')]
