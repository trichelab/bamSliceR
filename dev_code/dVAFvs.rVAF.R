# Define observed counts and totals
N_RM <- 2      # Number of RNA mutant reads
N_RT <- 28     # Total RNA reads
N_DM <- 31     # Number of DNA mutant reads (assuming 100 total reads for DNA)
N_DT <- 100    # Total DNA reads

# Define the prior parameters (uniform prior)
alpha <- 1
beta <- 1

# Calculate Beta-Binomial components
BB_RNA <- beta(N_RT + alpha, N_RM + beta) / beta(alpha, beta)
BB_DNA <- beta(N_DT + alpha, N_DM + beta) / beta(alpha, beta)
BB_combined <- beta(N_RT + N_DT + alpha, N_RM + N_DM + beta) / beta(alpha, beta)

# Calculate binomial coefficients (log scale for stability)
log_binom_combined <- lchoose(N_RT + N_DT, N_RM + N_DM)
log_binom_RNA <- lchoose(N_RT, N_RM)
log_binom_DNA <- lchoose(N_DT, N_DM)

# Bayes Factor calculation
BF <- (BB_RNA * BB_DNA / BB_combined) * exp(log_binom_combined - log_binom_RNA - log_binom_DNA)

# Print the Bayes Factor
print(BF)
