#contains all recurrent variants that compatible with more than one transcripts
GvsT_recurrent_multiTxs_list = readRDS("/varidata/research/projects/triche/Peter/BamSlicing/CMD_check/bamSliceR/inst/extdata/GvsT_recurrent_multiTxs_list.rds")
#one patient one variant
demo = GvsT_recurrent_multiTxs_list[[16]][1:4,]

plotHistAsDist = function(data, binwidth = 0.001, i)
{
  data = data.frame(d = 1:nrow(data), average = data[,i])
  ggplot(data, aes(x = average)) +
    geom_histogram(binwidth = binwidth, fill = "blue", color = "black", alpha = 0.7) +
    labs(title = "Distribution of True Averages",
         x = "True Average",
         y = "Frequency") +
    theme_minimal()
}
##### apply to our data #####
#Example:V = "chr15:90085321:90085321:T:C" patient:"02H060"
#g_VAF = 0.0967423
#g_altDepth = 98
#g_refDepth = 911
#EC:"ENST00000330062.8", "ENST00000560061.1", "ENST00000559482.5", "ENST00000540499.2"

#model 1: VAF(ENST00000330062.8) == g_VAF(V)
#t_totalDepth(ENST00000330062.8) == 258
#t_altDepth(ENST00000330062.8) == 10
#t_refDepth(ENST00000330062.8) == 248
#1) assume that in model 1 the prior distribution for theta(VAF)
# is Beta(g_altDepth, g_refDepth) == Beta(98, 911) 
# The observed data for "ENST00000330062.8" tells us that there
# are 258 "bernoulli trails"(txs total depth), and 10 "success" (txs alt depth).

#Use simulation to approximate the probability of observing 10 "success" 
#in 258 "bernoulli trails" given that model 1 is correct.
library(ggplot2)
Nrep = 1000000
g_VAF = rbeta(Nrep, 98, 911)
theta_df = data.frame(theta = g_VAF, no = 1:length(g_VAF))
plotHistAsDist(theta_df, 
               binwidth = 0.001, 
               i = which(colnames(theta_df) == "theta"))
y = rbinom(Nrep, 258, g_VAF)
sum(y == 10) / Nrep

## Approximation of integrate over the parameter space for the possible evidence
theta = seq(0,1,0.00001)
#p(theta|M) = p(g_VAF|beta dist) 
prior_of_theta = dbeta(theta, 98, 911)

#likelihood p(X | theta) == p(alt_txs_reads | g_VAF)
likelihood = dbinom(10 , 258, theta)

# integral of p(X| theta, M) * p(theta | M) dtheta
mg_txs1 = sum(likelihood * prior_of_theta * theta)


