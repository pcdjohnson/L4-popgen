##################################################################################################
# The aim of this simulation is to give you a sense for how the complex trait model works,       #
# that is, how the small effects of a large number of discrete factors such as genetic variants  #
# - plus the effect of the environment - can combine to produce smooth, normally distributed     #
# quantitative variation.                                                                        #
##################################################################################################

# This script simulates data similar to that in Figure 11.1 of
# Conservation and the Genomics of Populations (3rd edn), 
# https://doi.org/10.1093/oso/9780198856566.003.0011
# "Simulated distributions of phenotypes resulting from one, four, and 10 loci, 
# all with VA = 1 and VE = 0.05 showing how relatively few loci and a small 
# amount of environmental variance can generate an approximately normal distribution
# of phenotypes within a population. From Coop (2019)."


# Clear objects from memory
rm(list = ls())

# Choose the characteristics of the simulated data:

# No of individuals
n <- 1000

# No of SNP loci
n.loci <- 10

# Allele frequencies are drawn from a uniform distribution
q <- runif(n.loci)
p <- 1 - q
p

# Now we choose the effect of each locus on the trait

# Assume that the loci have different effects on the trait, could be 
# positive or negative
loci.effects <- rnorm(n.loci, sd = sqrt(3/n.loci))
loci.effects

# I've set up the simulation to give additive genetic variance Va = 1,
# on average (although the actual Va could be randomly quite different 
# for very small number of loci - try for yourself).

# Set "expected Va":
Va.exp  <- 1

# Now choose the environmental variance
Ve <- 0.05

# Trait (narrow sense) heritability, Hn, is therefore 
# expected to be
Hn.exp <- Va.exp / (Va.exp + Ve)
Hn.exp

# ...we can check whether Hn matches this expectation once we've simulated the data

# Simulate the genotypes as the number of "effect" alleles carried by each individual,
# so every column is a locus, with values of 0, 1 or 2 effect alleles:
dat <- 
  data.frame(sapply(1:n.loci, function(i) {
    geno.freq <- c(p[i]^2, 2*p[i]*q[i], q[i]^2) # genotype frequencies assuming HWE
    sample(0:2, n, replace = TRUE, prob = geno.freq)
  }))
head(dat)


# We can now calculate the effect of the loci on the trait for each individual
# by multiplying the number of effect alleles carried...
# (e.g. here's the first individual)
dat[1, ]
# ...by the effect for each locus
dat[1, ] * loci.effects
# and then summing across loci
sum(dat[1, ] * loci.effects)

# ...this was just for illustration, the following line does this for all individuals:
dat$genetic.effect <- (as.matrix(dat) %*% loci.effects)[, 1]

# We now have a column of genetic effects for each individual
dat$genetic.effect
head(dat)

# Check that the effect on the trait for individual 1 is the same as the value we 
# just calculated:
sum(dat[1, 1:n.loci] * loci.effects)
dat$genetic.effect[1]

# Add environmental variation to the data frame, using the value for Ve we set above:
dat$environment.effect <- rnorm(n, sd = sqrt(Ve))
head(dat)

# Sum the genetic and environmental effects to get the trait value
dat$trait <- dat$genetic.effect + dat$environment.effect
head(dat)

# Look at the trait distribution. Can you produce something similar to Figure 11.1?
hist(dat$trait, breaks = max(10, n/20))

# calculate Va, Ve, and Hn
Va <- var(dat$genetic.effect)
Ve <- var(dat$environment.effect)
Vp <- Va + Ve
Hn <- Va/Vp
Hn

# Check that the expected Hn is the close to realised Hn (for small numbers of loci, 
# n.loci < 10, they could be quite different)
Hn.exp
Hn
p
