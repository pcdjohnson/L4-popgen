# Clear memory, close plotting devices
rm(list = ls())
graphics.off()

# Load packages:
library(scales)
library(viridis)


# choose plotting options

# plot a graphic of the genotypes in each sub-population at the final generation?
plot.genotypes <- TRUE
# slow down genetic drift plotting to see one popn at a time
plot.delay <- 0.2 # change to e.g. 0.2 to slow down

# Set up population parameters

# no of diploid individuals
n <- 10
# number of generations to simulate
n.gen <- 150
# frequency of A1 allele at generation 0
p <- 0.5 
# frequency of A2 allele at generation 0
q <- 1 - p
# number of sub-populations to simulate
n.subpop <- 16
# selection coefficient against A1A1 (selection against A1A2 is s/2)
s <- 0 # try e.g. 0.05 
# additional selection against heterozygotes
s.het <- 0 

# probability of fixation for a new allele
1 / (2 * n)
# probability of fixation for a new nearly neutral allele:
-2 * s / (1 - exp(4 * n * s))

# create a template data frame to store genotypes
pop <- expand.grid(id = 1:n, gen = 0:n.gen, allele1 = NA, allele2 = NA)

# initialise founder population with exact allele frequency
pop[pop$gen == 0, c("allele1", "allele2")] <- 
  sample(c(rep("A1", round(p * n * 2)), rep("A2", round((1 - p) * n * 2))))

# have a look
pop[pop$gen == 0, ]

# simulate evolution forwards in time
start.time <- Sys.time()
simres <-
  lapply(1:n.subpop, function(j) {

    print("Generation 0", quote = FALSE)
    print(pop[pop$gen == 0, c("allele1", "allele2")], row.names = FALSE)
    
    for(generation in 1:n.gen) {
      print(paste("Generation", generation), quote = FALSE, row.names = FALSE)
      # pick mate pairs at random
      n.A1 <- rowSums(pop[pop$gen == (generation - 1), c("allele1", "allele2")] == "A1")
      w <- 1 -  n.A1/2 * s
      w.het <- 1 - (n.A1 == 1) * s.het
      parents <- replicate(n, sample(1:n, 2, replace = FALSE, prob = w * w.het))
      this.gen <-
        sapply(1:n, function(i) {
          allele.tab <- pop[pop$gen == (generation - 1), ][parents[, i], c("allele1", "allele2")]
          apply(allele.tab, 1, sample, 1)
        })
      pop[pop$gen == generation, c("allele1", "allele2")] <- t(this.gen)
      print(pop[pop$gen == generation, c("allele1", "allele2")], row.names = FALSE)
    }
    
    freq <- 
      sapply(0:n.gen, function(generation) {
        mean(unlist(pop[pop$gen == generation, c("allele1", "allele2")]) == "A1")
      })
    
    pop$genotype <- apply(pop[, c("allele1", "allele2")], 1, function(x) paste0(sort(x), collapse = ""))
    pop$genotype <- factor(pop$genotype, c("A1A1", "A1A2", "A2A2"))
    frequencies <- data.frame(gen = 0:n.gen)
    frequencies$A1 <- freq
    frequencies$A2 <- 1 - frequencies$A1
    frequencies$A1A1.exp <- frequencies$A1^2
    frequencies$A1A2.exp <- 2 * frequencies$A1 * frequencies$A2
    frequencies$A2A2.exp <- frequencies$A2^2
    frequencies[, c("A1A1.obs", "A1A2.obs", "A2A2.obs")] <- 
      do.call("rbind", (tapply(pop$genotype, pop$gen, function(x) prop.table(table(x)))))
    
    rowSums(frequencies[, c("A1A1.exp", "A1A2.exp", "A2A2.exp")])
    rowSums(frequencies[, c("A1A1.obs", "A1A2.obs", "A2A2.obs")])
    
    #message_parallel(paste0(100 * (j/n.subpop), "%"))
    
    return(list(pop = pop, frequencies = frequencies))
  })
finish.time <- Sys.time()

par(mar = c(5.1, 4.1, 4.1, 5.1))
gap <- ceiling(n.gen/5)
lwd <- 1.7
plot(x = 0:n.gen, y = 0:n.gen, ylim = 0:1, type = "n", xlab = "Generation", 
     ylab = "Freq. of 'A1' allele",
     xlim = c(0, n.gen * 2 + gap), axes = FALSE)
axis.labels <- pretty(c(0:n.gen))
axis(1, at = axis.labels, las = 3, pos = 0)
axis(1, at = axis.labels + n.gen + gap, labels = axis.labels, las = 3, pos = 0)
axis(2, pos = 0)
axis(2, pos = n.gen + gap)
axis(4)
mtext("Freq. of 'A1A2' genotype", side = 4, line = 2.5)
title(paste("Sub-population size =", n))
legend("topright", c("Expected", "Observed"), lty = 1:2,
       col = grey(0.7), bty = "n", lwd = lwd)

cols <- viridis(n.subpop + 1, alpha = 0.7, option = "turbo")

lapply(1:n.subpop, function(i) {
  x <- simres[[i]]
  Sys.sleep(min(plot.delay, 2/n.subpop))
  lines(x = 0:n.gen, y = x$frequencies$A1, 
        col = cols[i], lwd = lwd)
  lines(x = (0:n.gen) + gap + n.gen, y = x$frequencies$A1A2.exp, 
        col = cols[i], lwd = lwd)
  lines(x = (0:n.gen) + gap + n.gen, y = x$frequencies$A1A2.obs, 
        col = cols[i], lty = 2, lwd = lwd)
})

# final frequencies of A1 allele
final.freq <- sapply(simres, function(x) x$frequencies$A1[x$frequencies$gen == n.gen])
final.freq
finish.time - start.time



if(plot.genotypes) {
  old.par <- 
    par(mfrow = c(ceiling(n.subpop/ceiling(sqrt(n.subpop))), 
                  ceiling(sqrt(n.subpop))), 
        mar = c(0.3, 0.3, 0.3, 0.3))
  final.genotypes.list <- 
    lapply(1:n.subpop, function(i) {
      pop <- simres[[i]]$pop
      out <- pop[pop$gen == n.gen, c("allele1", "allele2")]
      plot(1, type = "n", xlim = 0:1, ylim = 0:1, axes = FALSE)
      rect(par("usr")[1], par("usr")[3],
           par("usr")[2], par("usr")[4],
           col = cols[i])
      
      box()
      gap <- 0.014 + 0.01 * sqrt(n.subpop)
      coords <- 
        expand.grid(x = (1:ceiling(n/2))/(ceiling(n/2)+1), 
                    y = (1:ceiling(n/2))/(ceiling(n/2)+1))
      coords <- coords[sample(nrow(coords), n), ] 
      out <- cbind(out, coords)
      apply(out, 1, function(x) {
        xcoord <- jitter(as.numeric(x[3]), factor = 1.5)
        ycoord <- jitter(as.numeric(x[4]), factor = 2)
        points(xcoord, ycoord, pch = 21, cex = 7, bg = "grey20", col = "white")
        alleles <- sort(unlist(x[1:2]))
        cols <- c(A1 = "yellow", A2 = "green")
        text(x = c(xcoord - gap/2, xcoord + gap/2), y = ycoord, labels = alleles, 
             col = cols[alleles], font = 2, cex = 1.1)
      })
      out
    })
  par(old.par)
  final.genotypes <- do.call("rbind.data.frame", final.genotypes.list)[, c("allele1", "allele2")]
  
  
  # Now we can estimate some population genetic parameters, and use them to answer
  # questions relevant to conservation:
  #  1. Are the sub-populations genetically isolated from each other?
  #  2. If so, 
  #     (a) how long have they been isolated for and 
  #     (b) what is the (effective) size of sub-populations?

  # Frequency of the A1 allele across the whole population
  p.est <- mean(unlist(final.genotypes) == "A1")
  print(paste("The frequency of the A1 allele across all", 
              n.subpop, "sub-populations is", round(p.est, 3)), 
        quote = FALSE)
  
  # Expected heterozygosity in the total population
  He.est <- 2 * p.est * (1 - p.est)
  print(paste("The expected heterozygosity across all", 
              n.subpop, "sub-populations is", round(He.est, 3)), 
        quote = FALSE)
  
  # Observed homozygosity
  n.homozygote <- sum(apply(final.genotypes, 1, function(x) length(unique(x))) == 2)
  Ho.est <- n.homozygote / (n * n.subpop)
  print(paste("The observed heterozygosity across all", 
              n.subpop, "sub-populations is", round(Ho.est, 3)), 
        quote = FALSE)
  
  # Calculate expected heterozygosity within each population
  Hs.est <- 
    mean(sapply(final.genotypes.list, function(fg) {
      p.est <- mean(unlist(fg[, c("allele1", "allele2")]) == "A1")
      2 * p.est * (1 - p.est)
    }))
  
  # Now we can calculate an estimate of Fit, which measures the total 
  # heterozygote deficit, potentially due to population structure and/or inbreeding
  
  # If Ho is less than He, then F should be > 0.
  Fit.est <- 1 - Ho.est / He.est
  print(paste("Estimated Fit:", round(Fit.est, 3)), 
        quote = FALSE)
  
  # Here Fit.est means "an estimate of F". F simply gauges the extent of 
  # homozygote excess (how many more homozygotes are there compared to the HWE expectation?).
  
  # We can also calculate Fst (see chapter 9 of CGP, eqn 9.3)
  Fst.est <-  1 - Hs.est / He.est
  print(paste("Estimated Fst:", round(Fst.est, 3)), 
        quote = FALSE)
  

  # We can also rearrange CGP eqn 9.4 to calculate Fis, the inbreeding coefficient,
  # which gauges the portion of heterozygote deficit (Fit) due to inbreeding rather than
  # population structure
  Fis.est <- (Fit.est - Fst.est) / (1 - Fst.est)
  print(paste("Estimated Fis:", round(Fis.est, 3)), 
        quote = FALSE)
  
  
  
  # To summarise, we break Fit down into Fis, the portion of F due to inbreeding 
  # (between individuals within sub-populations) and Fst, the portion of Fit due to divergence 
  # between sub-populations (population structure).  
  # 
  # Why do we expect isolation between sub-populations to lead to fewer homozygotes 
  # than expected and therefore high Fst? Perhaps this is easier to imagine if we take
  # an extreme case. Make the following changes near the top of the code:
  #   Set plot.genotypes to be TRUE, not FALSE
  #   Change the number generations (n.gen) to 100
  # Now run the code by clicking on "Source".

  # You should now see a diagram of all the sub-populations as different coloured 
  # squares, with individuals as black rectangles, and with each individual showing its
  # genotype. Probably most of the sub-populations will be full of either A1A1 homozygotes, 
  # or A2A2 homozygotes, with very few or no heterozygotes. This is because mating 
  # is not random among the whole population but restricted to within sub-populations
  
  
  # We know this assumption is true only because we have simulated the data. 
  # 
  # What does Fst mean, in practical terms? I.e. what can we learn about 
  # the "real" sub-populations from our Fst estimate?
  # To do that we have to make some assumptions about the population,
  # which taken together will form a model.
  # Let's assume that: 
  #   the sub-populations are diverging independently from a common ancestor with no migration
  #   all the sub-populations are the same constant size
  #   no mutation
  #   no selection
  #   the sub-populations have not reached equilibrium (think about why this assumption is necessary)
  #   random mating
  #   non-overlapping generations
  # If all of the above apply, the following equation holds:
  #
  #   Fst = 1 - (1 - 1/(2*n))^t
  # 
  # where n is the size of the sub-populations and t is the number of generations since divergence.
  #
  # My source for this equation was equations 9.9 of CGP.
  #
  # We know both n and t, so we can use this equation to calculate what Fst should be:

  print("***********************************", quote = FALSE)
  print(paste("Estimated Fst:", round(Fst.est, 3)), quote = FALSE)
  print(paste("'True' Fst:", round(1 - (1 - 1/(2*n))^n.gen, 3)), quote = FALSE)
  print("***********************************", quote = FALSE)
  
  # If this was real data, and we knew the (effective) size of each sub-population,
  # we could solve the above equation for t (time in generations), and estimate
  # the age of the sub-populations in no of generations since divergence:
  t.est <- log(1 - Fst.est) / log((1 - 1/(2*n)))
  print("***********************************", quote = FALSE)
  print(paste("Estimated no of generations since splitting:", 
              round(t.est, 1)), 
        quote = FALSE)
  print(paste("True number of generations:", n.gen), 
        quote = FALSE)
  print("***********************************", quote = FALSE)
  #
  # Alternatively, if we know the number of generations since splitting
  # (which is probably more realistic) we can estimate the effective size
  # of the sub-populations.
  # Solve the above equation for n:
  n.est <- 1/(2*(1 - (1 - Fst.est)^(1/n.gen)))
  print("***********************************", quote = FALSE)
  print(paste("Estimated effective sub-population size (Ne):", 
              round(n.est, 1)), 
        quote = FALSE)
  print(paste("True effective sub-population size (Ne):", n), quote = FALSE)
  print("***********************************", quote = FALSE)
  
}






