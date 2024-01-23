# Clear memory, close plotting devices
rm(list = ls())
graphics.off()

# Load packages:
library(scales)
library(viridis)


# choose plotting options

# plot a graphic of the genotypes in each population at the final generation?
plot.genotypes <- TRUE
# slow down genetic drift plotting to see one popn at a time
plot.delay <- 0.2 # change to e.g. 0.2 to slow down

# Set up population parameters

# no of diploid individuals
n <- 10 
# number of generations to simulate
n.gen <- 15 
# frequency of A1 allele at generation 0
p <- 0.5 
# frequency of A2 allele at generation 0
q <- 1 - p
# number of populations to simulate
n.pop <- 16 
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

# initialise populations with exact allele frequency
pop[pop$gen == 0, c("allele1", "allele2")] <- 
  sample(c(rep("A1", round(p * n * 2)), rep("A2", round((1 - p) * n * 2))))

# have a look
pop[pop$gen == 0, ]

# simulate evolution forwards in time
start.time <- Sys.time()
simres <-
  lapply(1:n.pop, function(j) {

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
    
    #message_parallel(paste0(100 * (j/n.pop), "%"))
    
    return(list(pop = pop, frequencies = frequencies))
  })
finish.time <- Sys.time()

par(mar = c(5.1, 4.1, 4.1, 5.1))
gap <- ceiling(n.gen/5)
lwd <- 1.7
plot(x = 0:n.gen, y = 0:n.gen, ylim = 0:1, type = "n", xlab = "Generation", 
     ylab = "Freq. of 'A1' allele",
     xlim = c(0, n.gen * 2 + gap), axes = FALSE)
#box()
#abline(v = n.gen + gap)
#abline(v = n.gen)
axis.labels <- pretty(c(0:n.gen))
axis(1, at = axis.labels, las = 3, pos = 0)
axis(1, at = axis.labels + n.gen + gap, labels = axis.labels, las = 3, pos = 0)
axis(2, pos = 0)
axis(2, pos = n.gen + gap)
axis(4)
mtext("Freq. of 'A1A2' genotype", side = 4, line = 2.5)
title(paste("Population size =", n))
legend("topright", c("Expected", "Observed"), lty = 1:2,
       col = grey(0.7), bty = "n", lwd = lwd)

cols <- viridis(n.pop + 1, alpha = 0.7, option = "turbo")

lapply(1:n.pop, function(i) {
  x <- simres[[i]]
  Sys.sleep(min(plot.delay, 2/n.pop))
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
#hist(final.freq[final.freq > 0 & final.freq < 1])
finish.time - start.time



if(plot.genotypes) {
  old.par <- par(mfrow = c(ceiling(sqrt(n.pop)), ceiling(sqrt(n.pop))), mar = c(0.3, 0.3, 0.3, 0.3))
  final.genotypes.list <- 
    lapply(1:n.pop, function(i) {
      pop <- simres[[i]]$pop
      out <- pop[pop$gen == n.gen, c("allele1", "allele2")]
      plot(1, type = "n", xlim = 0:1, ylim = 0:1, axes = FALSE)
      rect(par("usr")[1], par("usr")[3],
           par("usr")[2], par("usr")[4],
           col = cols[i]) # Color
      
      box()
      gap <- 0.05 * sqrt(n.pop)/3
      coords <- 
        expand.grid(x = (1:ceiling(n/2))/(ceiling(n/2)+1), 
                    y = (1:ceiling(n/2))/(ceiling(n/2)+1))
      coords <- coords[sample(nrow(coords), n), ] 
      out <- cbind(out, coords)
      apply(out, 1, function(x) {
        xcoord <- jitter(as.numeric(x[3]), factor = 1.5)
        ycoord <- jitter(as.numeric(x[4]), factor = 2)
        points(c(xcoord, xcoord + gap), c(ycoord, ycoord), pch = 15, cex = 4, 
               col = "grey20")
        alleles <- sort(unlist(x[1:2]))
        cols <- c(A1 = "yellow", A2 = "green")
        text(x = c(xcoord, xcoord + gap), y = ycoord, labels = alleles, col = cols[alleles],
             font = 2, cex = 1.1)
      })
      out
    })
  par(old.par)
  final.genotypes <- do.call("rbind.data.frame", final.genotypes.list)[, c("allele1", "allele2")]
  
  
  # now we can estimate some ppopulation genetic parameters
  
  # frequency of the p allele
  p.est <- mean(unlist(final.genotypes) == "A1")
  
  # expected heterozygosity
  He <- 2 * p.est * (1 - p.est)
  
  # observed homozygosity
  n.homozygote <- sum(apply(final.genotypes, 1, function(x) length(unique(x))) == 2)
  Ho <- n.homozygote / (n * n.pop)
  
  # now we can calculate an estimate of F, the inbreeding coefficient
  F.est <- 1 - Ho / He
  
  # Here F.est means "an estimate of F". F simply gauges the extent of 
  # homozygote excess (how many more homozygotes are there compared to the HWE expectation?).
  # Normally we would break F down into Fis, the portion of F due to inbreeding 
  # (between individuals within populations) and Fst, the portion of F due to divergence 
  # between populations (population structure). For the sake of simplicity, here I am 
  # assuming no inbreeding (Fis = 0) and assuming all of F is due to population structure
  # (F = Fst). We know this assumption is true only because we have simulated the data. 
  # *** You would never normally assume F = Fst, I'm only doing this to simplify this code ***
  # 
  # What does Fst mean, in practical terms? I.e. what can we learn about 
  # the "real" populations from our Fst estimate?
  # To do that we have to make some assumptions about the populations,
  # which taken together will form a model.
  # Let's assume that: 
  #   the populations are diverging independently from a common ancestor with no migration
  #   all the populations are the same constant size
  #   no mutation
  #   no selection
  #   the populations have not reached equilibrium (think about why this assumption is necessary)
  #   random mating
  #   non-overlapping generations
  # If all of the above apply, the following equation holds:
  #
  #   Fst = 1 - (1 - 1/(2*n))^t
  # 
  # where n is the size of the populations and t is the number of generations since divergence.
  #
  # (My source for this equation was equations s1 and s2 of this paper:
  # https://doi.org/10.1101%2Fgr.154831.113
  # the original sources are older and are cited above equation s1.)
  #
  # We know both n and t, so we can use this equation to calculate what F (which here = Fst)
  # should be:

  print("***********************************")
  print(paste("Estimated F:", round(F.est, 3)))
  print(paste("'True' F:", round(1 - (1 - 1/(2*n))^n.gen, 3)))
  print("***********************************")
  
  # If this was real data, and we knew the (effective) size of each population,
  # we could solve the above equation for t (time in generations), and estimate
  # the age of the populations in no of generations since divergence:
  t.est <- log(1 - F.est) / log((1 - 1/(2*n)))
  print("***********************************")
  print(paste("Estimated no of generations since splitting:", round(t.est, 1)))
  print(paste("True number of generations:", n.gen))
  print("***********************************")
  #
  # Alternatively, if we know the number of generations since splitting
  # (which is probably more realistic) we can estimate the effective size
  # of each population.
  # Solve the above equation for n:
  n.est <- 1/(2*(1 - (1 - F.est)^(1/n.gen)))
  print("***********************************")
  print(paste("Estimated effective population size (Ne):", round(n.est, 1)))
  print(paste("True effective population size (Ne):", n))
  print("***********************************")
  
}






