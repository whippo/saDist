#' RAD Null Models
#' 
#' This function randomly populates community matrices from various species abundance distributions (SADs),
#' based on input parameters chosen by the user. 
#' 
#' @param spp_count number of species to use in simulation
#' @param Nsamp number of individuals to sample in each plot 
#' @param Nplot number of 'plots' sampled per round
#' @param Nsite number of 'sites' to sample, essentially number of simulations to run
#' @param dist distribution to use for calculating RAD
#' @param gamma_step integer for stepwise increase of gamma with each progressive site sampled
#' @return Uses linear model to fit line to rank abundance distribution, plots and provides slope. 
#' @author Ross Whippo
#' @examples 
#' #description of examples here
#' #actual R code to execute examples
#' @export 
#' @importFrom VGAM pzipf
#' 
sad_beta <- function(spp_count, Nsamp, Nplot = 16,  Nsite = 9, dist = "qlnorm", gamma_step = 0)
{
  DISTS <- c("qlnorm", "qzipf") 
  dist <- pmatch(dist, DISTS)
  ind <- DISTS[dist]
  if (is.na(dist))
    stop("invalid distribution")
  
  rad.totals <- matrix(0, length(1:Nsite), length(1:(spp_count + (gamma_step * Nsite)))) #matrix to hold colsums for each round
  total.comm <- matrix(0, 0, length(1:(spp_count + (gamma_step * Nsite)))) #matrix to hold ALL samples
  colnames(total.comm) <- c(1:(spp_count + (gamma_step * Nsite))) #rename columns to match species count
  colnames(rad.totals) <- c(1:(spp_count + (gamma_step * Nsite))) #rename columns to match species count
  #beta.mat <- matrix(0, length(1:Nsite), 5) #matrix to hold beta values for each round
  for (i in 1:Nsite) #loop for populating abundance matrix and beta measurements
  {
    comm.mat <- matrix(0, 0, length(1:(spp_count + (gamma_step * Nsite)))) #empty matrix to populate with draws
    colnames(comm.mat) <- c(1:(spp_count + (gamma_step * Nsite))) #rename matrix columns for ease of reading
    for(j in 1:Nplot) #loop for plot sampling
    {
      sp <-   (1:(spp_count + ((gamma_step * i) - gamma_step))) #vector of species observed
      quantile_vect <- seq(0, 1, length.out = (spp_count + ((gamma_step * i) - gamma_step) + 2)) #evenly spaced sequence of numbers 0>n>1
      quantile_vect <- quantile_vect[2:(spp_count + ((gamma_step * i) - gamma_step) +1 )] #evenly spaced draws from distribution
      
      if (dist == 1) #conditional, determine abundance from log normal dist
        abundance <- qlnorm(quantile_vect)
      
      if (dist == 2) #conditional, determine abundance from zipf dist
        abundance <- qzipf(quantile_vect, length(sp), 1)
      
      probs <- abundance/sum(abundance) #converts abundance to probability
      probs <- rev(probs) #reverse order of probs to make sp 1 > sp 2 > sp 3....etc.
      comm <- sample(sp, size = Nsamp, replace = TRUE, prob = probs) #sample community
      comm <- t(as.matrix(table(comm))) #turn draws into matrix
      samecols <- intersect(colnames(comm.mat), colnames(comm)) #vector to align column names between matrices
      comm.mat <- merge(comm.mat,comm, by = samecols, all =TRUE) #put draws into matrix 
    }
    comm.mat <- comm.mat[,sort.list(as.numeric(colnames(comm.mat)))] #reorders columns by species ID
    comm.mat[is.na(comm.mat)] <- 0 #replace NAs with 0
    comm.cols <- colSums(comm.mat) #sum species totals for site[i]
    comm.cols <- t(as.matrix(comm.cols)) #turn sums into matrix
    rad.totals[i,] <- comm.cols #add summed species abund totals to site matrix
  }
  colnames(rad.totals) <- c(as.numeric(1:(spp_count + (gamma_step * Nsite)))) #rename matrix columns for ease of reading
  rad.quantiles <- apply(rad.totals, 2, quantile)
  rad.means <- apply(rad.totals, 2, mean)
  new.rad <- rbind(rad.means, sp)
  plot(log10(new.rad[1, ]) ~ new.rad[2,], xlab = "spp", ylab = "log abundance")
  rad.lm <- lm(log10(new.rad[1, ]) ~ new.rad[2,])
  lm.summ <- summary(rad.lm)
  abline(rad.lm)
  
  results <- lm.summ
  return(results)
  
} 
