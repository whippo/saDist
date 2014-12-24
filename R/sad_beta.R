#' Beta Diversity Null Models
#' 
#' This function randomly populates community matrices from various species abundance distributions (SADs),
#' based on input parameters chosen by the user for beta diversity comparison. 
#' 
#' @param spp_count number of species to use in simulation
#' @param Nsamp number of individuals to sample in each plot 
#' @param Nplot number of 'plots' sampled per round
#' @param Nsite number of 'sites' to sample, essentially number of simulations to run
#' @param dist distribution to use for calculating RAD
#' @param method beta measurement used
#' @return Output includes beta diversity plots of all permuted communities. 
#' @author Ross Whippo
#' @examples 
#' #description of examples here
#' #actual R code to execute examples
#' @export 
#' @importFrom VGAM pzipf
#' @importFrom vegan vegdist 
#' 
sad_beta <- function(spp_count, Nsamp, Nplot = 16,  Nsite = 9, dist = "lnorm", method = 'bray')
{
  DISTS <- c("lnorm", "geom", "zipf") 
  dist <- pmatch(dist, DISTS)
  ind <- DISTS[dist]
  if (is.na(dist))
    stop("invalid distribution")
  
  rad.totals <- matrix(0, length(1:Nsite), length(1:(spp_count))) #matrix to hold colsums for each round
  total.comm <- matrix(0, 0, length(1:(spp_count))) #matrix to hold ALL samples
  colnames(total.comm) <- c(1:(spp_count)) #rename columns to match species count
  colnames(rad.totals) <- c(1:(spp_count)) #rename columns to match species count
  beta.mat <- matrix(0, length(1:Nsite), 5) #matrix to hold beta values for each round
  for (i in 1:Nsite) #loop for populating abundance matrix and beta measurements
  {
    comm.mat <- matrix(0, 0, length(1:(spp_count))) #empty matrix to populate with draws
    colnames(comm.mat) <- c(1:(spp_count)) #rename matrix columns for ease of reading
    for(j in 1:Nplot) #loop for plot sampling
    {
      sp <-   (1:(spp_count)) #vector of species observed
      quantile_vect <- seq(0, 1, length.out = (spp_count + 2)) #evenly spaced sequence of numbers 0>n>1
      quantile_vect <- quantile_vect[2:(spp_count + 1)] #evenly spaced draws from distribution
      
      if (dist == 1) #conditional, determine abundance from log normal dist
        abundance <- qlnorm(quantile_vect)
      
      if (dist == 2) #conditional, determine abundance from geometric dist
        abundance <- qgeom(quantile_vect, (1/length(sp)))
      
      if (dist == 3) #conditional, determine abundance from zipf dist
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
    
    METH <- c("bray", "whittaker", "additive", "altGower")
    meth <- pmatch(method, METH)
    inm <- METH[meth]
    if (is.na(method))
      stop("invalid method")
    
    if (meth == 1){
      beta.dist <- vegdist(comm.mat, method = "bray") #calculates bray curtis pairwise distances
      beta.quan <- quantile(beta.dist) #creates list of quartiles for sample
      beta.mat[i,] <- beta.quan
    }
    
    if (meth == 2){
      comm.mat[comm.mat > 0] <- 1
      beta.dist <- apply(comm.mat, 1, sum)
      beta.mean <- mean(beta.dist)  
      beta.value <- (spp_count/beta.mean)
      beta.mat[i,1] <- beta.value
    }
    
    if (meth == 3){
      comm.mat[comm.mat > 0] <- 1
      beta.dist <- apply(comm.mat, 1, sum)
      beta.mean <- mean(beta.dist)  
      beta.value <- (spp_count - beta.mean)
      beta.mat[i,1] <- beta.value
    }
    
    if (meth == 4){
      beta.dist <- vegdist(comm.mat, method = "altGower") #calculates altGower pairwise distances
      beta.quan <- quantile(beta.dist) #creates list of quartiles for sample
      beta.mat[i,] <- beta.quan
    }
  }
  
    if (meth %in% c(1, 4)){
      beta.result <- t(beta.mat)
      colnames(beta.result) <- 1:Nsite
      rownames(beta.result) <- rownames(1:spp_count)
      boxplot(t(beta.mat), xlab = "Site(simulation) #", ylab = c(dQuote(inm),"distance to centroid"))
      } 
  
      if (meth %in% c(2,3)){
        beta.result <- as.matrix(beta.mat[,1])
        plot(beta.result, cex = 2, type = 'p', pch = 15, xlab = "Site(simulation) #", ylab = c(dQuote(inm),"beta"))
        }
    }

 
