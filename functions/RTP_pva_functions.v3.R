###################################################################################################################################################
##
## Population Projection Function
## This function runs the simulation model, and requires the xx input parameters as arguments (each has a default value, see Table 1 in the paper)
##
####################################################################################################################################################


# Define the population based simulation function
RTP.projection.func <- function(nYears = 80, #years of sim. 80 'years' = 20 years as each 'year' is a 3-month period
                                initialClass = 4, #age class of initial release i.e. 2:5
                                nInitial = 100, #initial release size
                                releasedSexRatio = 0.5, #initial sex ratio of 1st release
                                translocatedSexRatio = 0.5, #sex ratio for all subsequent translocations
                                survivalRates = list(
                                  maleDepJuvenile = 0.6, maleJuvenile = 0.6, maleSubadult = 0.5, maleAdult = 0.4,
                                  femaleDepJuvenile = 0.6, femaleJuvenile = 0.7, femaleSubadult = 0.6, femaleAdult = 0.5
                                ), #survival rates for different classes
                                sexRatio = 0.62, fecundityRateMean = 7.5, fecundityRateSD = 0.69, maxOffspring = 8, #reproductive parameters for RTP
                                breedingSuccess = 1, #probability of females breeding. i.e. 1 = all females breed
                                environmentStochasticity = 0.1,
                                carryingCapacity = 500, 
                                numSimulations = 1, 
                                translocations = list(
                                  list(tYear = 4, numIndividuals = 50),
                                  list(tYear = 8, numIndividuals = 50))
                                ) {
  
  # Set up a data frame to store results. initialise all as zero, so no NA's if population goes extinct
  res.mat <- matrix(NA,nr=numSimulations,ncol=31)
  colnames(res.mat) <- c('nInitial', 'initialClass', 'releasedSexRatio', 
                         'sexRatio', 'fecundityRateMean', 'fecundityRateSD', 'maxOffspring',
                         'breedingSuccess', 
                         'maleDepJuvenileSR', 'maleJuvenileSR', 'maleSubadultSR', 'maleAdultSR', 'femaleDepJuvenileSR', 'femaleJuvenileSR', 'femaleSubadultSR', 'femaleAdultSR', 
                         'translocatedSexRatio', 'tYear1', 'nYear1', 'tYear2', 'nYear2',
                         'environmentStochasticity', 'carryingCapacity', 
                         'nYears', 'numSimulations',
                         'Nfinal', 'Extant', 'Nmin', 'N5yr', 'lambda', 'r')
  
  #results for each year (for function testing)
  results <- data.frame(Generation = integer(nYears * numSimulations),
                        PopulationSize = integer(nYears * numSimulations),
                        Simulation = integer(nYears * numSimulations),
                        season = integer(nYears * numSimulations))
  
  # 2. Run simulation for each replicate (each iteration)
  for (sim in 1:numSimulations) {
    # 2.1 Initialize variables for the simulation
    males <- rep(0, 5)  # Males in age classes: 1 = dependent, 2 = independent, 3 = subadult, 4 = adult
    females <- rep(0, 13) # Females in age classes: 1 = dependent, 2 = independent, 3 = subadult, 4 = adult (4+ years), 5-13 = old adults (5-13 years)
    males[initialClass] <- floor(nInitial * releasedSexRatio)  # adjust by releasedSexRatio
    females[initialClass] <- nInitial - males[initialClass] 
    populationHistory <- numeric(nYears)
    if (initialClass == 5){ #if initial class = 5, i.e. pregnant females need to add in juveniles as well, and for it to be season 1
      males <- males[1:4] #delete males
      females[5] <- nInitial
      pouchYoung <- round(rnorm(nInitial, mean = fecundityRateMean, sd = fecundityRateSD))
      pouchYoung <- pmin(pouchYoung, maxOffspring)
      pouchYoungM <-  sapply(pouchYoung, function(x) sum(rbinom(x, 1, sexRatio)))
      pouchYoungF <- pouchYoung - pouchYoungM
      males[1] <- males[1] + sum(pouchYoungM)
      females[1] <- females[1] + sum(pouchYoungF)
    }
    
    # Initialize the season vector (1 = sep-oct-nov, 2 = dec-jan-feb, 3 = mar-apr-may, 4 = june-july-aug)
    season <- rep(1:4, length.out = nYears)  # Repeat the sequence 1, 2, 3, 4 for the number of years
    season <- ((1:nYears) + (initialClass - 1)) %% 4    
    season[season == 0] <- 4  # Ensure 0 maps to 4
    
    # initialise offspring vector
    offspringMales <- 0
    offspringFemales <- 0
    
    # 2.2 Run the simulation for each year
    for (year in 1:nYears) {
      
      # 4. Breeding: Only adult females (4 years and older) breed
      if (males[4]>0 & (females[4]>0 | females[8]>0 | females[12]>0)) {  # Check if there are any adult females and adult males to breed
        breedingFemales <- round((females[4]+females[8]+females[12]) * breedingSuccess)  # Fraction of adult females breeding
        if (breedingFemales > 0){ #need this conditional in here if breedingSuccess is low, breedingFemales can round to zero
        totalOffspring <- round(rnorm(breedingFemales, mean = fecundityRateMean, sd = fecundityRateSD)) # Generate the total number of offspring (considering both male and female offspring)
        totalOffspring <- pmin(totalOffspring, maxOffspring)  # Ensure max offspring
        malesBorn <- sapply(totalOffspring, function(x) sum(rbinom(x, 1, sexRatio)))# Vectorized approach to generate number of males in total offspring using Bernoulli trials - 0.62 chance of being male
        femalesBorn <- totalOffspring - malesBorn  # Remaining offspring are female
        
        # Assign offspring numbers
        offspringMales <- malesBorn
        offspringFemales <- femalesBorn
        }
      }
      
      # 5. Mortality
      
      #Apply survival to the respective age classes with environmental stochasticity pmax(pmin(..., 1), 0) to clamp values between 0 and 1
      offspringFemales <- round(sum(offspringFemales) * pmax(pmin(rnorm(1, mean = survivalRates$femaleDepJuvenile, sd = environmentStochasticity), 1), 0))
      females[1] <- round(females[1] * pmax(pmin(rnorm(1, mean = survivalRates$femaleDepJuvenile, sd = environmentStochasticity), 1), 0))
      females[2] <- round(females[2] * pmax(pmin(rnorm(1, mean = survivalRates$femaleJuvenile, sd = environmentStochasticity), 1), 0))
      females[3] <- round(females[3] * pmax(pmin(rnorm(1, mean = survivalRates$femaleSubadult, sd = environmentStochasticity), 1), 0))
      females[4:13] <- round(females[4:13] * pmax(pmin(rnorm(1, mean = survivalRates$femaleAdult, sd = environmentStochasticity), 1), 0))  # All adult females (4-13) have the same survival rate
      
      offspringMales <- round(sum(offspringMales) * pmax(pmin(rnorm(1, mean = survivalRates$maleDepJuvenile, sd = environmentStochasticity), 1), 0))
      males[1] <- round(males[1] * pmax(pmin(rnorm(1, mean = survivalRates$maleDepJuvenile, sd = environmentStochasticity), 1), 0))
      males[2] <- round(males[2] * pmax(pmin(rnorm(1, mean = survivalRates$maleJuvenile, sd = environmentStochasticity), 1), 0))
      males[3] <- round(males[3] * pmax(pmin(rnorm(1, mean = survivalRates$maleSubadult, sd = environmentStochasticity), 1), 0))
      males[4] <- round(males[4] * pmax(pmin(rnorm(1, mean = survivalRates$maleAdult, sd = environmentStochasticity), 1), 0))
      
      # 6. Aging: Age all individuals by one year
      males <- c(0, males[1:4])  # Shift males to the next age class
      males <- males[1:4] #drop overage individuals/extra age classes i.e. males[5]
      females <- c(0, females[1:13])  # Shift females to the next age class (including aging all adults 4+)
      females <- females[1:13] #drop overage individuals/extra age classes i.e. females[14]
      
      # Add offspring to the dependent age classes (juveniles)
      males[1] <- males[1] + offspringMales  # Add the number of male offspring to the juvenile age class
      females[1] <- females[1] + offspringFemales  # Add the number of female offspring to the juvenile age class
      
      #after addition to population vector, revert offspring to zero
      offspringFemales <- 0
      offspringMales <- 0
      
      # 8. Translocations: Add individuals to the population from external sources
      if (length(translocations) > 0) {
        # Look for any translocations at this year, if year in any element of the list == year, then that translocation list is selected
        translocationEvent <- translocations[which(sapply(translocations, function(x) x$tYear == year))]
        # If there is a translocation event for this year, add individuals
        if (length(translocationEvent) > 0) {
          translocationEvent <- translocationEvent[[1]]  # Get the first translocation event for this year
          numToAdd <- translocationEvent$numIndividuals
          if (numToAdd > 0) {
            # Add individuals to the age class appropriate for the current season
            if (season[year] == 1) {
              males[2] <- males[2] + round(numToAdd * translocatedSexRatio)
              females[2] <- females[2] + round(numToAdd * (1 - translocatedSexRatio))
            } else if (season[year] == 2) {
              males[3] <- males[3] + round(numToAdd * translocatedSexRatio)
              females[3] <- females[3] + round(numToAdd * (1 - translocatedSexRatio))
            } else if (season[year] == 3) {
              males[4] <- males[4] + round(numToAdd * translocatedSexRatio)
              females[4] <- females[4] + round(numToAdd * (1 - translocatedSexRatio))
            } else if (season[year] == 4) {
              females[5] <- females[5] + round(numToAdd * translocatedSexRatio)
              pouchYoung <- round(rnorm(numToAdd, mean = fecundityRateMean, sd = fecundityRateSD))
              pouchYoung <- pmin(pouchYoung, maxOffspring)
              pouchYoungM <-  sapply(pouchYoung, function(x) sum(rbinom(x, 1, sexRatio)))
              pouchYoungF <- pouchYoung - pouchYoungM
              males[1] <- males[1] + sum(pouchYoungM)
              females[1] <- females[1] + sum(pouchYoungF)
            }
          }
        }
      }
      
      # 9. Carrying Capacity: If the population exceeds capacity, reduce it via probabilistic truncation, as in Vortex
      #truncation in this way is probably realistic for phascogales, which are territorial. This basically simulates density dependent dispersal outside of the study area.
      totalPopulation <- sum(males) + sum(females)
      if (totalPopulation > carryingCapacity) {
        excess <- totalPopulation - carryingCapacity
        adjustmentFactor <- carryingCapacity / totalPopulation  # 1 / (1 + (excess / totalPopulation))
        
        # Apply the survival adjustment to all age classes (i.e., probabilistic truncation)
        # Survival probability is adjusted for age classes individually
        males[1:3] <- round(males[1:3] * adjustmentFactor)
        females[1:3] <- round(females[1:3] * adjustmentFactor)
        males[4] <- round(males[4] * adjustmentFactor)
        females[4:13] <- round(females[4:13] * adjustmentFactor)
      }
      
      # Calculate the total population size
      populationSize <- sum(males) + sum(females)
      populationHistory[year] <- populationSize
      
      # Store results of population for each year, for function testing
      results$Generation[(sim - 1) * nYears + year] <- year
      results$PopulationSize[(sim - 1) * nYears + year] <- populationSize
      results$Simulation[(sim - 1) * nYears + year] <- sim
      results$season[(sim - 1) * nYears + year] <- season[year]
      
      
      # Stop the simulation if the population goes extinct
      if (sum(males) == 0 | sum(females) == 0) {
        break
      }
    }
    
    ## calculate mean lambda and r
    (lambdas <- populationHistory[2:(nYears+1)]/populationHistory[1:nYears])
    (lambda <- mean(lambdas,na.rm=T))
    (r <- log(lambda))
    if(r%in%c('-Inf','Inf')) {r <- NA}
    
    ## store summary outputs of each simulation run
    res.mat[sim, 'nInitial'] <- nInitial
    res.mat[sim, 'initialClass'] <- initialClass
    res.mat[sim, 'releasedSexRatio'] <- releasedSexRatio
    res.mat[sim, 'sexRatio'] <- sexRatio
    res.mat[sim, 'fecundityRateMean'] <- fecundityRateMean
    res.mat[sim, 'fecundityRateSD'] <- fecundityRateSD
    res.mat[sim, 'maxOffspring'] <- maxOffspring
    res.mat[sim, 'breedingSuccess'] <- breedingSuccess
    res.mat[sim, 'maleDepJuvenileSR'] <- survivalRates$maleDepJuvenile
    res.mat[sim, 'maleJuvenileSR'] <- survivalRates$maleJuvenile
    res.mat[sim, 'maleSubadultSR'] <- survivalRates$maleSubadult
    res.mat[sim, 'maleAdultSR'] <- survivalRates$maleAdult
    res.mat[sim, 'femaleDepJuvenileSR'] <- survivalRates$femaleDepJuvenile
    res.mat[sim, 'femaleJuvenileSR'] <- survivalRates$femaleJuvenile
    res.mat[sim, 'femaleSubadultSR'] <- survivalRates$femaleSubadult
    res.mat[sim, 'femaleAdultSR'] <- survivalRates$femaleAdult
    res.mat[sim, 'translocatedSexRatio'] <- translocatedSexRatio
    res.mat[sim, 'tYear1'] <- translocations[[1]][["tYear"]]
    res.mat[sim, 'nYear1'] <- translocations[[1]][["numIndividuals"]]
    res.mat[sim, 'tYear2'] <- translocations[[2]][["tYear"]]
    res.mat[sim, 'nYear2'] <- translocations[[2]][["numIndividuals"]]
    res.mat[sim, 'environmentStochasticity'] <- environmentStochasticity
    res.mat[sim, 'carryingCapacity'] <- carryingCapacity
    res.mat[sim, 'nYears'] <- nYears
    res.mat[sim, 'numSimulations'] <- numSimulations
    res.mat[sim, 'Nfinal'] <- populationHistory[nYears]
    res.mat[sim, 'Extant'] <- ifelse((sum(males) == 0 | sum(females) == 0),0,1)
    res.mat[sim, 'Nmin'] <- min(populationHistory)
    res.mat[sim, 'N5yr'] <- populationHistory[24]
    res.mat[sim, 'lambda'] <- lambda
    res.mat[sim, 'r'] <- r
    
    
  }
  
  # Return the results
  return(res.mat)
  #return(results)
}

## compile the function so it runs faster
RTP.projection.func.comp <- cmpfun(RTP.projection.func)




##########################################################################################
## Sensitivity Analysis Function
## This function requires the package 'lhs', and specification of the following arguments:
## SAvars - vector of parameter names to be included
## nSims - total number of simulations
## nSamples - number of parameter samples to draw
## type - type of sampling to implement, currently must be 'random' or 'latin'
## NOTE: the number of model iterations per sample (nIter) is calculated as nSims/nSample
##References list parRange.list which is a list of variables included in the model, and their ranges for a sensitivity analysis i.e. parRange.list$nStart <- c(1,101)
##Also references list parDefault.list with default values for these variables i.e. parDefault.list$nStart <- 50
##########################################################################################


SA.func <- function(SAvars, nSims, nSamples, type, ...) {
  
  ## load packages
  if(type == 'latin') library(lhs)
  
  ## calculate number of iterations per sample
  numSimulations <- nSims / nSamples
  
  ####################################################
  ## Generate parameter samples
  ####################################################
  
  ## first set up template with default values
  samples <- expand.grid(parDefault.list)  # assuming parDefault.list is defined elsewhere
  samples <- samples[rep(seq_len(nrow(samples)), each = nSamples), ]
  
  nVars <- length(SAvars)
  
  ## generate uniform samples between 0 and 1 for required parameters
  if (type == 'random') {
    raw.samples <- matrix(NA, nrow = nSamples, ncol = nVars)
    for(i in 1:nVars) {
      raw.samples[,i] <- runif(nSamples, 0, 1)
    }
  } else if (type == 'latin') {
    raw.samples <- randomLHS(n = nSamples, k = nVars)
  }
  
  ## transform using required ranges, and convert values that need to be integers to integers.
  for(i in 1:length(SAvars)) {
    if (SAvars[i]%in%c('initialClass','nInitial', 'tYear1', 'nYear1', 'tYear2', 'nYear2' )) {
      samples[,SAvars[i]] <- round(qunif(raw.samples[,i], min=parRange.list[[SAvars[i]]][1], max=parRange.list[[SAvars[i]]][2]))
    } else {
      samples[,SAvars[i]] <- qunif(raw.samples[,i], min=parRange.list[[SAvars[i]]][1], max=parRange.list[[SAvars[i]]][2])
      }
    }

  #make 20% of reintroductions start with zero animals. i.e. get more variation in reintroduction strategies, by making some simulations have only 1 or no follow up translocations
  set.seed(123)
  zero.indices <- sample(seq_len(nrow(samples)), size = 0.2*nrow(samples))
  samples[zero.indices, 'nYear1'] <- 0
  set.seed(234)
  zero.indices <- sample(seq_len(nrow(samples)), size = 0.2*nrow(samples))
  samples[zero.indices, 'nYear2'] <- 0
  
  samples$nSims <- nSims
  samples$nSamples <- nSamples
  samples$numSimulations <- numSimulations
  
  ## run the RTP projection function for each parameter set
  result <- foreach(row = 1:nrow(samples), .combine = rbind, .packages=c('MASS'), .export = c('RTP.projection.func.comp')) %dopar% {

    # Run the RTP simulation
    res.mat <- RTP.projection.func.comp(
      nYears = samples$nYears[row],
      nInitial = samples$nInitial[row],
      initialClass = samples$initialClass[row],
      releasedSexRatio = samples$releasedSexRatio[row],
      translocatedSexRatio = samples$translocatedSexRatio[row],
      survivalRates = list(
        maleDepJuvenile = samples$maleDepJuvenileSR[row],
        maleJuvenile = samples$maleJuvenileSR[row],
        maleSubadult = samples$maleSubadultSR[row],
        maleAdult = samples$maleAdultSR[row],
        femaleDepJuvenile = samples$femaleDepJuvenileSR[row],
        femaleJuvenile = samples$femaleJuvenileSR[row],
        femaleSubadult = samples$femaleSubadultSR[row],
        femaleAdult = samples$femaleAdultSR[row]),
      sexRatio = samples$sexRatio[row],
      fecundityRateMean = samples$fecundityRateMean[row],
      fecundityRateSD = samples$fecundityRateSD[row],
      maxOffspring = samples$maxOffspring[row],
      breedingSuccess = samples$breedingSuccess[row],
      environmentStochasticity = samples$environmentStochasticity[row],
      carryingCapacity = samples$carryingCapacity[row],
      numSimulations = samples$numSimulations[row],
      translocations = list(
        list(tYear = samples$tYear1[row], numIndividuals = samples$nYear1[row]),
        list(tYear = samples$tYear2[row], numIndividuals = samples$nYear2[row])) 
    )
    
    # Return the results
    return(res.mat)
  }
  
  ## return the results
  return(result)
}


################################################################################################################################
## Emulation Function
## This function emulates the sensitivity analysis output with boosted regression trees with different interaction depths
## It takes the following arguments:
## data - the sensitivity analysis output to use (produced by the function SA.func above)
## SAvars - vector of parameter names to be included
## resp - the focal response variable
## subsample - vector of subsamples (i.e., number of data rows) for which emulation will be performed
## tree.complexities - vector of tree complexities (interaction depths) to test
################################################################################################################################

#Note, 19/2/2025 After getting the emulator to run correctly on the HPC, i have found that cross-validation deviance is increasing as sample size increases, not stabilizing
#I need to tune the gbms to this particular dataset, altering bag fraction, lr, n.folds etc so that over and under fitting is constrained.
#original lr = 0.01, n.folds = 5, bag.fraction = 0.5
#20/2/2025 using 10 folds, lr = 0.01, bag.fraction 0.5 for 10,000 samples CVdev = 0.22. Potentially, CVdev goes up before it comes down, i.e. 20,000 samples might stabilise

emulation.func <- function(data, SAvars, resp, subsample, tree.complexities, ...) {
  
  require(dismo)
  
  ## Subset required simulation results
  dataset <- data.frame(data[1:subsample, ])
  dataset$r[dataset$r %in% c('-Inf', 'Inf')] <- NA
  
  ## Statistical distribution for fitting BRTs
  brt.dist <- ifelse(resp == 'Extant', 'bernoulli', 'gaussian')
  
  ## Store BRT models in a list
  brt.models <- list()
  
  ## Fit BRT emulators of different tree complexities for given response variable
  for (i in seq_along(tree.complexities)) {
    
    tc <- tree.complexities[i]
    brt.fit <- NULL
    
    if ((resp == 'Extant' & length(unique(dataset$Extant)) > 1) | resp == 'r') {
      x.col <- which(names(dataset) %in% SAvars)
      y.col <- which(names(dataset) == resp)
      brt.dataset <- dataset[!is.na(dataset[, y.col]), ]
      
      brt.fit <- try(dismo::gbm.step(
        data = brt.dataset, gbm.x = x.col, gbm.y = y.col, learning.rate = 0.01, bag.fraction = 0.5,
        family = brt.dist, tree.complexity = tc, n.folds = 10,
        tolerance.method = 'auto', max.trees = 200000
      ))
      
      ## Try again with a decreased learning rate if necessary
      indic <- 1
      while ((indic <= 49) & (inherits(brt.fit, 'try-error') | is.null(brt.fit))) {
        lr <- max(0.001 - 0.0002 * indic, 0.0001)  # adaptive learning rate that decreases as indic increases (i.e. previous attempt failed, decrease learning rate). Never goes below 0.0001 though, and cuts out after 50 attempts.
        brt.fit <- try(dismo::gbm.step(
          data = brt.dataset, gbm.x = x.col, gbm.y = y.col,
          family = brt.dist, tree.complexity = tc, n.folds = 5,
          tolerance.method = 'auto', tolerance = tol,
          learning.rate = lr, step.size = ss, max.trees = 200000
        ))
        indic <- indic + 1
      }
      
      ## If model fits successfully, store in list
      if (!inherits(brt.fit, 'try-error') & !is.null(brt.fit)) {
        brt.models[[paste0("brt_tc_", tc, "_", resp)]] <- brt.fit
      }
    }
  }
  
  return(brt.models)  # Return the list of fitted BRT models
}


 ###################################################################################################################################
## Emulation Summary Function
## This function summarises the emulation results produced by 'emulation.func' above
## It returns a list of 2 data frames that store:
## the cross-validation deviance for each emulation
## the stability of relative influence metrics results as the subsample size is increased. Beta-diversity (see Toms paper)
###################################################################################################################################

emulation.summary.func <- function(brt.models=brt.models, resp=resp, subsamples=subsamples, tree.complexities=tree.complexities) {
  
  #create ed function from MDM (which has been pulled from CRAN) that calculates beta diversities
  ed <- function (x, q = 1, w = 1, retq = TRUE){
    bsums <- function(x, q = 1, w = 1) {
      if (all(x == 0))
        return(0)
      if (q == 0)
        sum(x > 0)
      else if (q == 1)
        sum(-w * x * log(x), na.rm = TRUE)
      else sum(w * x^q)
    }
    rs <- rowSums(x)
    if (length(w) != 1 & length(w) != nrow(x))
      cat("Warning: n weights NE n rows of x !!")
    if (any(rs == 0)) {
      drops <- (1:nrow(x))[rs == 0]
      rs <- rs[rs != 0]
      x <- x[rs != 0,]
      w <- w[rs != 0]
      cat("Dropping zero sum rows: ", drops, "\n")
    }
    x <- x/rs
    wa <- w^q/mean(w^q)
    a <- bsums(x, q = q, w = wa)/nrow(x)
    wg <- w/mean(w)
    g <- bsums(colMeans(wg * x, na.rm = TRUE), q = q, w = 1)
    if (retq) {
      if (q == 1) {
        a <- exp(a)
        g <- exp(g)
      }
      else if (q!= 1) {
        a <- a^(1/(1 - q))
        g <- g^(1/(1 - q))
      }
      c(alpha = a, beta = g/a, gamma = g)
    }
    else c(absums=a, gbsums=g)
  }
  
  #Create data frames to store cross-validation deviance and the stability of relative influence metrics results as the subsample size is increased
  cvDev.df <- betaDiv.df <- data.frame(matrix(NA, nrow=length(subsamples), ncol=length(tree.complexities)+1))
  names(cvDev.df) <- names(betaDiv.df) <- c('subsample',paste('tc',tree.complexities,sep='.'))
  cvDev.df$subsample <- betaDiv.df$subsample <- subsamples
  ri.df <- NULL
  
# Initialize empty list to store data frames for each tree complexity
ri.list <- list()

for (j in 1:length(tree.complexities)) {  
  tc <- tree.complexities[j]  
  ri.temp <- data.frame(Parameter = sort(SAvars))  
  ri.temp$tc <- tc
  
  for (i in 1:length(brt.models)) {  
    subsample <- subsamples[i]  
    brt <- brt.models[[i]][[j]]  
    ri <- brt$contributions  
    ri <- ri[order(ri$var), ]  
    
    # Store relative influence under the correct subsample column
    col_name <- paste("subsample", subsample, sep = ".")  
    ri.temp[[col_name]] <- ri[, 2]
    
    ## calculate cross-validation deviance
    dev <- brt$cv.statistics$deviance.mean
    cvDev.df[i,j+1] <- dev #add deviance to dataframe in correct row and column (+1 avoids subsample)
  }
  ri.list[[j]] <- ri.temp # Store each tree complexity as a separate element in a list
}
ri.df <- do.call(rbind, ri.list)  # Combine all tree complexities into a single data frame

#calculate beta diversity of pairs (Same way it was done in Tom's paper. Due to pair comparisons, no beta diversity calculated for first subsample)
for (i in 1:length(tree.complexities)) {
  tc <- tree.complexities[i]
  df <- ri.df[ri.df$tc==tc,]
  betaDiv.dummy <- data.frame(subsample=subsamples,resp=NA)
  for (j in 2:length(subsamples)) {
    beta.div <- as.numeric(ed(t(as.matrix(df[,j:(j+1)])),q=1,retq=T)['beta'])
    betaDiv.dummy[j,2] <- beta.div    
  }
  betaDiv.df[,paste('tc',tc,sep='.')] <- betaDiv.dummy[,2]
}

return(list(cvDev = cvDev.df, betaDiv=betaDiv.df))
}


##########################################################################################################
#Function to extract fitted values from gbm objects, for use in ggplot for partial residuals plots
##########################################################################################################

gbm_to_dataframe <- function(gbm.object,                # a gbm object - could be one from gbm.step
                             variable.no = 0,            # the var to extract - if zero then extracts all
                             n.plots = length(pred.names), # extract the first n most important preds
                             show.contrib = TRUE,        # include contribution percentage in the output
                             ...                         # other arguments to pass, e.g., for customization
) {
  if (!requireNamespace('gbm')) { stop('You need to install the gbm package to run this function') }
  requireNamespace('splines')
  
  # Extract gbm object components
  gbm.call <- gbm.object$gbm.call
  gbm.x <- gbm.call$gbm.x
  pred.names <- gbm.call$predictor.names
  response.name <- gbm.call$response.name
  data <- gbm.call$dataframe
  
  # Setup the number of plots to extract
  max.vars <- length(gbm.object$contributions$var)
  if (n.plots > max.vars) {
    n.plots <- max.vars
    warning("Reducing number of extracted predictors to maximum available (", max.vars, ")")
  }
  
  predictors_list <- list()
  responses_list <- list()
  contrib_list <- list()
  predictor_names_list <- list()  # To store the predictor names
  
  # Loop to extract data for the top n predictors
  for (j in 1:n.plots) {
    if (n.plots == 1) {
      k <- variable.no
    } else {
      k <- match(gbm.object$contributions$var[j], pred.names)
    }
    
    var.name <- gbm.call$predictor.names[k]
    pred.data <- data[, gbm.call$gbm.x[k]]
    
    # Get the response matrix from the gbm plot function
    response.matrix <- gbm::plot.gbm(gbm.object, k, return.grid = TRUE)
    
    # Store predictor values, responses, and contributions
    predictors_list[[j]] <- response.matrix[, 1]
    responses_list[[j]] <- response.matrix[, 2] - mean(response.matrix[, 2])
    
    # Optionally include contribution in the dataframe
    if (show.contrib) {
      contrib_list[[j]] <- rep(round(gbm.object$contributions[j, 2], 1), length(predictors_list[[j]]))
    } else {
      contrib_list[[j]] <- rep(NA, length(predictors_list[[j]]))
    }
    
    # Store the predictor names
    predictor_names_list[[j]] <- rep(var.name, length(predictors_list[[j]]))
  }
  
  # Combine all extracted data into a dataframe
  results_df <- data.frame(
    Predictor = unlist(predictor_names_list),
    Predictor_Value = unlist(predictors_list),
    Response = unlist(responses_list),
    Contribution = unlist(contrib_list)
  )
  
  # Return the dataframe
  return(results_df)
}
