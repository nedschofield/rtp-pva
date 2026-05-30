#Individual based RTP PVA v1

  profvis({
################################################################
#functions to compute He from population data frame, with format
#locus1_a1, locus1_a2, locus2_a1, locus2_a2 etc
################################################################

#function to compute He for each locus
He.func <- function(a1, a2) { 
  alleles <- c(a1, a2)  # Combine allele data
  allele_freqs <- table(alleles) / length(alleles)  # Compute allele frequencies
  He <- 1 - sum(allele_freqs^2)  # He formula
  return(He)
}

#compute single population level He
pop.He.func <- function(pop_data) {
  allele_cols <- grep("^locus", names(pop_data), value = TRUE)  # Select loci columns
  locus_names <- unique(gsub("_(a1|a2)$", "", allele_cols)) #get unique locus names
  He_values <- sapply(locus_names, function(locus) {
    He.func(pop_data[[paste0(locus, "_a1")]], pop_data[[paste0(locus, "_a2")]]) # Compute He for each locus
  })
  global_He <- mean(He_values, na.rm = TRUE) # Compute global He as the mean of per-locus He
  return(global_He)
}


#################################################################################################
#Individual based simulation function, with genetics, and multiple paternity
#################################################################################################

#26/2
#need to add in environmental stochasticity
#multiple paternity works, but not accurate to species. Need to vary multiple paternity according to literature.
#need to convert function to use data.table and Rcpp, to increase speed.




##Create the population projection function.

RTP.IBM.func <- function(N, #initial release size
                                initialClass, #age class of initial release i.e. 2:5
                                time_steps = 80, #years of sim. 80 'years' = 20 years as each 'year' is a 3-month period
                                releasedSexRatio = 0.5, #initial sex ratio of 1st release
                                sexRatio = 0.62, fecundityRateMean = 7.5, fecundityRateSD = 0.69, maxOffspring = 8, #reproductive parameters for RTP
                                breedingSuccess = 1, #probability of females breeding. i.e. 1 = all females breed
                                environmentStochasticity = 0.1,
                                K = 500, #carrying capacity
                                translocations = list(
                                  list(tYear = 4, numIndividuals = 50),
                                  list(tYear = 8, numIndividuals = 50)),
                                translocatedSexRatio = 0.5, #sex ratio for all subsequent translocations
                                allele_frequencies,
                                maleSR = c(0.9, 0.85, 0.8, 0),  # Male survival rates for ages 1-4
                                femaleSR = c(0.9, 0.88, 0.83, 0.8),  # Female survival rates for ages 1-3, and 4-13 same rate
                                numSimulations = 1
                                ) {
  
  #load packages
  library(data.table)
  
  # Set up a data frame to store results. initialise all as zero, so no NA's if population goes extinct
  res.mat <- matrix(NA,nr=numSimulations,ncol=33)
  colnames(res.mat) <- c('nInitial', 'initialClass', 'releasedSexRatio', 'sim',
                         'sexRatio', 'fecundityRateMean', 'fecundityRateSD', 'maxOffspring',
                         'breedingSuccess', 
                         'maleDepJuvenileSR', 'maleJuvenileSR', 'maleSubadultSR', 'maleAdultSR', 'femaleDepJuvenileSR', 'femaleJuvenileSR', 'femaleSubadultSR', 'femaleAdultSR', 
                         'translocatedSexRatio', 'tYear1', 'nYear1', 'tYear2', 'nYear2',
                         'environmentStochasticity', 'carryingCapacity', 
                         'nYears',
                         'Nfinal', 'Extant', 'Nmin', 'N5yr', 'lambda', 'r', 'initialHe', 'finalHe')
  
  #results for each year (for function testing)
  results <- data.frame(Generation = integer(time_steps * numSimulations),
                        PopulationSize = integer(time_steps * numSimulations),
                        Simulation = integer(time_steps * numSimulations),
                        season = integer(time_steps * numSimulations))
  
  #run iterations
  for (sim in 1:numSimulations) {
  # Initialize population
  populationHistory <- numeric(time_steps)
  age <- rep(initialClass, N)
  pop <- data.frame(
    id = 1:N,
    age = age,
     #if pregnant female
    sex = ifelse(age == 5, "F", (sample(c("M", "F"), N, replace = TRUE, prob = c(releasedSexRatio, 1-releasedSexRatio)))),
    offspringM = 0,
    offspringF = 0,
    mother_id = NA
  )
  
  # Initialize the season vector (1 = sep-oct-nov, 2 = dec-jan-feb, 3 = mar-apr-may, 4 = june-july-aug)
  season <- rep(1:4, length.out = time_steps)  # Repeat the sequence 1, 2, 3, 4 for the number of years
  season <- ((1:time_steps) + (initialClass - 1)) %% 4    
  season[season == 0] <- 4  # Ensure 0 maps to 4
  
  # Assign genetic information to starting population
  num_loci <- length(allele_frequencies)
  loci_data <- lapply(1:num_loci, function(locus) {
    num_alleles <- length(allele_frequencies[[locus]])
    alleles <- 1:num_alleles
    data.frame(
      a1 = sample(alleles, N, replace = TRUE, prob = allele_frequencies[[locus]]),
      a2 = sample(alleles, N, replace = TRUE, prob = allele_frequencies[[locus]])
    )
  })
  loci_df <- do.call(cbind, loci_data)
  colnames(loci_df) <- unlist(lapply(1:num_loci, function(locus) c(paste0("locus", locus, "_a1"), paste0("locus", locus, "_a2"))))
  pop <- cbind(pop, loci_df)
  
  #create offspring if release was for pregnant females
  if (initialClass == 5){
    num_births <- pmin(maxOffspring, round(rnorm(NROW(pop), mean = fecundityRateMean, sd = fecundityRateSD))) * 
      rbinom(NROW(pop), 1, breedingSuccess)                                                                      #calculate litter sizes, clamp to maxOffspring, then multiply by bernoulli trial 
    pop$offspringM <- sapply(num_births, function(x) sum(rbinom(x, 1, sexRatio)))                                #generate number of M/F in litters using Bernoulli trials
    pop$offspringF <- num_births - pop$offspringM
    
    #create offspring dataframe
    if(sum(num_births) > 0){
      offspringM <- data.frame(
        id = 0, 
        age = 1, #start one, as season initialised as 1. i.e. skip breeding time_step as preg females reintroduced in sep-oct-nov
        sex = "M", 
        offspringM = 0,
        offspringF = 0,
        mother_id = rep(pop$id, pop$offspringM)  # Repeat each mother ID for the number of male offspring
      )
      offspringF <- data.frame(
        id = 0, 
        age = 1, 
        sex = "F",  
        offspringM = 0,
        offspringF = 0,
        mother_id = rep(pop$id, pop$offspringF)  # Repeat each mother ID for the number of female offspring
      )
      offspring <- rbind(offspringM, offspringF)      # Combine male and female offspring into one dataframe
      offspring <- offspring[order(offspring$mother_id), ] #so litters sit together
      offspring$id <- (max(pop$id, na.rm = TRUE) + 1):(max(pop$id, na.rm = TRUE) + sum(num_births)) #create unique ID
    }
    
    ##### Assign genetic information to offspring
    if (NROW(offspring) > 0) {
      mother_ids <- offspring$mother_id
      
      #get an allele from mother, and a random allele from allele_frequences
      loci_inheritance <- lapply(1:num_loci, function(locus) {
        mother_alleles <- mapply(function(id) {pop[id == pop$id, paste0("locus", locus, "_a", sample(1:2, 1))]}, mother_ids)
        num_alleles <- length(allele_frequencies[[locus]])
        alleles <- 1:num_alleles
        father_alleles <- sample(alleles, length(mother_alleles), replace = TRUE, prob = allele_frequencies[[locus]])
        data.frame(
          a1 = mother_alleles,
          a2 = father_alleles
        )
      })
      
      loci_offspring <- do.call(cbind, loci_inheritance)
      colnames(loci_offspring) <- unlist(lapply(1:num_loci, function(locus) c(paste0("locus", locus, "_a1"), paste0("locus", locus, "_a2"))))
      offspring <- cbind(offspring, loci_offspring)
      
      pop <- rbind(pop, offspring)
    }
  }
  
  #Calculate year 1 He, and add to res.mat
  res.mat[sim, 'initialHe'] <- pop.He.func(pop)
  
#### Run the simulation ####  
  for (t in 1:time_steps) {
    
  #### Reproduction process ####
    if (season[t] == 4) {
    breeding_females <- pop[pop$sex == "F" & pop$age %in% c(4,8,12), ]                                                        #get all females that can breed, 
    num_births <- pmin(maxOffspring, round(rnorm(NROW(breeding_females), mean = fecundityRateMean, sd = fecundityRateSD))) * 
      rbinom(NROW(breeding_females), 1, breedingSuccess)                                                                      #calculate litter sizes, clamp to maxOffspring, then multiply by bernoulli trial 
    breeding_females$offspringM <- sapply(num_births, function(x) sum(rbinom(x, 1, sexRatio)))                                #generate number of M/F in litters using Bernoulli trials
    breeding_females$offspringF <- num_births - breeding_females$offspringM
    
    #match_indices <- match(breeding_females$id, pop$id) # Find matching indices in pop for breeding females
    #pop$offspringM[match_indices] <- breeding_females$offspringM # Update offspring counts in pop
    #pop$offspringF[match_indices] <- breeding_females$offspringF
    
    
    #create offspring dataframe
    if(sum(num_births) > 0){
    offspringM <- data.frame(
      id = 0, 
      age = 0,  
      sex = "M", 
      offspringM = 0,
      offspringF = 0,
      mother_id = rep(breeding_females$id, breeding_females$offspringM)  # Repeat each mother ID for the number of male offspring
    )
    offspringF <- data.frame(
      id = 0, 
      age = 0,  
      sex = "F",  
      offspringM = 0,
      offspringF = 0,
      mother_id = rep(breeding_females$id, breeding_females$offspringF)  # Repeat each mother ID for the number of female offspring
    )
    offspring <- rbind(offspringM, offspringF)      # Combine male and female offspring into one dataframe
    offspring <- offspring[order(offspring$mother_id), ] #so litters sit together
    offspring$id <- (max(pop$id, na.rm = TRUE) + 1):(max(pop$id, na.rm = TRUE) + sum(num_births)) #unique ID
    }
    
  ##### Assign genetic information from parents ####
    if (NROW(offspring) > 0) {
      mother_ids <- offspring$mother_id
      fathers <- pop[pop$sex == "M", ]
      
      # Create a father assignment for each offspring
      father_ids <- unlist(lapply(breeding_females$id, function(mother_id) {
        num_offspring <- breeding_females[breeding_females$id == mother_id, "offspringM"] + breeding_females[breeding_females$id == mother_id, "offspringF"]         # Get the number of offspring for each mother
        primary_father <- sample(fathers$id, 1)  # Randomly pick a primary father
        
        num_secondary <- round(0.2 * num_offspring) #assign most offspring to the primary father, and a small fraction to a secondary father
        father_ids <- c(rep(primary_father, num_offspring - num_secondary),
                              sample(fathers$id, num_secondary, replace = TRUE))  # Select secondary fathers randomly for a fraction of the offspring
        sample(father_ids)         # Shuffle the order of fathers (so the distribution is random across offspring)
      }))
      
      
      loci_inheritance <- lapply(1:num_loci, function(locus) {
        mother_alleles <- mapply(function(id) {pop[id == pop$id, paste0("locus", locus, "_a", sample(1:2, 1))]}, mother_ids)
        father_alleles <- mapply(function(id) {pop[id == pop$id, paste0("locus", locus, "_a", sample(1:2, 1))]}, father_ids)
        data.frame(
          a1 = mother_alleles,
          a2 = father_alleles
        )
      })
      
      loci_offspring <- do.call(cbind, loci_inheritance)
      colnames(loci_offspring) <- unlist(lapply(1:num_loci, function(locus) c(paste0("locus", locus, "_a1"), paste0("locus", locus, "_a2"))))
      offspring <- cbind(offspring, loci_offspring)
      
      pop <- rbind(pop, offspring)
      }
    }
    
    
    
  #### Survival process ####
    #first create masks for different age classes
    pouch_young_mask <- pop$age == 0    #pouch young mask -> in period 4 = june-july-aug, no mortality
    female_dep_juv_mask <- pop$sex == "F" & pop$age == 0
    male_mask <- pop$sex == "M" & pop$age > 0 & pop$age <= 4          #all males mask
    female_young_mask <- pop$sex == "F" & pop$age > 0 & pop$age <= 3      #not adult females mask
    female_adult_mask <- pop$sex == "F" & pop$age >= 4 & pop$age <= 13     #adult females mask
    
    #create a vector of survival probabilities 
    survival_probs <- numeric(NROW(pop))
    survival_probs[pouch_young_mask] <- 1 #all pouch young survive breeding period, no need for extra check here
    survival_probs[male_mask] <- maleSR[pop$age[male_mask]] 
    survival_probs[female_young_mask] <- femaleSR[pop$age[female_young_mask]]
    survival_probs[female_adult_mask] <- femaleSR[4]
    
    #adjust all survival probabilities by environmental stochasticity

    #bernoulli trial for survival
    survived <- runif(NROW(pop)) < survival_probs
    
    #remove dependent offspring as well, even if they survived
    dead_mothers_ids <- pop$id[!survived & pop$id %in% pop$mother_id] #get females who died who have dependent juveniles
    juvenile_mask <- !(pop$id %in% pop[pop$age %in% c(0, 1), ]$id & pop$mother_id %in% dead_mothers_ids) # Remove dependent juveniles whose mother_id matches dead mother's id

    # Apply survival to the rest of the population
    pop <- pop[juvenile_mask & survived, ]

    
  #### Aging process ####
    pop$age <- pop$age + 1
    
  
  #### Translocation process ####
    if (length(translocations) > 0) {        # Look for any translocations at this year, if year in any element of the list == t, then that translocation list is selected
      translocationEvent <- translocations[which(sapply(translocations, function(x) x$tYear == t))]
      if (length(translocationEvent) > 0) {        # If there is a translocation event for this year, add individuals
        translocationEvent <- translocationEvent[[1]]  # Get the first translocation event for this year
        numToAdd <- translocationEvent$numIndividuals
        if (numToAdd > 0) {
          age <- rep(ifelse(season[t] == 1, 2, ifelse(season[t] == 2, 3, ifelse(season[t] == 3, 4, 5))), numToAdd) #choose age according to season
          trans <- data.frame(
            id = ((max(pop$id, na.rm = TRUE) + 1):(max(pop$id, na.rm = TRUE) + numToAdd)), #unique ID,
            age = age,
            #if pregnant female
            sex = ifelse(age == 5, "F", (sample(c("M", "F"), N, replace = TRUE, prob = c(translocatedSexRatio, 1-translocatedSexRatio)))),
            offspringM = 0,
            offspringF = 0,
            mother_id = NA
          )
          
          # Assign genetic information from allele frequency list
          num_loci <- length(allele_frequencies)
          loci_data <- lapply(1:num_loci, function(locus) {
            num_alleles <- length(allele_frequencies[[locus]])
            alleles <- 1:num_alleles
            data.frame(
              a1 = sample(alleles, numToAdd, replace = TRUE, prob = allele_frequencies[[locus]]),
              a2 = sample(alleles, numToAdd, replace = TRUE, prob = allele_frequencies[[locus]])
            )
          })
          loci_df <- do.call(cbind, loci_data)
          colnames(loci_df) <- unlist(lapply(1:num_loci, function(locus) c(paste0("locus", locus, "_a1"), paste0("locus", locus, "_a2"))))
          trans <- cbind(trans, loci_df)
          
            if (season[t] == 4) { #this checks if season == 4, i.e. age = 5. preg female translocation
              num_births <- pmin(maxOffspring, round(rnorm(NROW(trans), mean = fecundityRateMean, sd = fecundityRateSD))) * 
                rbinom(NROW(trans), 1, breedingSuccess)                                                                      #calculate litter sizes, clamp to maxOffspring, then multiply by bernoulli trial 
              trans$offspringM <- sapply(num_births, function(x) sum(rbinom(x, 1, sexRatio)))                                #generate number of M/F in litters using Bernoulli trials
              trans$offspringF <- num_births - trans$offspringM
              
              #create offspring dataframe
              if(sum(num_births) > 0){
                offspringM <- data.frame(
                  id = 0, 
                  age = 1, #as ageing step has already taken place 
                  sex = "M", 
                  offspringM = 0,
                  offspringF = 0,
                  mother_id = rep(trans$id, trans$offspringM)  # Repeat each mother ID for the number of male offspring
                )
                offspringF <- data.frame(
                  id = 0, 
                  age = 1, #as ageing step has already taken place  
                  sex = "F",  
                  offspringM = 0,
                  offspringF = 0,
                  mother_id = rep(trans$id, trans$offspringF)  # Repeat each mother ID for the number of female offspring
                )
                offspring <- rbind(offspringM, offspringF)      # Combine male and female offspring into one dataframe
                offspring <- offspring[order(offspring$mother_id), ] #so litters sit together
                offspring$id <- (max(trans$id, na.rm = TRUE) + 1):(max(trans$id, na.rm = TRUE) + sum(num_births)) #create unique ID
              }
              
          
              ##### Assign genetic information to offspring, from mother and then random from allele frequencies?
              if (NROW(offspring) > 0) {
                mother_ids <- offspring$mother_id
                
                #get an allele from each parent
                loci_inheritance <- lapply(1:num_loci, function(locus) {
                  mother_alleles <- mapply(function(id) {trans[id == trans$id, paste0("locus", locus, "_a", sample(1:2, 1))]}, mother_ids)
                  num_alleles <- length(allele_frequencies[[locus]])
                  alleles <- 1:num_alleles
                  father_alleles <- sample(alleles, length(mother_alleles), replace = TRUE, prob = allele_frequencies[[locus]])
                  data.frame(
                    a1 = mother_alleles,
                    a2 = father_alleles
                  )
                })
                
                loci_offspring <- do.call(cbind, loci_inheritance)
                colnames(loci_offspring) <- unlist(lapply(1:num_loci, function(locus) c(paste0("locus", locus, "_a1"), paste0("locus", locus, "_a2"))))
                offspring <- cbind(offspring, loci_offspring)
                
                trans <- rbind(trans, offspring)
              }
            }
          
          
          pop <- rbind(pop, trans) #append to population
        }
      }
    }
          

    
  #### Enforce carrying capacity ####
  #If the population exceeds capacity, reduce it via probabilistic truncation, as in Vortex
  #truncation in this way is probably realistic for phascogales, which are territorial. This basically simulates density dependent dispersal outside of the study area.
    populationSize <- NROW(pop)
    if (populationSize > K) {
      excess <- populationSize - K
      adjustmentFactor <- K / populationSize  # 1 / (1 + (excess / totalPopulation))
      
      truncation <- runif(NROW(pop)) < adjustmentFactor
      pop <- pop[truncation, ]
    }
    
    # Calculate the total population size
    populationHistory[t] <- populationSize
    
    # Store results of population for each year, for function testing
    results$Generation[(sim - 1) * time_steps + t] <- t
    results$PopulationSize[(sim - 1) * time_steps + t] <- populationSize
    results$Simulation[(sim - 1) * time_steps + t] <- sim
    results$season[(sim - 1) * time_steps + t] <- season[t]
    
    # Stop the simulation if the population goes extinct
    males <- pop[pop$sex == "M",]
    females <- pop[pop$sex == "F",]
    if (NROW(males) == 0 | NROW(females) == 0) {
      break
      } 
    
  }
  
  ## calculate mean lambda and r
  (lambdas <- populationHistory[2:(time_steps+1)]/populationHistory[1:time_steps])
  (lambda <- mean(lambdas,na.rm=T))
  (r <- log(lambda))
  if(r%in%c('-Inf','Inf')) {r <- NA}
  
  ## store summary outputs of each simulation run
  res.mat[sim, 'sim'] <- sim
  res.mat[sim, 'nInitial'] <- N
  res.mat[sim, 'initialClass'] <- initialClass
  res.mat[sim, 'releasedSexRatio'] <- releasedSexRatio
  res.mat[sim, 'sexRatio'] <- sexRatio
  res.mat[sim, 'fecundityRateMean'] <- fecundityRateMean
  res.mat[sim, 'fecundityRateSD'] <- fecundityRateSD
  res.mat[sim, 'maxOffspring'] <- maxOffspring
  res.mat[sim, 'breedingSuccess'] <- breedingSuccess
  
  res.mat[sim, 'maleDepJuvenileSR'] <- maleSR[1]
  res.mat[sim, 'maleJuvenileSR'] <- maleSR[2]
  res.mat[sim, 'maleSubadultSR'] <- maleSR[3]
  res.mat[sim, 'maleAdultSR'] <- maleSR[4]
  res.mat[sim, 'femaleDepJuvenileSR'] <- femaleSR[1]
  res.mat[sim, 'femaleJuvenileSR'] <- femaleSR[2]
  res.mat[sim, 'femaleSubadultSR'] <- femaleSR[3]
  res.mat[sim, 'femaleAdultSR'] <- femaleSR[4]
  
  res.mat[sim, 'translocatedSexRatio'] <- translocatedSexRatio
  res.mat[sim, 'tYear1'] <- translocations[[1]][["tYear"]]
  res.mat[sim, 'nYear1'] <- translocations[[1]][["numIndividuals"]]
  res.mat[sim, 'tYear2'] <- translocations[[2]][["tYear"]]
  res.mat[sim, 'nYear2'] <- translocations[[2]][["numIndividuals"]]
  res.mat[sim, 'environmentStochasticity'] <- environmentStochasticity
  res.mat[sim, 'carryingCapacity'] <- K
  res.mat[sim, 'nYears'] <- time_steps
  res.mat[sim, 'Nfinal'] <- populationHistory[time_steps]
  res.mat[sim, 'Extant'] <- ifelse((NROW(males) == 0 | NROW(females) == 0), 0, 1)
  res.mat[sim, 'Nmin'] <- min(populationHistory)
  res.mat[sim, 'N5yr'] <- populationHistory[24]
  res.mat[sim, 'lambda'] <- lambda
  res.mat[sim, 'r'] <- r
  res.mat[sim, 'finalHe'] <- pop.He.func(pop) 

  }
  # Return the results
  return(res.mat)
  #return(results)
  
}





# Example run
maleSR <- c(0.5, 0.5, 0.5, 0)  # Male survival rates for ages 1-4
femaleSR <- c(0.5, 0.58, 0.53, 0.5)  # Female survival rates for ages 1-3, and 4-13 same rate
K <- 500  # Carrying capacity

#test values to run code inside function
N <- 100  # Change as needed
initialClass <- 4
time_steps <- 80  
releasedSexRatio <- 0.5  
sexRatio <- 0.62  
fecundityRateMean <- 7.5  
fecundityRateSD <- 0.69  
maxOffspring <- 8  
breedingSuccess <- 1  
environmentStochasticity <- 0.1  
K <- 500  
translocations <- list(
  list(tYear = 4, numIndividuals = 50),
  list(tYear = 8, numIndividuals = 50)
)
translocatedSexRatio <- 0.5  
maleSR <- c(0.9, 0.85, 0.8, 0)  # Male survival rates for ages 1-4
femaleSR <- c(0.9, 0.88, 0.83, 0.8)  # Female survival rates for ages 1-3, and same rate for 4-13

# Example allele frequencies
allele_frequencies <- list(
  c(0.33, 0.33, 0.34),  # Locus 1
  c(0.25, 0.25, 0.25, 0.25),  # Locus 2
  c(0.5, 0.5)  # Locus 3
)


  final_population <- RTP.IBM.func(initialClass = 2, N = 100, time_steps = 80, 
                                        maleSR = maleSR, 
                                        femaleSR = femaleSR,
                                        K = K, 
                                        allele_frequencies = allele_frequencies)



