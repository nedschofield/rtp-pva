#Individual based RTP PVA v3

#10/3
# tessa gave me instructions for improving the multiple paternity part to make it more realistic
#don't worry about EV at this stage, as variation is already incorporated in sensitivity analysis design.


#load packages
library(data.table)
library(compiler)
library(iterators)
library(snow)
library(doSNOW)
library(foreach)
library(MASS)
library(dismo)
library(gbm) #not sure this is actually required, gbm.step wouldn't run on pawsey without gbm installed, even though gbm.step is from dismo??



  #test values to run with prof vis without running as function
  numSimulations=1
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
  
  # Define the function to generate the list with 1000 vectors of 2 random alleles
  # Set the number of loci and number of vectors
  num_loci <- 1000
  
  # Generate 1000 random vectors with 2 alleles each (Dirichlet distribution)
  allele_frequencies <- lapply(1:num_loci, function(x) rdirichlet(1, c(1, 1)))
  
  # Display the first few frequency vectors
  head(allele_frequencies)
  
  
  
  
  
#########################################################
#import actual allele frequencies, and convert to list
#########################################################
  
ONE_pop_n.loci_Vortex1000L_new <- read.table("~/Research/Bush Heritage/PVA_secret-rocks-RTP/data/ONE_pop_n.loci_Vortex1000L_new.txt", quote="\"", comment.char="")
allele_frequencies <- lapply(split(ONE_pop_n.loci_Vortex1000L_new, seq(nrow(ONE_pop_n.loci_Vortex1000L_new))), as.double)
  
  
################################################################
#functions to compute He from population data frame, with format
#locus1_a1, locus1_a2, locus2_a1, locus2_a2 etc
################################################################

#function to compute He for each locus
  He.func <- function(a1, a2) { 
    alleles <- c(a1, a2)  
    allele_counts <- tabulate(alleles, nbins = max(alleles, na.rm = TRUE))  # Faster than table()
    allele_freqs <- allele_counts / sum(allele_counts)  # Normalize to get frequencies
    He <- 1 - sum(allele_freqs^2, na.rm = TRUE)  # He formula
    return(He)
  }

#compute single population level He
  pop.He.func <- function(pop_data) {
    allele_cols <- grep("^locus", names(pop_data), value = TRUE)  # Find all allele columns
    locus_names <- unique(gsub("_(a1|a2)$", "", allele_cols))  # Extract locus names
    He_values <- pop_data[, lapply(locus_names, function(locus) {
      He.func(get(paste0(locus, "_a1")), get(paste0(locus, "_a2")))  # Compute He for each locus
    })]
    global_He <- mean(unlist(He_values), na.rm = TRUE)  # Compute mean He across loci
    return(global_He)
  }

  
#################################################################################################
#Individual based simulation function, with genetics, and multiple paternity
#################################################################################################

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
                         maleSR = c(0.9, 0.85, 0.8, 0),  # Male survival rates for ages 1-4
                         femaleSR = c(0.9, 0.88, 0.83, 0.8),  # Female survival rates for ages 1-3, and 4-13 same rate
                         numSimulations = 1,
                       allele_frequencies
) { 
  
  
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
  #results <- data.frame(Generation = integer(time_steps * numSimulations),
  #                      PopulationSize = integer(time_steps * numSimulations),
  #                      Simulation = integer(time_steps * numSimulations),
  #                      season = integer(time_steps * numSimulations),
  #                      tHe = integer(time_steps*numSimulations))
  
  #run iterations
  for (sim in 1:numSimulations) {
    
    #initialise storage vectors
    populationHistory <- numeric(time_steps)
    males <- numeric(time_steps)
    females <- numeric(time_steps)
    
    num_loci <- length(allele_frequencies)
    #memory preallocation for reproduction
    #max_possible_offspring <- 5000
    #mother_allele_choices_R <- matrix(0L, nrow = max_possible_offspring, ncol = num_loci)
    #father_allele_choices_R <- matrix(0L, nrow = max_possible_offspring, ncol = num_loci)
    
    # Initialize population
    age <- rep(initialClass, N)
    pop <- data.table(
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
    allele_matrix <- matrix(nrow = N, ncol = num_loci * 2)  #create empty matrix, 2 columns per locus (a1, a2)
    for (locus in 1:num_loci) {                             #for loop is faster than lapply in this case, as it operates column-wise. lapply would have function overhead for column and row.
      num_alleles <- length(allele_frequencies[[locus]]) 
      alleles <- 1:num_alleles  # Define possible allele values
      allele_matrix[, (2*locus - 1)] <- sample(alleles, N, replace = TRUE, prob = allele_frequencies[[locus]])
      allele_matrix[, (2*locus)] <- sample(alleles, N, replace = TRUE, prob = allele_frequencies[[locus]])
    }
    loci_df <- as.data.table(allele_matrix)
    setnames(loci_df, 
             old = names(loci_df), 
             new = unlist(lapply(1:num_loci, function(locus) c(paste0("locus", locus, "_a1"), paste0("locus", locus, "_a2")))))
    
    # Merge with population
    pop[, (names(loci_df)) := loci_df]
    
    #create offspring if release was for pregnant females
    if (initialClass == 5){
      num_births <- pmin(maxOffspring, round(rnorm(NROW(pop), mean = fecundityRateMean, sd = fecundityRateSD))) * 
        rbinom(NROW(pop), 1, breedingSuccess)                                                                      #calculate litter sizes, clamp to maxOffspring, then multiply by bernoulli trial 
      pop[, offspringM := rbinom(.N, num_births, sexRatio)]  #generate number of M/F in litters using Bernoulli trials
      pop[, offspringF := num_births - offspringM]  
      
      #create offspring dataframe
      if(sum(num_births) > 0){
        offspringM <- data.table(
          id = 0, 
          age = 1, #start one, as season initialised as 1. i.e. skip breeding time_step as preg females reintroduced in sep-oct-nov
          sex = "M", 
          offspringM = 0,
          offspringF = 0,
          mother_id = rep(pop[,id], pop[,offspringM])  # Repeat each mother ID for the number of male offspring
        )
        offspringF <- data.table(
          id = 0, 
          age = 1, 
          sex = "F",  
          offspringM = 0,
          offspringF = 0,
          mother_id = rep(pop[,id], pop[,offspringF])  # Repeat each mother ID for the number of female offspring
        )
        offspring <- rbind(offspringM, offspringF)      # Combine male and female offspring
        setorder(offspring, mother_id) #sort so litters sit together
        offspring[, id := pop[,max(id)+1]:pop[,max(id) + sum(num_births)]] #create unique ID
        
        
        ##### Assign genetic information to offspring
        if (NROW(offspring) > 0) {
          #create breeding colony father genetics (1 per litter)
          mother_ids <- offspring[,mother_id]           #get mother id's
          unique_mothers <- unique(mother_ids)          # Get unique mothers
          
          # Create an empty matrix for father alleles (one per unique mother) and fill with alleles from allele_frequencies
          allele_matrix_fathers <- matrix(nrow = length(unique_mothers), ncol = num_loci * 2)           
          for (locus in 1:num_loci) {  
            num_alleles <- length(allele_frequencies[[locus]])  
            alleles <- 1:num_alleles  # Define possible allele values
            allele_matrix_fathers[, (2*locus - 1)] <- sample(alleles, length(unique_mothers), replace = TRUE, prob = allele_frequencies[[locus]])
            allele_matrix_fathers[, (2*locus)] <- sample(alleles, length(unique_mothers), replace = TRUE, prob = allele_frequencies[[locus]])
          }
          
          father_genetics <- as.data.table(allele_matrix_fathers)
          setnames(father_genetics, 
                   old = names(father_genetics), 
                   new = unlist(lapply(1:num_loci, function(locus) c(paste0("locus", locus, "_a1"), paste0("locus", locus, "_a2")))))
          
          father_loci_df <- father_genetics[match(mother_ids, unique_mothers), ]           # Expand father genetics to match all offspring. Each offspring gets the same father genetics as their littermates
          
          
          # Generate random allele selections for mothers and fathers
          mother_allele_choices <- matrix(sample(1:2, length(mother_ids) * num_loci, replace = TRUE), nrow = length(mother_ids), ncol = num_loci)
          father_allele_choices <- matrix(sample(1:2, length(father_loci_df) * num_loci, replace = TRUE), nrow = length(father_loci_df), ncol = num_loci)
          
          # Get mother and father rows using fast join
          mothers <- pop[J(mother_ids), on = "id"]
          fathers <- father_loci_df
          
          # Define column names for loci
          locus_columns_a1 <- paste0("locus", 1:num_loci, "_a1")
          locus_columns_a2 <- paste0("locus", 1:num_loci, "_a2")
          
          # Extract all allele values into separate data.tables
          mother_alleles_a1 <- mothers[, ..locus_columns_a1]   # Mothers' first allele
          mother_alleles_a2 <- mothers[, ..locus_columns_a2]   # Mothers' second allele
          father_alleles_a1 <- fathers[, ..locus_columns_a1]   # Fathers' first allele
          father_alleles_a2 <- fathers[, ..locus_columns_a2]   # Fathers' second allele
          
          # Select alleles dynamically based on earlier sample (1 = _a1, 2 = _a2)
          loci_inheritance_mother <- mother_alleles_a1 * (mother_allele_choices == 1) + 
            mother_alleles_a2 * (mother_allele_choices == 2)
          loci_inheritance_father <- father_alleles_a1 * (father_allele_choices == 1) + 
            father_alleles_a2 * (father_allele_choices == 2)
          
          # Convert results to data.table
          setnames(loci_inheritance_father, 
                   old = names(loci_inheritance_father), 
                   new = gsub("_a1$", "_a2", names(loci_inheritance_father)))
          loci_inheritance <- cbind(loci_inheritance_mother, loci_inheritance_father)
          
          # sort columns
          col_order <- names(loci_inheritance)             # Get the column names
          locus_numbers <- as.numeric(sub("locus(\\d+)_a[12]", "\\1", col_order))  # Extracts numbers
          allele_suffix <- sub(".*_(a[12])", "\\1", col_order)  # Extracts a1 or a2
          sorted_order <- order(locus_numbers, allele_suffix)             # Order by locus number, then by suffix (a1 before a2)
          setcolorder(loci_inheritance, col_order[sorted_order]) # Reorder columns
          
          offspring[,(names(loci_inheritance)) := loci_inheritance]
          pop <- rbind(pop, offspring)
        }
      }
    }
    
    #Calculate year 1 He, and add to res.mat
    res.mat[sim, 'initialHe'] <- pop.He.func(pop)
    
    
    
    ######################################################################
    # Run the simulation 
    ######################################################################
    for (t in 1:time_steps) {
      
      ##########################################
      # Reproduction process
      ##########################################
      
      # check if reproduction possible
      if (season[t] == 4){
        males_age4 <- pop[sex == "M" & age == 4]   # Extract all relevant males
        breeding_females <- pop[sex == "F" & age %in% c(4,8,12), ] #get all females that can breed, 
        
        if (NROW(males_age4) > 0 & NROW(breeding_females) > 0) {
        num_births <- pmin(maxOffspring, round(rnorm(NROW(breeding_females), mean = fecundityRateMean, sd = fecundityRateSD))) * 
          rbinom(NROW(breeding_females), 1, breedingSuccess)                                                                      #calculate litter sizes, clamp to maxOffspring, then multiply by bernoulli trial 
        breeding_females[, offspringM := rbinom(.N, num_births, sexRatio)]  #generate number of M/F in litters using Bernoulli trials
        breeding_females[, offspringF := num_births - offspringM]  
        
        #create offspring dataframe
        offspringM <- data.table(
          id = 0, 
          age = 0, 
          sex = "M", 
          offspringM = 0,
          offspringF = 0,
          mother_id = rep(breeding_females[,id], breeding_females[,offspringM])  # Repeat each mother ID for the number of male offspring
        )
        offspringF <- data.table(
          id = 0, 
          age = 0, 
          sex = "F",  
          offspringM = 0,
          offspringF = 0,
          mother_id = rep(breeding_females[,id], breeding_females[,offspringF])  # Repeat each mother ID for the number of female offspring
        )
        offspring <- rbindlist(list(offspringM, offspringF))      # Combine male and female offspring
        setorder(offspring, mother_id) #sort so litters sit together
        offspring[, id := pop[,max(id)+1]:pop[,max(id) + sum(num_births)]] #create unique ID
        
        
        ##### Assign genetic information from parents ####
        if (NROW(offspring) > 0) {
          #get mother_ids
          mother_ids <- offspring[,mother_id]
          
          
          #get father_ids
          fathers <- pop[sex == "M" & age == 4, ]
          primary_fathers <- sample(fathers[,id], nrow(breeding_females), replace = TRUE)           # Assign a primary father to each litter/female
          secondary_fathers <- sample(fathers[, id], nrow(breeding_females), replace = TRUE)        # Assign a potential secondary father to each litter/female
          #num_secondary <- round(0.2 * num_births)           
          num_secondary <- rpois(num_births, lambda = 1) #compute portion of litter that is fathered by secondary father           #multiple paternity. between 45 and 80% of litters have a second father.second father number of offspring sampled from poisson dist. usually 1, up to 7 offspring.

          
          # Create father_id list
          father_ids_list <- mapply(function(n, primary, secondary, sec_count) {
            assigned_fathers <- c(rep(primary, n - min(n, sec_count)), rep(secondary, min(n, sec_count)))
            sample(assigned_fathers)  # Shuffle within litter
          }, num_births, primary_fathers, secondary_fathers, num_secondary, SIMPLIFY = FALSE)
          father_ids <- unlist(father_ids_list)           # Unlist to get final vector of father IDs
          
          # Generate random allele selections for mothers and fathers
          mother_allele_choices <- matrix(fifelse(runif(length(mother_ids) * num_loci) > 0.5, 2, 1), nrow = length(mother_ids), ncol = num_loci)
          father_allele_choices <- matrix(fifelse(runif(length(father_ids) * num_loci) > 0.5, 2, 1), nrow = length(father_ids), ncol = num_loci)
          
          # Get mother and father rows using fast join
          mothers <- pop[J(mother_ids), on = "id"]
          fathers <- pop[J(father_ids), on = "id"]
          
          
          # Define column names for loci
          locus_columns_a1 <- paste0("locus", 1:num_loci, "_a1")
          locus_columns_a2 <- paste0("locus", 1:num_loci, "_a2")
          
          # Extract all allele values into separate data.tables
          mother_alleles_a1 <- mothers[, ..locus_columns_a1]   # Mothers' first allele
          mother_alleles_a2 <- mothers[, ..locus_columns_a2]   # Mothers' second allele
          father_alleles_a1 <- fathers[, ..locus_columns_a1]   # Fathers' first allele
          father_alleles_a2 <- fathers[, ..locus_columns_a2]   # Fathers' second allele
          
          # Select alleles dynamically based on earlier sample (1 = _a1, 2 = _a2)
          loci_inheritance_mother <- mother_alleles_a1 * (mother_allele_choices == 1) + 
            mother_alleles_a2 * (mother_allele_choices == 2)
          loci_inheritance_father <- father_alleles_a1 * (father_allele_choices == 1) + 
            father_alleles_a2 * (father_allele_choices == 2)
          
          # Convert results to data.table
          setnames(loci_inheritance_father, 
                   old = names(loci_inheritance_father), 
                   new = gsub("_a1$", "_a2", names(loci_inheritance_father)))
          loci_inheritance <- loci_inheritance_mother[, names(loci_inheritance_father) := loci_inheritance_father]

          
          # sort columns
          col_order <- names(loci_inheritance)             # Get the column names
          locus_numbers <- as.numeric(sub("locus(\\d+)_a[12]", "\\1", col_order))  # Extracts numbers
          allele_suffix <- sub(".*_(a[12])", "\\1", col_order)  # Extracts a1 or a2
          sorted_order <- order(locus_numbers, allele_suffix)             # Order by locus number, then by suffix (a1 before a2)
          setcolorder(loci_inheritance, col_order[sorted_order]) # Reorder columns
          
          offspring[,(names(loci_inheritance)) := loci_inheritance]
          pop <- rbindlist(list(pop, offspring))
        }
      }
    }
      
      
      #############################################
      # Survival process 
      #############################################
      
      #first create masks for different age classes
      pouch_young_mask <- pop$age == 0    #pouch young mask -> in period 4 = june-july-aug, no mortality
      female_dep_juv_mask <- pop$sex == "F" & pop$age == 0
      male_mask <- pop$sex == "M" & pop$age > 0 & pop$age <= 4          #all males mask
      female_young_mask <- pop$sex == "F" & pop$age > 0 & pop$age <= 3      #not adult females mask
      female_adult_mask <- pop$sex == "F" & pop$age >= 4 & pop$age <= 13     #adult females mask
      
      #create a vector of survival probabilities - all adult males get removed
      survival_probs <- numeric(NROW(pop))
      survival_probs[pouch_young_mask] <- 1 #all pouch young survive breeding period, no need for extra check here
      survival_probs[male_mask] <- maleSR[pop$age[male_mask]] 
      survival_probs[female_young_mask] <- femaleSR[pop$age[female_young_mask]]
      survival_probs[female_adult_mask] <- femaleSR[4]
      
      #bernoulli trial for survival
      survived <- runif(NROW(pop)) < survival_probs
      
      #remove dependent offspring as well, even if they survived
      dead_mothers_ids <- pop$id[!survived & pop$sex == "F" & pop$id %in% pop$mother_id]
      juvenile_mask <- !pop$id %in% pop[pop$age %in% c(0, 1) & pop$mother_id %in% dead_mothers_ids, id]
      
      # Apply survival to the rest of the population
      pop <- pop[survived & juvenile_mask, ]
      
      
      
      #################################
      # Aging process 
      ################################# 
      pop[, age := age +1]

      
      
      
      #########################################
      # Translocation process 
      ######################################### 
      if (length(translocations) > 0) {        # Look for any translocations at this year, if year in any element of the list == t, then that translocation list is selected
        translocationEvent <- translocations[which(sapply(translocations, function(x) x$tYear == t))]
        if (length(translocationEvent) > 0) {        # If there is a translocation event for this year, add individuals
          translocationEvent <- translocationEvent[[1]]  # Get the first translocation event for this year
          numToAdd <- translocationEvent$numIndividuals
          if (numToAdd > 0) {
            age <- rep(ifelse(season[t] == 1, 2, ifelse(season[t] == 2, 3, ifelse(season[t] == 3, 4, 5))), numToAdd) #choose age according to season
            trans <- data.table(
              id = pop[,max(id)+1]:pop[,max(id) + numToAdd], #unique ID,
              age = age,
              #if pregnant female
              sex = ifelse(age == 5, "F", (sample(c("M", "F"), numToAdd, replace = TRUE, prob = c(translocatedSexRatio, 1-translocatedSexRatio)))),
              offspringM = 0,
              offspringF = 0,
              mother_id = NA
            )
            
            # Assign genetic information from allele frequency list
            transN <- trans[,.N]
            allele_matrix <- matrix(nrow = transN, ncol = num_loci * 2)  #create empty matrix, 2 columns per locus (a1, a2)
            for (locus in 1:num_loci) {                             #for loop is faster than lapply in this case, as it operates column-wise. lapply would have function overhead for column and row.
              num_alleles <- length(allele_frequencies[[locus]]) 
              alleles <- 1:num_alleles  # Define possible allele values
              allele_matrix[, (2*locus - 1)] <- sample(alleles, transN, replace = TRUE, prob = allele_frequencies[[locus]])
              allele_matrix[, (2*locus)] <- sample(alleles, transN, replace = TRUE, prob = allele_frequencies[[locus]])
            }
            loci_df <- as.data.table(allele_matrix)
            setnames(loci_df, 
                     old = names(loci_df), 
                     new = unlist(lapply(1:num_loci, function(locus) c(paste0("locus", locus, "_a1"), paste0("locus", locus, "_a2")))))
            
            # Merge with population
            trans[, (names(loci_df)) := loci_df]
            
            if (season[t] == 4) { #this checks if season == 4, i.e. age = 5. preg female translocation
              num_births <- pmin(maxOffspring, round(rnorm(NROW(trans), mean = fecundityRateMean, sd = fecundityRateSD))) * 
                rbinom(NROW(trans), 1, breedingSuccess)                                                                      #calculate litter sizes, clamp to maxOffspring, then multiply by bernoulli trial 
              trans[, offspringM := rbinom(.N, num_births, sexRatio)]       #generate number of M/F in litters using Bernoulli trials
              trans[, offspringF := num_births - offspringM]  
              
              #create offspring dataframe
              if(sum(num_births) > 0){
                offspringM <- data.table(
                  id = 0, 
                  age = 1, #as ageing step has already taken place 
                  sex = "M", 
                  offspringM = 0,
                  offspringF = 0,
                  mother_id = rep(trans[,id], trans[,offspringM]))  # Repeat each mother ID for the number of male offspring
                
                offspringF <- data.table(
                  id = 0, 
                  age = 1, #as ageing step has already taken place  
                  sex = "F",  
                  offspringM = 0,
                  offspringF = 0,
                  mother_id = rep(trans[,id], trans[,offspringF]))  # Repeat each mother ID for the number of female offspring
                
                offspring <- rbind(offspringM, offspringF)      # Combine male and female offspring into one dataframe
                setorder(offspring, mother_id) #sort so litters sit together
                offspring[, id := trans[,max(id)+1]:trans[,max(id) + sum(num_births)]] #create unique ID
              }
              
              
              ##### Assign genetic information to offspring, from mother and from 1 father
              if (NROW(offspring) > 0) {
                #create breeding colony father genetics (1 per litter)
                mother_ids <- offspring[,mother_id]           #get mother id's
                unique_mothers <- unique(mother_ids)          # Get unique mothers
                
                # Create an empty matrix for father alleles (one per unique mother) and fill with alleles from allele_frequencies
                allele_matrix_fathers <- matrix(nrow = length(unique_mothers), ncol = num_loci * 2)           
                for (locus in 1:num_loci) {  
                  num_alleles <- length(allele_frequencies[[locus]])  
                  alleles <- 1:num_alleles  # Define possible allele values
                  allele_matrix_fathers[, (2*locus - 1)] <- sample(alleles, length(unique_mothers), replace = TRUE, prob = allele_frequencies[[locus]])
                  allele_matrix_fathers[, (2*locus)] <- sample(alleles, length(unique_mothers), replace = TRUE, prob = allele_frequencies[[locus]])
                }
                
                father_genetics <- as.data.table(allele_matrix_fathers)
                setnames(father_genetics, 
                         old = names(father_genetics), 
                         new = unlist(lapply(1:num_loci, function(locus) c(paste0("locus", locus, "_a1"), paste0("locus", locus, "_a2")))))
                
                father_loci_df <- father_genetics[match(mother_ids, unique_mothers), ]           # Expand father genetics to match all offspring. Each offspring gets the same father genetics as their littermates
                
                
                #get mother and father alleles
                # Generate random allele selections for mothers and fathers
                mother_allele_choices <- matrix(sample(1:2, length(mother_ids) * num_loci, replace = TRUE), nrow = length(mother_ids), ncol = num_loci)
                father_allele_choices <- matrix(sample(1:2, length(father_loci_df) * num_loci, replace = TRUE), nrow = length(father_loci_df), ncol = num_loci)
                
                # Get mother and father rows using fast join
                mothers <- trans[J(mother_ids), on = "id"]
                fathers <- father_loci_df
                
                # Define column names for loci
                locus_columns_a1 <- paste0("locus", 1:num_loci, "_a1")
                locus_columns_a2 <- paste0("locus", 1:num_loci, "_a2")
                
                # Extract all allele values into separate data.tables
                mother_alleles_a1 <- mothers[, ..locus_columns_a1]   # Mothers' first allele
                mother_alleles_a2 <- mothers[, ..locus_columns_a2]   # Mothers' second allele
                father_alleles_a1 <- fathers[, ..locus_columns_a1]   # Fathers' first allele
                father_alleles_a2 <- fathers[, ..locus_columns_a2]   # Fathers' second allele
                
                # Select alleles dynamically based on earlier sample (1 = _a1, 2 = _a2)
                loci_inheritance_mother <- mother_alleles_a1 * (mother_allele_choices == 1) + 
                  mother_alleles_a2 * (mother_allele_choices == 2)
                loci_inheritance_father <- father_alleles_a1 * (father_allele_choices == 1) + 
                  father_alleles_a2 * (father_allele_choices == 2)
                
                # Convert results to data.table
                setnames(loci_inheritance_father, 
                         old = names(loci_inheritance_father), 
                         new = gsub("_a1$", "_a2", names(loci_inheritance_father)))
                #loci_inheritance <- cbind(loci_inheritance_mother, loci_inheritance_father)
                loci_inheritance <- loci_inheritance_mother[, names(loci_inheritance_father) := loci_inheritance_father]
                
                # sort columns
                col_order <- names(loci_inheritance)             # Get the column names
                locus_numbers <- as.numeric(sub("locus(\\d+)_a[12]", "\\1", col_order))  # Extracts numbers
                allele_suffix <- sub(".*_(a[12])", "\\1", col_order)  # Extracts a1 or a2
                sorted_order <- order(locus_numbers, allele_suffix)             # Order by locus number, then by suffix (a1 before a2)
                setcolorder(loci_inheritance, col_order[sorted_order]) # Reorder columns
                
                offspring[,(names(loci_inheritance)) := loci_inheritance]
                trans <- rbind(trans, offspring)
                
                mother_ids <- offspring$mother_id
              }
            }
            pop <- rbind(pop, trans) #append to population
          }
        }
      }
      
      
      ###################################
      # Enforce carrying capacity 
      ###################################
      #If the population exceeds capacity, reduce it via probabilistic truncation, as in Vortex
      #truncation in this way is probably realistic for phascogales, which are territorial. This basically simulates density dependent dispersal outside of the study area.
      
      populationSize <- NROW(pop)
      if (populationSize > K) {
        excess <- populationSize - K
        adjustmentFactor <- K / populationSize  # 1 / (1 + (excess / totalPopulation))
        
        truncation <- runif(NROW(pop)) < adjustmentFactor
        pop <- pop[truncation, ]
      }
      
      
      ###################################
      #save time-step information
      ###################################
      populationHistory[t] <- populationSize  # Calculate the total population size
      
      # Store results of population for each year, for function testing
      #results$Generation[(sim - 1) * time_steps + t] <- t
      #results$PopulationSize[(sim - 1) * time_steps + t] <- populationSize
      #results$Simulation[(sim - 1) * time_steps + t] <- sim
      #results$season[(sim - 1) * time_steps + t] <- season[t]
      #results$tHe[(sim - 1) * time_steps + t] <- pop.He.func(pop)
      
      
      ####################################
      # Check pop extant
      ####################################
      males[t] <- pop[sex == "M", .N]
      females[t] <- pop[sex == "F", .N]
      if (males[t] == 0 | females[t] == 0) {
        break
      } 
      
    }
    
    ########################################
    # Add values to run storage matrix
    ########################################
    
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
    res.mat[sim, 'Extant'] <- ifelse((males[time_steps] == 0 | females[time_steps] == 0), 0, 1)
    res.mat[sim, 'Nmin'] <- min(populationHistory)
    res.mat[sim, 'N5yr'] <- populationHistory[24]
    res.mat[sim, 'lambda'] <- lambda
    res.mat[sim, 'r'] <- r
    res.mat[sim, 'finalHe'] <- pop.He.func(pop) 
    
  }
  # Return the results
  #print(res.mat)
  return(res.mat)
  #return(results)
  
} #RTP.IBM.func wrapper

finalpop <- RTP.IBM.func(N=100, initialClass = 4, allele_frequencies = allele_frequencies)





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
  
  ##############################################################
  ## Run the RTP projection function for each parameter set
  ##############################################################
  
  result <- foreach(row = 1:nrow(samples), .combine = rbind, .packages=c('MASS', 'data.table'), .export = c('RTP.IBM.func', 'allele_frequencies', 'He.func', 'pop.He.func')) %dopar% {
    
    # Run the RTP simulation
    res.mat <- RTP.IBM.func(
      N = samples$nInitial[row],
      initialClass = samples$initialClass[row],
      time_steps = samples$nYears[row],
      releasedSexRatio = samples$releasedSexRatio[row],
      sexRatio = samples$sexRatio[row],
      fecundityRateMean = samples$fecundityRateMean[row],
      fecundityRateSD = samples$fecundityRateSD[row],
      maxOffspring = samples$maxOffspring[row],
      breedingSuccess = samples$breedingSuccess[row],
      environmentStochasticity = samples$environmentStochasticity[row],
      K = samples$carryingCapacity[row],
      maleSR = c(samples$maleDepJuvenileSR[row], 
                 samples$maleJuvenileSR[row],
                 samples$maleSubadultSR[row],
                 samples$maleAdultSR[row]),
      femaleSR = c(samples$femaleDepJuvenileSR[row],
                   samples$femaleJuvenileSR[row],
                   samples$femaleSubadultSR[row],
                   samples$femaleAdultSR[row]),
      translocations = list(
        list(tYear = samples$tYear1[row], numIndividuals = samples$nYear1[row]),
        list(tYear = samples$tYear2[row], numIndividuals = samples$nYear2[row])),
      translocatedSexRatio = samples$translocatedSexRatio[row],
      numSimulations = samples$numSimulations[row],
      allele_frequencies = allele_frequencies #allele freq is constant throughout sensitivity analysis
    )
    
    # Return the results
    return(res.mat)
  }
  
  ## return the results
  return(result)
}


####################################
#Run the model
####################################

#import allele frequencies and convert to list
ONE_pop_n.loci_Vortex1000L_new <- read.table("~/Research/Bush Heritage/PVA_secret-rocks-RTP/data/ONE_pop_n.loci_Vortex1000L_new.txt", quote="\"", comment.char="")
allele_frequencies <- lapply(split(ONE_pop_n.loci_Vortex1000L_new, seq(nrow(ONE_pop_n.loci_Vortex1000L_new))), as.double)

## Set parameter defaults and ranges for testing with sensitivity analyses
## Parameter defaults 
parDefault.list <- list()
parDefault.list$nYears <- 80
parDefault.list$nInitial <- 40
parDefault.list$initialClass <- 4
parDefault.list$releasedSexRatio <- 0.5
parDefault.list$translocatedSexRatio <- 0.5 
parDefault.list$tYear1 <- 3 
parDefault.list$nYear1 <- 20
parDefault.list$tYear2 <- 7 
parDefault.list$nYear2 <- 20
parDefault.list$maleDepJuvenileSR <- 0.6
parDefault.list$maleJuvenileSR <- 0.6
parDefault.list$maleSubadultSR <- 0.6
parDefault.list$maleAdultSR <- 0
parDefault.list$femaleDepJuvenileSR <- 0.6
parDefault.list$femaleJuvenileSR <- 0.6
parDefault.list$femaleSubadultSR <- 0.6
parDefault.list$femaleAdultSR <- 0.8
parDefault.list$sexRatio <- 0.62
parDefault.list$fecundityRateMean <- 7.5
parDefault.list$fecundityRateSD <- 0.69
parDefault.list$maxOffspring <- 8
parDefault.list$breedingSuccess <- 1
parDefault.list$environmentStochasticity <- 0.1
parDefault.list$carryingCapacity <- 500
parDefault.list$numSimulations <- 1


## Parameter ranges that are varied in sensitivity analysis
parRange.list <- list()
parRange.list$nInitial <- c(10,200)
parRange.list$initialClass <- c(2,5)
parRange.list$releasedSexRatio <- c(0.35,0.65)
parRange.list$translocatedSexRatio <- c(0.35,0.65) 
parRange.list$tYear1 <- c(1,4) 
parRange.list$nYear1 <- c(5,100)
parRange.list$tYear2 <- c(5,8) 
parRange.list$nYear2 <- c(5,100)
parRange.list$maleDepJuvenileSR <- c(0.2,1)
parRange.list$maleJuvenileSR <- c(0.2,1)
parRange.list$maleSubadultSR <- c(0.2,1)
parRange.list$maleAdultSR <- c(0,0)
parRange.list$femaleDepJuvenileSR <- c(0.2,1)
parRange.list$femaleJuvenileSR <- c(0.2,1)
parRange.list$femaleSubadultSR <- c(0.2,1)
parRange.list$femaleAdultSR <- c(0.2,1)
parRange.list$breedingSuccess <- c(0.2,1)
parRange.list$environmentStochasticity <- c(0.01,0.2)

## Set up parallel processing
## This will reduce simulation time and emulation time
## nproc is the number of processing cores you want to use

nproc <- 3
cl.tmp = makeCluster(rep('localhost',nproc), type='SOCK')
registerDoSNOW(cl.tmp)
getDoParWorkers()


## Sensitivity analysis

#Parameters to include:
SAvars <- c('nInitial','initialClass','releasedSexRatio',
            'translocatedSexRatio','tYear1','nYear1','tYear2', 'nYear2', 
            'maleDepJuvenileSR', 'maleJuvenileSR', 'maleSubadultSR', 'maleAdultSR', 'femaleDepJuvenileSR', 'femaleJuvenileSR', 'femaleSubadultSR', 'femaleAdultSR', 
            'breedingSuccess', 'environmentStochasticity')

## number of simulations (nSims) and number of parameter samples (nSamples)
## note that this code assumes you are running a single simulation iteration per parameter sample (as recommended in Toms paper)
## hence nSims = nSamples
nSims <- 100
nSamples <- nSims

## run the global sensitivity analysis
sa <- SA.func(SAvars=SAvars, nSims=nSims, nSamples=nSamples, type='latin')

#add a categorical variable for translocation strategy to sa
sa <- as.data.frame(sa)  # Convert matrix to data frame so you can add a factor
get_age_class <- function(x) {ifelse(x == 2, "juv", ifelse(x == 3, "sub", ifelse(x == 4, "adult", ifelse(x == 5, "pregF", "NA"))))}
sa <- cbind(sa, transStrat = paste(
  get_age_class(sa[, "initialClass"]), ".init_", 
  ifelse(sa[, "nYear1"] > 0, get_age_class(sa[, "tYear1"]), "NA"), ".y1_",
  ifelse(sa[, "nYear2"] > 0, get_age_class(sa[, "tYear2"]), "NA"), ".y2", 
  sep = ""
))
sa$transStrat <- as.factor(sa$transStrat)

#respecify SAvars with transStrat as a categorical variable, excluding intialClass, tYear1 and tYear2. i.e. only one variable describing translocation strategy
SAvars <- c('nInitial','releasedSexRatio',
            'translocatedSexRatio','nYear1', 'nYear2', 'transStrat',
            'maleDepJuvenileSR', 'maleJuvenileSR', 'maleSubadultSR', 'maleAdultSR', 'femaleDepJuvenileSR', 'femaleJuvenileSR', 'femaleSubadultSR', 'femaleAdultSR', 
            'breedingSuccess', 'environmentStochasticity')

