#' Unfortunately, I am afraid that incident only cases cannot 
#' prevent Berkson bias.

#' Say individuals have 3 independent genetic variants 
#' that confere risk for 3 independent diseases.
#' Can we prevent berkson bias by looking only at incident
#' cases? That is, by only using the disease that first occured
#' can we control this bias? 
 
#' In this setting, two independent risk factors will become 
#' associated with each other by chance and a false positive
#' association with a disease will arise.

#' Genetic background of each indinvidual:
#' 
#' * N    - number of individuals
#' * Phi1 - frequency of mutation 1
#' * Phi2 - frequency of mutation 2
#' * Phi3 - frequency of mutation 3 
N    = 2000
Phi1 = .1
Phi2 = .2
Phi3 = .3
G1 = rbinom(n = N,prob = Phi1, size = 1)
G2 = rbinom(n = N,prob = Phi2, size = 1)
G3 = rbinom(n = N,prob = Phi3, size = 1)
Genotype = cbind(G1,G2,G3)


#' Life time and simulation time quanta
years_longevity = 1               # Assume in this world, individuals die at 1 year old
lifetime = 365 * years_longevity  # Using 24 time quanta for the simulation

#' Assume constant hazard rates, i.e. the probabilities of 
#' disease incidence in 24h is a constant, irrespective of age, etc:
#' 
#' * The probability of developing disease D1 given variant G1 is P1
#' * The probability of developing disease D2 given variant G1 is P2
#' * The probability of developing disease D3 given variant G1 is P3
#' 
#' Assume identical hazard rates for all diseases
P1 = .3 / lifetime
P2 = .3 / lifetime
P3 = .3 / lifetime


#' Keep a record of who went to the hospital with what disease
#' We assume that hospitalization rate is 100% and that hospitalization.
#' That is, disease occurance = hospitalization.
#' As berkson notes, if hospitalization is 100% bias
#' will be smaller.
#' 
#' We assume that hospitalization will cure the disease, but the 
#' disease may occur again, e.g. flu,  bronchitis, drug depency.
#' This doesent matter for the results because we only care 
#' about the first disease the patient ever had. Altough to be precise
#' the simulation is based on recurrent and remission disease episodes,
#' the actual treatment of the simulation makes the resultidentical for 
#' acute and chronic diseases.
hospital = matrix(data = NA,nrow =  N,ncol =  lifetime)
for (iday in 1:lifetime) {
  D1 = as.logical(rbinom(n = N,prob = P1, size = 1)*G1)
  D2 = as.logical(rbinom(n = N,prob = P2, size = 1)*G2)
  D3 = as.logical(rbinom(n = N,prob = P3, size = 1)*G3)
  state = as.data.frame(cbind(D1,D2,D3))
  
  hospital[,iday] = apply(X = state, MARGIN = 1, FUN = function(x) ifelse( sum(x) > 0, max(which(x)),0 ))
}

#' Hospitalizations per individual
hist(rowSums(hospital > 0,na.rm = T),xlab = "Hospitalizations/Individual",main = "")


#' Hospitalizations per day 
hist(colSums(hospital > 0,na.rm = T),xlab = "Hospitalizations/Day in the population",main = "")

#' Produce a disease variable using incident cases
#' i.e. select the disease that manifested first
disease = matrix(data = 0,nrow = N,ncol = 1)
for (i in 1:N) {
  case       = hospital[i,]
  disease[i] = case[which.max(case > 0)] #Retrieve only first incident case
}

#' Altough all genotypes were i.i.d there are cases of 
#' individuals that have more than 2 risk-confering mutations,
#' see below the number in this simulation:
sum(rowSums(Genotype) > 1)
#' What I submit to you is that even if you just have one disease in your lifetime and then die, 
#' the individuals that have risk factors for multiple diseases are more frequent among incident cases
#' and through this mechanism spurious associations between diseases and all independent risk factors
#' are created. This is the actual mechanism by which berkson bias occur in the hospitalized
#' population.It is not a matter of differential selection. Berkson paper uses constant selection
#' for all diseases, but some individuals will have by change more than one disease.
#' 
#' Notice that even if the risk factors are completelly independent from each other, some individuals will still have
#' multiple risk factors with probability P1 x P2 x P3.  
#' 
#' 
#' Proportion of hospitalizations among carriers of multiple risk factors
#' and single carriers
barplot( c(mean(disease[rowSums(Genotype) > 1] > 0), mean(disease[rowSums(Genotype) == 1] > 0)),names.arg = c("Carrier of multiple mutations","Carrier of single mutation"),ylab="Proportion of hospitalized individuals", col= "cyan" )


#' Let's analyse the whole population (with and without disease) and 
#' the hospitalized population. 
population_study = data.frame(disease = as.factor(disease), G1,G2,G3)
hospital_records = subset(population_study, disease != 0 ) 

results = 
  rbind(
    
#Whole population analysis

cbind(
  rbind(
    coef(summary(glm( I(disease==1) ~ G1,  data =  population_study)))[2,c(1,4)],
    coef(summary(glm( I(disease==1) ~ G2,  data =  population_study)))[2,c(1,4)],
    coef(summary(glm( I(disease==1) ~ G3,  data =  population_study)))[2,c(1,4)]
  ),
  
  rbind(
    coef(summary(glm( I(disease==2) ~ G1, data =  population_study)))[2,c(1,4)],
    coef(summary(glm( I(disease==2) ~ G2, data =  population_study)))[2,c(1,4)],
    coef(summary(glm( I(disease==2) ~ G3, data =  population_study)))[2,c(1,4)]
  ),
  
  rbind(
    coef(summary(glm( I(disease==3) ~ G1, data =  population_study)))[2,c(1,4)],
    coef(summary(glm( I(disease==3) ~ G2, data =  population_study)))[2,c(1,4)],
    coef(summary(glm( I(disease==3) ~ G3, data =  population_study)))[2,c(1,4)]
  )
),

#Hospitalized population analysis

cbind(
  rbind(
    coef(summary(glm( I(disease==1) ~ G1, data =  hospital_records)))[2,c(1,4)],
    coef(summary(glm( I(disease==1) ~ G2, data =  hospital_records)))[2,c(1,4)],
    coef(summary(glm( I(disease==1) ~ G3, data =  hospital_records)))[2,c(1,4)]
  ),
  
  rbind(
    coef(summary(glm( I(disease==2) ~ G1, data =  hospital_records)))[2,c(1,4)],
    coef(summary(glm( I(disease==2) ~ G2, data =  hospital_records)))[2,c(1,4)],
    coef(summary(glm( I(disease==2) ~ G3, data =  hospital_records)))[2,c(1,4)]
  ),
  
  rbind(
    coef(summary(glm( I(disease==3) ~ G1, data =  hospital_records)))[2,c(1,4)],
    coef(summary(glm( I(disease==3) ~ G2, data =  hospital_records)))[2,c(1,4)],
    coef(summary(glm( I(disease==3) ~ G3, data =  hospital_records)))[2,c(1,4)]
  )
)
)

rownames(results) = c( paste0( rep("Population D=",3),1:3),
                       paste0( rep("Hospital D=",3),1:3) )

P_results = results[,seq(2,ncol(results),2)]
colnames(P_results) = paste0( rep("G",3),1:3 )

#' The population data only shows associations between disease and their respective risk factor. 
#' On the other hand, the risk factors are associated with all diseases 
#' in the hospitalized population 
print(P_results)
