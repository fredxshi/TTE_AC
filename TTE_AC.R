library(R2OpenBUGS)
library(survival)
library(survminer)
library(gg.gap)
library(dplyr)

#Analyze the simulated data - good bump with a good bump
dd_ac <- read.csv("C:/fakepath/Simulated_Data.csv", header = T) 
dd <- dd_ac %>%
  mutate(arm = case_when(Dose==1 ~ 1, Dose==2 ~ 2, Dose==3 ~ 3, Dose==8 ~ 4,
                         Dose==4 ~ 5, Dose==5 ~ 6, Dose==6 ~ 7,Dose==7~ 8),
         arm3 = case_when(Dose==1 ~ 1, Dose==2 ~ 2, Dose==3 ~ 2, Dose== 8 ~ 3,
                          Dose==4 ~ 2, Dose==5 ~ 2, Dose==6 ~ 2, Dose==7~ 2), #arm of three-group model
         arm2 = case_when(Dose==1 ~ 1, Dose==2 ~ 2, Dose==3 ~ 2, Dose==8 ~ 2,
                          Dose==4 ~ 2, Dose==5 ~ 2, Dose==6 ~ 2, Dose==7~ 2))  #arm of two-group model

##Read simulated data from FACTS
data.list = list(N = nrow(dd_ac),
                 death = dd_ac$Outcome,
                 t = dd_ac$Duration,
                 dose.ac = dd$Dose,
                 dose = dd$arm,
                 dose2 = dd$arm2,
                 dose3 = dd$arm3)

##Create BUGS model
select.model <- function(model){
  cat(" model{
  c[1]<-0
  c[2]<-2
  c[3]<-52
  \n")
  
  #Select arm
  if ( model %in% c("EMAX","HEMAX", "INDEP") ) { cat("for (i in 1:N) { arm[i]<-dose[i] } \n") }
  else if ( model == "TWOGROUP" ) { cat("for (i in 1:N) { arm[i]<-dose2[i] } \n") }
  else if ( model == "THREEGROUP" ) { cat("for (i in 1:N) { arm[i]<-dose3[i] } \n") }
  else { cat("for (i in 1:N) { arm[i]<-dose.ac[i] } \n") }
  
  #Poisson process
  cat(" 
  for (i in 1:N) { #pt index
    for (k in 1:2){ #segment index
    event[i, k] <- (death[i]) * step(t[i] - c[k]) * step(c[k+1] - t[i]) #event indicator
    delta[i, k] <- (min(t[i], c[k+1]) - c[k]) * step(t[i] - c[k]) #time lived in this segment
    hazard[i, k] <- lambda[k] * exp(theta[arm[i]]) #hazard function
    }
  }
  for (i in 1:N) { #pt index
    for (k in 1:2){ #segment index
    #likelihood
    event[i, k] ~ dpois(mu[i, k])
    mu[i, k] <- delta[i,k] * hazard[i,k]
    }
  }

  theta[1] <- 0
  \n")
  
  if ( model %in% c("EMAX","HEMAX","EMAXAC","HEMAXAC") ) {
    if ( model %in% c("EMAXAC","HEMAXAC") ) { 
      cat("for (d in 2:7) {\n") 
    }
    else { 
      cat("for (d in 2:8) {\n")
    } 
    if ( model %in% c("HEMAX","HEMAXAC") ){
      cat("theta[d] <- a[1] + a[2] * nu[d] / (nu[d] + a[3]) + sig[d]\n")
      cat("sig[d] ~ dnorm(0, invtau2)\n") # prior of error term for HEMAX
    } else {
      cat("theta[d] <- a[1] + a[2] * nu[d] / (nu[d] + a[3])\n")
    }
    cat("}\n")
    if ( model %in% c("HEMAX","HEMAXAC") ) { cat("invtau2 ~ dgamma(0.1, 0.1)\n") } # hyperprior of error term for HEMAX
    if ( model %in% c("EMAXAC","HEMAXAC") ) { cat("theta[8] ~ dnorm(0, 0.01)\n") } # active comparator
  }
  else if (model == "TWOGROUP") { cat("theta[2] ~ dnorm(0,0.01)\n") }
  else if (model == "THREEGROUP") { cat("theta[2] ~ dnorm(0,0.01)
                                        theta[3] ~ dnorm(0,0.01)\n") }
  else { cat("for (d in 2:8) { theta[d] ~ dnorm(0,0.01) }\n") }
  
  if (model == "TWOGROUP") { cat("for (d in 1:2) {\n")}
  else if (model == "THREEGROUP") { cat("for (d in 1:3) {\n")}
  else { cat("for (d in 1:8) {\n") }
  
  cat("
    hr[d] <- exp(theta[d])
    for (k in 1:2){
      haz[d,k] <- lambda[k] * hr[d]
    }
    m1[d] <- -log(0.5) / haz[d,1]
    m2[d] <- (-log(0.5)-2*haz[d,1] + 2*haz[d,2])/haz[d,2]
    mediansurv[d] <- m1[d]*(1-step(m1[d]-2)) + m2[d]*(step(m1[d]-2)) #PE median ST
    rmst[d] <- (1-exp(-2*haz[d,1]))/haz[d,1] + (exp(-2*haz[d,1])-exp(-24*haz[d,2]-2*haz[d,1]))/haz[d,2] #RMST
    meansurv[d] <- 1/haz[d,1] + (1/haz[d,2]-1/haz[d,1])*exp(-2*haz[d,1]) #PE mean ST
  }
  \n")
  
  if ( model %in% c("EMAX","HEMAX","EMAXAC","HEMAXAC") ) {
    cat("
  nu[2]<-2.6
  nu[3]<-4.17
  \n")
    
    if ( model %in% c("EMAXAC","HEMAXAC") ) {
      cat("
  nu[4]<-5.92
  nu[5]<-6.2
  nu[6]<-7.76
  nu[7]<-9.52
  \n")
    }
    else {
      cat("
  nu[4]<-5.40
  nu[5]<-5.92
  nu[6]<-6.20
  nu[7]<-7.76
  nu[8]<-9.52
  \n")
    }
    
    cat("
  a[1]~dnorm(0,0.01) #1/10^2
  a[2]~dnorm(0,0.04) #1/5^2
  a[3]~dnorm(3,0.04)T(0,)
  \n")
  }
  
  cat("
  lambda[1]~dgamma(1,5)
  lambda[2]~dgamma(1,20)
  }\n") 
}

########################################################################################################################
sink("C:/fakepath/model.txt")
select.model(model="EMAXAC") #Select your dose response model: EMAX, EMAXAC, HEMAX, HEMAXAC, INDEP, TWOGROUP, THREEGROUP
sink()
########################################################################################################################

##Posteriors
set.seed(3500)
model.output <- bugs(data.list, inits = NULL, 
                     parameters.to.save = c("hr","lambda", "mediansurv", "rmst", "meansurv"), 
                     working.directory=" C:/fakepath/", model.file = "model.txt", 
                     n.chains = 1, n.iter = 55000, n.burnin = 5000, n.thin = 1, debug = TRUE)
round(model.output$summary,2)