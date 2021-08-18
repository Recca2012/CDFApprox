############################################################
###Example file for CDF estimation in parallel
############################################################


############################################################
###Requiring packages#######################################
############################################################

###Package used to calculate the multivariate Normal distribution using Quasi Monte Carlo
require(mvPot)
###Package used to run functions in parallel
require(doParallel)
require(foreach)
require(snow)
require("Rmpi",lib="/storage/work/mfn120/") 
require(doMPI)
###Package used to manipulate strings
require(stringr)
###Internal functions used to calculate the likelihood in parallel
require(CDFNormalAproxPackNoC)

############################################################
###Defining areas locations and randmoly rearrange them#####
############################################################
A <- expand.grid(1:100,1:100)
A <- A[sample(nrow(A)),]

############################################################
###Calculating covariance and correlation matrix using######
############################################################
###Exponential Function#####################################
############################################################
###Creating distance matrix
D <- as.matrix(dist(A,upper = T,diag = T))
###Defining the variance parameter
sigma <- 1
###Defining the correlation parameter
rho <- 10
###Creating the covariance matrix using the exponential function
Sigma <- sigma^2*exp(-D/rho)
###Creating the correlation matrix, faster to create outside and feed to the function
Corr <- cov2cor(Sigma)

############################################################
###Creating mean vector#####################################
############################################################
M <- rep(0,nrow(A))

############################################################
###Creating upper limit to calculate the multivariate CDF###
############################################################
X <- M+1

############################################################
###Creating lattice used by the QuasiMonteCarlo approach####
############################################################

###Full distribution
p <- 3607
latticeRule <- genVecQMC(p, (nrow(A) ))
p <- 499
latticeRule1 <- genVecQMC(p, (nrow(A) ))
###Using 50 neighbors sequentially "R"
joint=50
use <- 50
latticelist50 <- lapply(X = 2:(use*joint+joint),FUN = function(x){
  genVecQMC(p, (x))})
###Using 50 neighbors sequentially "C"
latticeVec <- lapply(X = latticelist50,FUN = function(y){
  y$genVec
})
latticePrime<-sapply(X = latticelist50,FUN = function(y){
  y$primeP
},simplify = TRUE)

############################################################
###Running functions########################################
############################################################

###Creating array to save the final result
ResultsArray <- array(data = NA_real_,dim = c(7,7,15),dimnames = list(paste0("Neighbor",str_pad(c(2,5,10,20,30,50,100),3,pad="0")),#Neighbors are the number of neighbors used to estimate the joint distribution
                                                                      paste0("Joint",str_pad(c(2,5,10,20,30,50,100),3,pad="0")),#Join is the number of observations being jointly estimated
                                                                      paste0("Int",str_pad(c(1:15),2,pad="0"))))#Number of independent interations being used
###Creating array to save the running time
TimeArray <- array(data = NA_real_,dim = c(7,7,15),dimnames = list(paste0("Neighbor",str_pad(c(2,5,10,20,30,50,100),3,pad="0")),
                                                                      paste0("Joint",str_pad(c(2,5,10,20,30,50,100),3,pad="0")),
                                                                      paste0("Int",str_pad(c(1:15),2,pad="0"))))
############################################################
###Defining parameters used by the functions################
############################################################

###Upper value used to define the CDF
upper=X
###Defining the number of cores being used, one less than the number of available cores
cores=39
###Defining the lattice used for the multivariate Normal distribution
lattice = latticelist50
###Starting the cluster using MPI, there are other functions that do not use MPI
cl <- startMPIcluster(count = 40,verbose = TRUE)
registerDoMPI(cl = cl)

###Running the example 15 independent times
for(i in 1:15){
  ###Rearrange locations in other to estimate the same CDF with possibly different observations
  pos <- sample(nrow(Sigma))
  Sigma <- Sigma[pos,pos]
  Corr <- Corr[pos,pos]
  #Joint estimation using 2 observations at a time
  TimeArray[1,1,i]<-system.time(ResultsArray[1,1,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 002,cores = cores,
                                                              lattice = lattice,joint = 002))[[3]]
  TimeArray[2,1,i]<-system.time(ResultsArray[2,1,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 005,cores = cores,
                                                              lattice = lattice,joint = 002))[[3]]
  TimeArray[3,1,i]<-system.time(ResultsArray[3,1,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 010,cores = cores,
                                                              lattice = lattice,joint = 002))[[3]]
  TimeArray[4,1,i]<-system.time(ResultsArray[4,1,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 020,cores = cores,
                                                              lattice = lattice,joint = 002))[[3]]
  TimeArray[5,1,i]<-system.time(ResultsArray[5,1,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 030,cores = cores,
                                                              lattice = lattice,joint = 002))[[3]]
  TimeArray[6,1,i]<-system.time(ResultsArray[6,1,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 050,cores = cores,
                                                              lattice = lattice,joint = 002))[[3]]
  
  #Joint estimation using 5 observations at a time
  TimeArray[1,2,i]<-system.time(ResultsArray[1,2,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 002,cores = cores,
                                                              lattice = lattice,joint = 005))[[3]]
  TimeArray[2,2,i]<-system.time(ResultsArray[2,2,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 005,cores = cores,
                                                              lattice = lattice,joint = 005))[[3]]
  TimeArray[3,2,i]<-system.time(ResultsArray[3,2,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 010,cores = cores,
                                                              lattice = lattice,joint = 005))[[3]]
  TimeArray[4,2,i]<-system.time(ResultsArray[4,2,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 020,cores = cores,
                                                              lattice = lattice,joint = 005))[[3]]
  TimeArray[5,2,i]<-system.time(ResultsArray[5,2,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 030,cores = cores,
                                                              lattice = lattice,joint = 005))[[3]]
  TimeArray[6,2,i]<-system.time(ResultsArray[6,2,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 050,cores = cores,
                                                              lattice = lattice,joint = 005))[[3]]
  
  #Joint estimation using 10 observations at a time
  TimeArray[1,3,i]<-system.time(ResultsArray[1,3,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 002,cores = cores,
                                                              lattice = lattice,joint = 010))[[3]]
  TimeArray[2,3,i]<-system.time(ResultsArray[2,3,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 005,cores = cores,
                                                              lattice = lattice,joint = 010))[[3]]
  TimeArray[3,3,i]<-system.time(ResultsArray[3,3,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 010,cores = cores,
                                                              lattice = lattice,joint = 010))[[3]]
  TimeArray[4,3,i]<-system.time(ResultsArray[4,3,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 020,cores = cores,
                                                              lattice = lattice,joint = 010))[[3]]
  TimeArray[5,3,i]<-system.time(ResultsArray[5,3,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 030,cores = cores,
                                                              lattice = lattice,joint = 010))[[3]]
  TimeArray[6,3,i]<-system.time(ResultsArray[6,3,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 050,cores = cores,
                                                              lattice = lattice,joint = 010))[[3]]
  
  #Joint estimation using 20 observations at a time
  TimeArray[1,4,i]<-system.time(ResultsArray[1,4,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 002,cores = cores,
                                                              lattice = lattice,joint = 020))[[3]]
  TimeArray[2,4,i]<-system.time(ResultsArray[2,4,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 005,cores = cores,
                                                              lattice = lattice,joint = 020))[[3]]
  TimeArray[3,4,i]<-system.time(ResultsArray[3,4,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 010,cores = cores,
                                                              lattice = lattice,joint = 020))[[3]]
  TimeArray[4,4,i]<-system.time(ResultsArray[4,4,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 020,cores = cores,
                                                              lattice = lattice,joint = 020))[[3]]
  TimeArray[5,4,i]<-system.time(ResultsArray[5,4,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 030,cores = cores,
                                                              lattice = lattice,joint = 020))[[3]]
  TimeArray[6,4,i]<-system.time(ResultsArray[6,4,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 050,cores = cores,
                                                              lattice = lattice,joint = 020))[[3]]
  
  #Joint estimation using 30 observations at a time
  TimeArray[1,5,i]<-system.time(ResultsArray[1,5,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 002,cores = cores,
                                                              lattice = lattice,joint = 030))[[3]]
  TimeArray[2,5,i]<-system.time(ResultsArray[2,5,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 005,cores = cores,
                                                              lattice = lattice,joint = 030))[[3]]
  TimeArray[3,5,i]<-system.time(ResultsArray[3,5,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 010,cores = cores,
                                                              lattice = lattice,joint = 030))[[3]]
  TimeArray[4,5,i]<-system.time(ResultsArray[4,5,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 020,cores = cores,
                                                              lattice = lattice,joint = 030))[[3]]
  TimeArray[5,5,i]<-system.time(ResultsArray[5,5,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 030,cores = cores,
                                                              lattice = lattice,joint = 030))[[3]]
  TimeArray[6,5,i]<-system.time(ResultsArray[6,5,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 050,cores = cores,
                                                              lattice = lattice,joint = 030))[[3]]
  
  #Joint estimation using 50 observations at a time
  TimeArray[1,6,i]<-system.time(ResultsArray[1,6,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 002,cores = cores,
                                                              lattice = lattice,joint = 050))[[3]]
  TimeArray[2,6,i]<-system.time(ResultsArray[2,6,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 005,cores = cores,
                                                              lattice = lattice,joint = 050))[[3]]
  TimeArray[3,6,i]<-system.time(ResultsArray[3,6,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 010,cores = cores,
                                                              lattice = lattice,joint = 050))[[3]]
  TimeArray[4,6,i]<-system.time(ResultsArray[4,6,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 020,cores = cores,
                                                              lattice = lattice,joint = 050))[[3]]
  TimeArray[5,6,i]<-system.time(ResultsArray[5,6,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 030,cores = cores,
                                                              lattice = lattice,joint = 050))[[3]]
  TimeArray[6,6,i]<-system.time(ResultsArray[6,6,i] <- pmvNorm4(M = M,Sigma = Sigma,Corr=Corr,upper = upper,
                                                              use = 050,cores = cores,
                                                              lattice = lattice,joint = 050))[[3]]
}

###Stopping cluster
closeCluster(cl)

###Saving results
save(ResultsArray,TimeArray,file = "result.RData")
