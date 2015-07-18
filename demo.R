###Here is a demonstration of the function
###	SEandBIAS
###
###As input, the function requires the arguments
###	sig.t = 	var(T)
###	sid.d = 	var(delta)
###	sig.u = 	var(U)
###	sig.e = 	var(epislon)
###	sig.r = 	var(r)
###	beta.1 = 	beta1
###	N	 =	Number of observations per individual
###	RHO	 =	covariance matrix of delta (I = Independent; T = Toeplitz; A = Autocorrelation)
###   rho_y  =  	correlation in RHO (when not independent)
###All Notation is consistent w/ Sampson et al. Journal of Applied Statistics.
###Output
###	3 x 3 matrix.
###	Rows correspond to the parameters beta.1, sig.e, and sig.r
###	Entries are the SE when the full model is fit, the bias when the reduced model is fit, the se when the reduced model is fit
###
###	This code requires the following four libraries 

	library(mvtnorm)
	library(stats4)
	library(lme4)
	library(rootSolve)

###EXAMPLE FROM SECTION 2.3
###source the functions5.R code included here
folder <- "H:/Projects/Matthews/estVar/"
source(paste(folder,"functions5.R",sep=""))

###Define the parameters; 
sig.t  <- 4
sig.d  <- 1.5
sig.u  <- 0.1
sig.e  <- 3.0
sig.r  <- 0.8
beta.1 <- 1.0
N      <- 4
nsub   <- 50
RHO 	 <- "I"
rho.y  <- 0.0

###Run the function
SEandBIAS(sig.t=sig.t,sig.d=sig.d,sig.u=sig.u,beta.1=beta.1,sig.r=sig.r,sig.e=sig.e,N=N,nsub=nsub,RHO=RHO,rho.y=rho.y)



