###All notation corresponds to summary.doc

###Calculates the trace of a matrix

tr    <- trace <- function(x) sum(diag(x))

####calculates the true correlation matrix of Tij (deriv = 0)
#### rho can be either 
#### "I": independent, 		\rho_{ij} = 0 for all i,j
#### "T": Toeplitz,    		\rho_{ij} = rho.y if |i-j| = 1
#### "A": Autocorrelation,    \rho_{ij} = rho.y^|i-j| 
#### or calculates the derivative of the true correlation matrix (relative to rho.y) if deriv = 1

calcT2 <- function(RHO="I",rho.y=0,nops=3,deriv=0){

	if (deriv == 0){
			  T2 <- diag(1,nrow=nops)
	if (RHO=="T") for (i in 1:nops) for (j in 1:nops) if (abs(i-j)==1) T2[i,j] <- rho.y
	if (RHO=="A") for (i in 1:nops) for (j in 1:nops) 	 		 T2[i,j] <- rho.y^abs(i-j)
	}

	if (deriv == 1){
			  T2 <- diag(0,nrow=nops)
	if (RHO=="T") for (i in 1:nops) for (j in 1:nops) if (abs(i-j)==1) T2[i,j] <- 1
	if (RHO=="A") for (i in 1:nops) for (j in 1:nops) 	 		 T2[i,j] <- max(abs(i-j),0)*rho.y^max(0,abs(i-j)-1)
	}
T2
}

### calculates the covariance matrix for the vector of c(Q,S) assuming a given correlation structure, RHO
### theta1:	\sigma^2_t + \sigma^2_{\gamma}
### sig.t:	\sigma^2_t
### sig.g:	\sigma^2_{\delta}
### beta.1:	beta_1
### sig.r:	\sigma^2_R
### sig.e:	\sigma^2_{\epsilon}
### nops:	number of observations per subject

calcV <- function(theta1=1,sig.t=0,sig.g=1,beta.1=1,sig.r=0,sig.e=1,nops=2,RHO="I",rho.y=0){
	T2	    <- calcT2(RHO,rho.y,nops)

	vQ        <- (beta.1^2)*theta1 + sig.r + sig.e
	vQ12      <- sig.r + (beta.1^2)*(sig.t+(theta1-sig.t)*T2)
	vS        <- theta1 + sig.g
	vS12      <- sig.t + (theta1-sig.t)*T2
	vQS       <- beta.1*theta1
	vQ1S2     <- beta.1*(sig.t + (theta1-sig.t)*T2)

	V11       <- vQ12
	V12       <- vQ1S2
	V22       <- vS12
	diag(V11) <- vQ
	diag(V22) <- vS
	diag(V12) <- vQS
	V         <- rbind(cbind(V11,V12),cbind(V12,V22))
V
}


###	calculates the true information matrix for the true parameters when Q and S 
###	Parameters are divided into groups: those related to the variance and those related to the mean
###   Parameters related to variance:	dnames = c("theta1","sig.e","sig.t","sig.g","beta.1","sig.r")
###		Note: any of these parameters can be left of the list if they are presumed to be known
###		(e.g. if sig.r = 0, or R is to be left of out of the model)
###		the vector d indicates which variables (as ordered above) are desired
### 	Parameters related to the mean are only included if mu = 1
###	The Information matrix is produced for the desired parameters (based on dnames and mu)
###	calculation based on theory given in wikipedia 
###		http://en.wikipedia.org/wiki/Fisher_information#Multivariate_normal_distribution

calcI  <- function(theta1=1,sig.t=0,sig.g=1,beta.1=1,sig.r=0,sig.e=1,d=c(1,2,3,4),nops=2,RHO="I",rho.y=0,mu=0){

	##	a list containing 7 matrices. Each matrix is a term x term derivative of the variance matrix for c(Q,S)
	##	matrices are ordered to correspond to dnames

	derivs <- list()
	dnames <- c("sig.y","sig.e","sig.t","sig.g","beta.1","sig.r","rho.y")

	T2	    <- calcT2(RHO,rho.y,nops)
	T2d	    <- calcT2(RHO,rho.y,nops,deriv=1)
	M1	    <- matrix(1,nrow=nops,ncol=nops)
	MI	    <- diag(1,nrow=nops)
	sig.y	    <- theta1-sig.t

	for (i in 1:7){

	V11 <- V12 <- V22 <- matrix(0,nrow=nops,ncol=nops)

	if (i==1) { V11 <- (beta.1^2)*T2
                  V12 <- beta.1*T2
			V22 <- T2}

	if (i==2)   V11 <- MI

	if (i==3) { V11 <- (beta.1^2)*M1
			V12 <- beta.1*M1
			V22 <- M1}

	if (i==4)   V22 <- MI

	if (i==5)  {V11 <- 2*beta.1*(sig.t+sig.y*T2)
			V12 <- (sig.t+sig.y*T2)}

	if (i==6)   V11 <- M1
	
	if (i==7)	{V11 <- (beta.1^2)*sig.y*T2d
			 V12 <- beta.1*sig.y*T2d
			 V22 <- sig.y*T2d}
		
		V         		 <- rbind(cbind(V11,V12),cbind(V12,V22))
		derivs[[i]] 	 <- V
	}

	names(derivs) <- dnames

	V   <- calcV(theta1,sig.t,sig.g,beta.1,sig.r,sig.e,nops,RHO,rho.y)
	V.i <- solve(V)
	ld  <- length(d)

	###I.mat is the information matrix
	I.mat <- matrix(nrow=ld,ncol=ld)
	###each term is calculated from the equations on the wikipedia page
	for (i in 1:ld) for (j in 1:ld) I.mat[i,j] <- 0.5*sum(diag(V.i%*%derivs[[d[i]]]%*%V.i%*%derivs[[d[j]]]))
	rownames(I.mat) <- dnames[d]
	I.out <- I.mat

	###If we need to include the intercepts (interest estimates are independent of all other parameters)
	du.dbq <- matrix(rep(c(1,0),each=nops),ncol=1)
	du.dbs <- matrix(rep(c(0,1),each=nops),ncol=1)
	I.mat.w.mean <- matrix(0,nrow=(ld+2),ncol=(ld+2))
	I.mat.w.mean[1:ld,1:ld] <- I.mat
	I.mat.w.mean[ld+1,ld+1] <- t(du.dbq)%*%V.i%*%du.dbq
	I.mat.w.mean[ld+2,ld+2] <- t(du.dbs)%*%V.i%*%du.dbs
	I.mat.w.mean[ld+1,ld+2] <- I.mat.w.mean[ld+2,ld+1] <- t(du.dbq)%*%V.i%*%du.dbs
	rownames(I.mat.w.mean)  <- c(dnames[d],"mean.S","mean.Q")
	if (mu == 1) I.out <- I.mat.w.mean

I.out
}

###	Calculates the standard error for the MLE estimates of the following parameters:
###	c("theta1","sig.e","sig.t","sig.g","beta.1","sig.r")[d] [add on estimates of E[Q] and E[S] if mu=1

calcSE  <- function(theta1=1,sig.t=0,sig.g=1,beta.1=1,sig.r=0,sig.e=1,d=c(1,2,3,4),nops=2,nsub=100,RHO="I",rho.y=0,mu=0,wantVar=0){
	ld	  <- length(d)
	I.mat   <- calcI(theta1,sig.t,sig.g,beta.1,sig.r,sig.e,d,nops,RHO,rho.y,mu)
	se      <- matrix(sqrt(diag(solve(nsub*(I.mat)))),ncol=1)
	rownames(se) <- rownames(I.mat)
	if (wantVar==0) out <- se
	if (wantVar==1) out <- solve(nsub*(I.mat))
out
}

####
####	Fitting the reduced model (where we assume that S is the truth)
####


###	Calculate the variance of the derivative of the log-likelihood at the limits of the naive solutions
###	(e.g. the var(\ell' (\Theta^{\dagger}))$ 
###	the order of the derivatives: \beta_0, \beta_1, \sigma^2_R, \sigma^2_{\epsilon}

calc.Vscore	  <- function(nops,VWZ,V.naive.i,nsub){

	###V.star is shorthand notation for var( \sigma_{red}^{-1} (Q - \beta_0 - S\beta_1) )

	V.star   <- VWZ[(nops+1):(2*nops),(nops+1):(2*nops)]

	###derivates of $\Sigma_{red}$ relative to $\sigma^2_R$ and $\Sigma^2_{\epsilon}
	dvdr	   <- matrix(1,nrow=nops,ncol=nops)
	dvde	   <- diag(1,nrow=nops)

	V.dlogLik <- matrix(0,nrow=4,ncol=4)
	W0	    <- matrix(1,nrow=nops,ncol=1)
	V.dlogLik[1,1] <- t(W0)%*%V.star%*%W0
	V.dlogLik.22   <- 0

	for (i in 1:nops) for (j in 1:nops) {
		V.dlogLik.22 <- V.dlogLik.22 + VWZ[i,i+nops]*VWZ[j,j+nops]+VWZ[i,j]*VWZ[i+nops,j+nops]+VWZ[i,j+nops]*VWZ[j,i+nops]
	}
	VWZ.od	   <- 0
	for (j in 1:nops) VWZ.od <- VWZ[j,j+nops]
	V.dlogLik[2,2] <- V.dlogLik.22 - (VWZ.od)^2


	V.dlogLik.33   <- V.dlogLik.44   <- V.dlogLik.34   <- 0
	for (i in 1:nops) for (j in 1:nops) for (k in 1:nops) for (l in 1:nops) {
		V.dlogLik.33 <- V.dlogLik.33 + (V.star[i,j]*V.star[k,l] + V.star[i,k]*V.star[j,l]+ V.star[i,l]*V.star[j,k])*dvdr[i,j]*dvdr[k,l]
		V.dlogLik.44 <- V.dlogLik.44 + (V.star[i,j]*V.star[k,l] + V.star[i,k]*V.star[j,l]+ V.star[i,l]*V.star[j,k])*dvde[i,j]*dvde[k,l]
		V.dlogLik.34 <- V.dlogLik.34 + (V.star[i,j]*V.star[k,l] + V.star[i,k]*V.star[j,l]+ V.star[i,l]*V.star[j,k])*dvdr[i,j]*dvde[k,l]
	}

	V.dlogLik[3,3] <- V.dlogLik.33 - (sum(V.star*dvdr))^2
	V.dlogLik[4,4] <- V.dlogLik.44 - (sum(V.star*dvde))^2
	V.dlogLik[3,4] <- V.dlogLik[4,3] <- V.dlogLik.34 - sum(V.star*dvde)*sum(V.star*dvdr)

	V.dlogLik.23   <- V.dlogLik.24   <- 0
	for (i in 1:nops) for (k in 1:nops) for (l in 1:nops) {
		V.dlogLik.23 <- V.dlogLik.23 + (VWZ[i,i+nops]*VWZ[k+nops,l+nops] + VWZ[i,k+nops]*VWZ[i+nops,l+nops]+ VWZ[i,l+nops]*VWZ[i+nops,k+nops])*dvdr[k,l]
		V.dlogLik.24 <- V.dlogLik.24 + (VWZ[i,i+nops]*VWZ[k+nops,l+nops] + VWZ[i,k+nops]*VWZ[i+nops,l+nops]+ VWZ[i,l+nops]*VWZ[i+nops,k+nops])*dvde[k,l]
	}

	V.dlogLik[2,3] <- V.dlogLik[3,2] <- V.dlogLik.23 - trace(V.naive.i%*%dvdr)*trace(VWZ[1:nops,(nops+1):(2*nops)])
	V.dlogLik[2,4] <- V.dlogLik[4,2] <- V.dlogLik.24 - trace(V.naive.i%*%dvde)*trace(VWZ[1:nops,(nops+1):(2*nops)])

V.dlogLik/nsub
}

####	calculate the expectation of the derivative of the score equations (evaluated at the limites of the naive solutions)
####	(e.g. E[\ell''(\Theta^{\dagger})]) 
####	rows correspond to the four score equations, 
####	columns correspond to \beta_0, \beta_1, \sigma^2_R, \sigma^2_{\epsilon}
####	V.naive.i:	Inverse of $\Sigma_{red}
####  VV:		Variance of (Q - \beta^{\dagger}_0 - \beta^{\dagger}_1 S)
####  VVW:		Variance of S and (Q - \beta^{\dagger}_0 - \beta^{\dagger}_1 S)


calc.L2 <- function(beta.1,theta1,sig.t,sig.g,sig.e,sig.r,nops,nsub,RHO,rho.y,V.naive.i,VV,VVW){

		###L2 is the matrix of second derivatives
		L2 		<- matrix(0,nrow=4,ncol=4)
		###SIG.W is the variance of S
		SIG.W   	<- sig.t*matrix(1,nrow=nops,ncol=nops)+(theta1-sig.t)*calcT2(RHO,rho.y,nops)+sig.g*diag(1,nrow=nops)
		vec.1   	<- matrix(1,nrow=nops,ncol=1)
		L2[1,1]   	<- -t(vec.1)%*%V.naive.i%*%vec.1 
		L2[2,2]   	<- -trace(V.naive.i%*%SIG.W) 

		dvdr	   <- matrix(1,nrow=nops,ncol=nops)
		dvde	   <- diag(1,nrow=nops)

		dvide <- dvidr <- matrix(nrow=nops,ncol=nops)
		for (i in 1:nops) for (j in 1:nops) {dvide[i,j] <- - (V.naive.i%*%dvde%*%V.naive.i)[i,j]
								 dvidr[i,j] <- - (V.naive.i%*%dvdr%*%V.naive.i)[i,j]}

		mmat.e    <- rbind( cbind(diag(1,nrow=nops),diag(0,nrow=nops)),cbind(diag(0,nrow=nops),dvide))	
		mmat.r    <- rbind( cbind(diag(1,nrow=nops),diag(0,nrow=nops)),cbind(diag(0,nrow=nops),dvidr))	
		VWZ.e     <- trace((t(mmat.e)%*%VVW%*%mmat.e)[1:nops,(nops+1):(2*nops)])
		VWZ.r     <- trace((t(mmat.r)%*%VVW%*%mmat.r)[1:nops,(nops+1):(2*nops)])
		L2[2,3]   <- VWZ.r
		L2[3,2]   <- 2*VWZ.r
		L2[4,2]   <- 2*VWZ.e
		L2[2,4]   <- VWZ.e

		dvidrvi.dr <- dvidr%*%dvdr%*%V.naive.i+ V.naive.i%*%dvdr%*%dvidr
		L2[3,3]    <- trace(dvidrvi.dr%*%VV) - trace(dvidr%*%dvdr)
		dvidevi.de <- dvide%*%dvde%*%V.naive.i+ V.naive.i%*%dvde%*%dvide
		L2[4,4]    <- trace(dvidevi.de%*%VV) - trace(dvide%*%dvde)
		dvidrvi.de <- dvide%*%dvdr%*%V.naive.i+ V.naive.i%*%dvdr%*%dvide
		dvidevi.dr <- dvidr%*%dvde%*%V.naive.i+ V.naive.i%*%dvde%*%dvidr
		L2[3,4]    <- trace(dvidrvi.de%*%VV) - trace(dvide%*%dvdr)
		L2[4,3]    <- trace(dvidevi.dr%*%VV) - trace(dvidr%*%dvde)
	L2
	}


####	calculates the standard error of estimates of \beta_0, \beta_1, \sigma^2_R, \sigma^2_{\epsilon} 
####	when the score equations for the reduced model are fit


calcSE.red <- function(beta.1,theta1,sig.t,sig.g,sig.e,sig.r,nops,nsub,RHO,rho.y,wantVar=0){

		###Calculate VV, the variance of (Q - \beta^{\dagger}_0 - \beta^{\dagger}_1 S)
		T2      <- calcT2(RHO,rho.y,nops)
		bias	  <- calcBsSE( theta1=theta1,sig.t=sig.t,sig.g=sig.g,sig.e=sig.e,sig.r=sig.r,nops=nops,nsub=nsub,onlyBias=1)
		Vres    <- diag(sig.e,nrow=nops)+matrix(sig.r,nrow=nops,ncol=nops)
		Vbdif   <- ((theta1-sig.t)*T2+matrix(sig.t,nrow=nops,ncol=nops))*bias[1,1]^2
		Vbdif2  <- (beta.1+bias[1,1])^2*diag(sig.g,nrow=nops)
		VV	  <- Vres+Vbdif+Vbdif2
		###calcVV.emp(beta.1,theta1,sig.t,sig.g,sig.e,sig.r,nops,nsub,RHO,rho.y)

		###Estimates the variance that is presumed under the naive model (\Sigma_{red})
		sig.e.naive 	<- sig.e+bias[2,1]
		sig.r.naive 	<- sig.r+bias[3,1]
		beta.1.naive	<- beta.1+bias[1,1]
		beta.0       	<- c(0,beta.1.naive,sig.r.naive,sig.e.naive)
		V.naive  		<- diag(sig.e.naive,nrow=nops)+matrix(sig.r.naive,nrow=nops,ncol=nops)
		V.naive.i 		<- solve(V.naive)

		###The true variance of S
		VW	  <- sig.t+(theta1-sig.t)*T2+sig.g*diag(1,nrow=nops)

		###Covariance between 	S and (Q - \beta^{\dagger}_0 - \beta^{\dagger}_1 S)
		CVW     <- -(sig.t+(theta1-sig.t)*T2)*bias[1,1]-diag(1,nrow=nops)*(beta.1+bias[1,1])*sig.g
		
		###Variance of S and  (Q - \beta^{\dagger}_0 - \beta^{\dagger}_1 S)
		VVW     <- rbind( cbind(VW,CVW),cbind(CVW,VV))
	
		###Variance of S and  \Sigma^{-1}_{red}(Q - \beta^{\dagger}_0 - \beta^{\dagger}_1 S)
		mmat    <- rbind( cbind(diag(1,nrow=nops),diag(0,nrow=nops)),cbind(diag(0,nrow=nops),V.naive.i))	
		VWZ     <- t(mmat)%*%VVW%*%mmat

		###var(\ell' (\Theta^{\dagger}))
		V.dlogLik <- calc.Vscore(nops,VWZ,V.naive.i,nsub)
		###calcVWZandV.dlogLik.emp(beta.1,theta1,sig.t,sig.g,sig.e,sig.r,nops,nsub,RHO,rho.y)

		####E[\ell''(\Theta^{\dagger})]
		L2        <- calc.L2(beta.1,theta1,sig.t,sig.g,sig.e,sig.r,nops,nsub,RHO,rho.y,V.naive.i,VV,VVW)
		L2.i	    <- solve(L2)
		SE	    <- sqrt(diag(signif(L2.i%*%V.dlogLik%*%t(L2.i),3)))
		####calcL2SE.emp(beta.1,theta1,sig.t,sig.g,sig.e,sig.r,nops,nsub,RHO,rho.y,beta.0,V.naive.i,L2,nsim=10,wantDeriv=1,wantlmat =1)

		if (wantVar==0) out <- SE
		if (wantVar==1) out <- L2.i%*%V.dlogLik%*%t(L2.i)

	out
	}



###	calculates both the bias and standard error of estimates produced when the reduced model is assumed

calcBsSE <- function(theta1=1,sig.t=0,sig.g=1,beta.1=1,sig.r=0,sig.e=1,d2=c(1:2),nops=2,nsub=100,onlyBias=0,RHO="I",rho.y=0){

	###obtains the naive estimates when presuming the reduced model
	###switching notation to Wang, Lin, Gutierrez, and Carroll (eq 17 and 18)

	theta.w <- sig.r
	beta.x  <- beta.1
	sig.xu  <- sig.t
	sig.x	  <- theta1-sig.t
	sig.u   <- sig.g
	sig	  <- sig.e
	n	  <- nops

	model <- function(xxx){
	###xxx <- beta.1.naive theta.naive and sig2.naive (don't need to worry about b0)

	V       <- diag(xxx[3],nrow=nops)+matrix(xxx[2],nrow=nops,ncol=nops)
	V.i     <- solve(V)
	SIG.X   <- sig.xu*matrix(1,nrow=nops,ncol=nops)+sig.x*calcT2(RHO,rho.y,nops)

	dV.dsig <- diag(1,nrow=nops)
	dV.dthe <- matrix(1,nrow=nops,ncol=nops)
	SIG.dV1 <- V.i%*%dV.dsig%*%V.i
	SIG.dV2 <- V.i%*%dV.dthe%*%V.i

	
	scoreEq <- function(V.i,SIG.dV,dV.d,SIG.X,beta.1,xxx,nops){
		(beta.1-xxx[1])*tr(SIG.dV%*%SIG.X)*(beta.1-xxx[1])+
		tr(SIG.dV%*%matrix(1,nrow=nops,ncol=nops)*theta.w)+
		tr(SIG.dV*sig)+
		xxx[1]*tr(SIG.dV*sig.u)*xxx[1]-
		tr(V.i%*%dV.d)
	}

	F1      <- tr(V.i%*%SIG.X)*(beta.1-xxx[1]) - tr(V.i*sig.u)*xxx[1]
	F2      <- scoreEq(V.i,SIG.dV1,dV.dsig,SIG.X,beta.1,xxx,nops)
	F3      <- scoreEq(V.i,SIG.dV2,dV.dthe,SIG.X,beta.1,xxx,nops)
	c(F1=F1,F2=F2,F3=F3)}

	naive.est <- multiroot(f=model,start=c(beta.1,sig.r,sig.e))$root


	###obtains the resulting biases

	bias.beta 	<- naive.est[1] - beta.1
	bias.sig.r 	<- naive.est[2] - sig.r
	bias.sig.e 	<- naive.est[3] - sig.e


	###obtains the standard errors of the naive estimates when presuming the reduced model
	se.red	<- c(0,0,0,0)
	if (onlyBias != 1) se.red  <- calcSE.red(beta.1,theta1,sig.t,sig.g,sig.e,sig.r,nops,nsub,RHO,rho.y)

	sum.res <- rbind( c(bias.beta, se.red[2]),
		     		c(bias.sig.e,se.red[4]),
		     		c(bias.sig.r,se.red[3]))
row.names(sum.res) <- c("beta.1","sig.e","sig.r")
sum.res
}


calcI.red  <- function(theta1=1,sig.t=0,sig.g=0.1,nops=2,nsub=10){
		derivs		<- list()	
		ld			<- 2
		V 			<- diag(theta1-sig.t+sig.g,nops)+matrix(sig.t,nrow=nops,ncol=nops)
		V.i			<- solve(V)
		derivs[[1]]  	<- matrix(1,nrow=nops,ncol=nops)
		derivs[[2]] 	<- diag(1,nrow=nops)
		###I.mat is the information matrix
		I.mat <- matrix(nrow=ld,ncol=ld)
		###each term is calculated from the equations on the wikipedia page
		for (i in 1:ld) for (j in 1:ld) I.mat[i,j] <- 0.5*sum(diag(V.i%*%derivs[[i]]%*%V.i%*%derivs[[j]]))
solve(I.mat)/nsub
}


###gets the MSE of the attenuation factor for both the full and reduced models
###set wantCor=0 to get the MSE of the correlation coefficient

calc.lambda.mse <- function( sig.t,sig.d,sig.u,sig.e,sig.r,beta.1,N,nsub,nsim=10,wantCor=0){
	
	calc.lambda2 <- function(beta.1,sig.e,sig.r,sig.t,sig.d){
			beta.1*(sig.t+sig.d)/( beta.1^2*(sig.t+sig.d) + sig.r + sig.e	)
	}
	calc.lambda <- function(tv) calc.lambda2(tv[1],tv[2],tv[3],tv[4],tv[5])

	calc.cor2 <- function(beta.1,sig.e,sig.r,sig.t,sig.d){
			beta.1*(sig.t+sig.d)/(sqrt( beta.1^2*(sig.t+sig.d) + sig.r + sig.e)*sqrt(sig.t+sig.d))
	}
	calc.cor <- function(tv) calc.cor2(tv[1],tv[2],tv[3],tv[4],tv[5])


	lambda.true <- calc.lambda2(beta.1,sig.e,sig.r,sig.t,sig.d)
	if (wantCor==1) lambda.true   <- calc.cor2(beta.1,sig.e,sig.r,sig.t,sig.d)
	if (wantCor==1) calc.lambda	<- calc.cor
	bias 		<- calcBsSE2( sig.t=sig.t,sig.d=sig.d,sig.u=sig.u,sig.e=sig.e,sig.r=sig.r,beta.1=beta.1,
						N=N,nsub=nsub,RHO=RHO,rho.y=0)
	V		<- calcSE.red(beta.1,(sig.t+sig.d),sig.t,sig.u,sig.e,sig.r,N,nsub,RHO="I",rho.y=0,wantVar=1)[c(2,4,3),c(2,4,3)]
	mu		<- bias[,1]+c(beta.1,sig.e,sig.r)

	V2		<- calcI.red(theta1=(sig.t+sig.d),sig.t=sig.t,sig.g=sig.u,nops=N,nsub=nsub)
	mu2		<- c(sig.t,sig.d+sig.u)

	X  		<- cbind(rmvnorm(nsim,mu,V),rmvnorm(nsim,mu2,V2))
	for (i in 2:5) X[,i] <- ifelse(X[,i]<0,0,X[,i])
	MSE.red	<- mean(  (apply(t(X),2,calc.lambda)-lambda.true)^2)

	V2		<- calcSE(sig.t+sig.d,sig.t,sig.u,beta.1,sig.r,sig.e,d=c(1,2,3,4,5,6),N,nsub,RHO="I",rho.y=0,mu=0,wantVar=1)[c(5,2,6,3,1),c(5,2,6,3,1)]
	mu2		<- c(beta.1,sig.e,sig.r,sig.t,sig.d)
	X2  		<- rmvnorm(nsim,mu2,V2)
	for (i in 2:5) X2[,i] <- ifelse(X2[,i]<0,0,X2[,i])
	MSE.full	<- mean(  (apply(t(X2),2,calc.lambda)-lambda.true)^2)

c(sqrt(MSE.full)/lambda.true,sqrt(MSE.red)/lambda.true)
}



####
####	Functions need to get emperical estimates of the standard errors and the bias
####

###generates the dataset
genData <- function(theta1=1,sig.t=0,sig.g=1,beta.1=1,sig.r=0,sig.e=1,nops=2,RHO="I",rho.y=0,nsub=100){
		V   <- calcV(theta1,sig.t,sig.g,beta.1,sig.r,sig.e,nops,RHO,rho.y)
		mu  <- rep(0,nrow(V))
		X   <- rmvnorm(nsub,mu,V)
X
}

###
### Want to Calculate the MLE when assuming the correct full model
### calculating the likelihood depents on RHO

###calculate the log-likelihood for the full model
logLik.I  <- function(theta1=1,sig.e=1,sig.g=1,beta.1=1,sig.r=0,tcor=0){
		V   <- calcV(theta1,tcor*theta1,sig.g,beta.1,sig.r,sig.e,nops,RHO,rho.y)
		mu  <- rep(0,nrow(V))
		ll <- sum(-log(dmvnorm(X,mu,V)))
ifelse(is.na(ll) | abs(ll) > 99999999,999999,ll)
}


###calculate the log-likelihood for the full model
logLik.AT  <- function(theta1=1,sig.e=1,sig.g=1,beta.1=1,sig.r=0,tcor=0,rho.y=0){
		V   <- calcV(theta1,tcor*theta1,sig.g,beta.1,sig.r,sig.e,nops,RHO,rho.y)
		mu  <- rep(0,nrow(V))
		ll <- sum(-log(dmvnorm(X,mu,V)))
ifelse(is.na(ll) | abs(ll) > 99999999,999999,ll)
}

###calculate the log-likelihood for the full model
logLik.A  <- function(theta1=1,sig.e=1,sig.g=1,beta.1=1,sig.r=0,tcor=0,rho.y=0){
		V   <- calcV(theta1,tcor*theta1,sig.g,beta.1,sig.r,sig.e,nops,RHO="A",rho.y)
		mu  <- rep(0,nrow(V))
		ll <- sum(-log(dmvnorm(X,mu,V)))
ifelse(is.na(ll) | abs(ll) > 99999999,999999,ll)
}

###calculate the log-likelihood for the full model
logLik.T  <- function(theta1=1,sig.e=1,sig.g=1,beta.1=1,sig.r=0,tcor=0,rho.y=0){
		V   <- calcV(theta1,tcor*theta1,sig.g,beta.1,sig.r,sig.e,nops,RHO="T",rho.y)
		mu  <- rep(0,nrow(V))
		ll <- sum(-log(dmvnorm(X,mu,V)))
ifelse(is.na(ll) | abs(ll) > 99999999,999999,ll)
}

###finds the MLE for the full model
###it is possible to presume certain parameters are known
###known parameters must be removed from the "start =" list and correctly specified as defaults in the logLik functions
###functions calculates the MLE estimates (and SE) for those variables listed in the "start = "
###		the list for "start =" must be in the same order as the input requested in logLik

calcMLE.I <- function(tn=6,nlow=5){
		tmle    <- suppressWarnings(mle(logLik.I,start =list(theta1=10,sig.e=1,sig.g=1,beta.1=1,sig.r=0.1,tcor=0.1) ,method="L-BFGS-B",
                                    lower=c(rep(0.01,nlow),   rep(-0.99,tn-nlow)),
                                    upper=c(rep(9999999,nlow),rep(0.99, tn-nlow))))
coef(summary(tmle))
}

calcMLE.AT <- function(tn=7,nlow=6){
		tmle    <- suppressWarnings(mle(logLik.AT,start =list(theta1=10,sig.e=1,sig.g=1,beta.1=1,sig.r=0.1,tcor=0.1,rho.y=0.1) ,method="L-BFGS-B",
                                    lower=c(rep(0.01,nlow),   rep(-0.99,tn-nlow)),
                                    upper=c(rep(9999999,nlow),rep(0.99, tn-nlow))))
coef(summary(tmle))
}

###
### (outdated) reduced model
###

#logLik2  <- function(sig.e2=1,sig.r=0){
#		Va            <- diag(sig.e2,nrow=nops)
#		V		  <- Va+sig.r
#		mu     <- rep(0,nrow(V))
#		ll     <- sum(-log(dmvnorm(X2,mu,V)))
#ifelse(is.na(ll) | abs(ll) > 99999999,999999,ll)
#}


#calcMLE2 <- function(){
#		g1	  <- cov(X2[,1],X2[,2])
#		g2      <- var(X2[,1])-g1
#		tmle    <- suppressWarnings(mle(logLik2,start =list(sig.e2=g2,sig.r=g1) ,method="L-BFGS-B",
#                                   lower=c(0.01,0.01),
#                                    upper=c(9999999,100)))
#coef(summary(tmle))
#}

###########################################################################################################
###
### When calculating the Variances of the parameters for the reduced models, it can be beneficial
### to compare the theoretical values w/ emperical values
###

	calcVV.emp <- function(beta.1,theta1,sig.t,sig.g,sig.e,sig.r,nops,nsub,RHO,rho.y){
		X       <- genData(theta1=theta1,sig.t=sig.t,sig.g=sig.g,beta.1=beta.1,sig.e=sig.e,sig.r=sig.r,nops=nops,RHO=RHO,rho.y=rho.y)
		Q  	    <- c(X[,1:nops])
		S  	    <- c(X[,(nops+1):(2*nops)])
		ID 	    <- rep(c(1:nsub),nops)
		tm 	    <- lmer(Q~S+(1|ID))
		time	    <- rep(c(1:nops),each=nsub)
		resid.tm <- Q-fixef(tm)[1]-(fixef(tm)[2])*S
		mtime <- function(t1,t2) mean(as.numeric(resid.tm)[time==t1]*as.numeric(resid.tm)[time==t2])
		VV.emp <- matrix(nrow=nops,ncol=nops)
		for (i in 1:nops) for (j in 1:nops) VV.emp[i,j] <- mtime(i,j)
	VV.emp
	}



	###Check that the estimated variance is correct
	calcVWZandV.dlogLik.emp <- function(beta.1,theta1,sig.t,sig.g,sig.e,sig.r,nops,nsub,RHO,rho.y){
		dvdr	   <- matrix(1,nrow=nops,ncol=nops)
		dvde	   <- diag(1,nrow=nops)
		l.mat  <- matrix(nrow=4,ncol=nsub)
		WZ.mat <- matrix(nrow=6,ncol=nsub)
		X       <- genData(theta1=theta1,sig.t=sig.t,sig.g=sig.g,beta.1=beta.1,sig.e=sig.e,sig.r=sig.r,nops=nops,RHO=RHO,rho.y=rho.y)		
		Q  	    <- c(X[,1:nops])
		S  	    <- c(X[,(nops+1):(2*nops)])
		ID 	    <- rep(c(1:nsub),nops)
		for (sub in 1:nsub){
			w0 <- matrix(1,nrow=nops,ncol=1)
			w1 <- matrix(S[ID==sub],ncol=1)
			Z  <- V.naive.i%*%matrix(Q[ID==sub]-w1*beta.1.naive,ncol=1)
			l.mat[,sub] <- matrix(c(t(w0)%*%Z,t(w1)%*%Z,t(Z)%*%dvdr%*%Z-trace(V.naive.i%*%dvdr),t(Z)%*%dvde%*%Z-trace(V.naive.i%*%dvde)),ncol=1)
			WZ.mat[,sub]      <- c(w1,Z)
		}

		l.mat.cent  <- l.mat-matrix(rowMeans(l.mat),byrow=F,nrow=nrow(l.mat),ncol=ncol(l.mat))
		WZ.mat.cent <- WZ.mat-matrix(rowMeans(WZ.mat),byrow=F,nrow=nrow(WZ.mat),ncol=ncol(WZ.mat))

	VWZ.emp	  <- round(WZ.mat.cent%*%t(WZ.mat.cent)/nsub,6)
	V.dlogLik.emp <- round(l.mat.cent%*%t(l.mat.cent)/nsub^2,6)
	list("VWZ.emp"= VWZ.emp,"V.dlogLik.emp"=  V.dlogLik.emp)
	}


	calcL2SE.emp   <- function(beta.1,theta1,sig.t,sig.g,sig.e,sig.r,nops,nsub,RHO,rho.y,beta.0,V.naive.i,L2,nsim=100,wantDeriv=0,wantlmat =0){
		dm      <- dm2 <- diag(0.01,nrow=4)
		sim.mat   <- matrix(nrow=4,ncol=nsim)
		lmat	    <- matrix(nrow=4,ncol=nsim)
		lmat.deriv   <- list()
		beta.0.naive <- beta.0[1]       
		beta.1.naive <- beta.0[2]       
		sig.r.naive  <- beta.0[3]       
		sig.e.naive  <- beta.0[4]       
		dvdr	   	 <- matrix(1,nrow=nops,ncol=nops)
		dvde	   	 <- diag(1,nrow=nops)

		for (ii in 1:4) lmat.deriv[[ii]] <- matrix(nrow=4,ncol=nsim)
		for (theSim in 1:nsim){
			X         <- genData(theta1=theta1,sig.t=sig.t,sig.g=sig.g,beta.1=beta.1,sig.e=sig.e,sig.r=sig.r,nops=nops,RHO=RHO,rho.y=rho.y)
			Q  	    <- c(X[,1:nops])
 			S  	    <- c(X[,(nops+1):(2*nops)])
			ID 	    <- rep(c(1:nsub),nops)
			tm 	    <- lmer(Q~S+(1|ID))
			s2e       <- (attr(VarCorr(tm), "sc"))^2
			s2r       <- VarCorr(tm)$ID[1]
			sim.mat[,theSim] <- c(fixef(tm),s2r,s2e)


			if (wantlmat == 1 | wantDeriv==1){
			l.mat2  <- l.mat3 <- matrix(nrow=4,ncol=nsub)
			for (sub in 1:nsub){
				w0 <- matrix(1,nrow=nops,ncol=1)
				w1 <- matrix(S[ID==sub],ncol=1)
				Z  <- V.naive.i%*%matrix(Q[ID==sub]-w1*beta.1.naive,ncol=1)
				l.mat2[,sub] <- matrix(c(t(w0)%*%Z,t(w1)%*%Z,t(Z)%*%dvdr%*%Z-trace(V.naive.i%*%dvdr),t(Z)%*%dvde%*%Z-trace(V.naive.i%*%dvde)),ncol=1)
			}
			lmat[,theSim] 		   <- rowMeans(l.mat2)
			}

			if (wantDeriv == 1){
			for (ii in 1:4) {
					V.naive2  <- diag(sig.e.naive+dm[ii,4],nrow=nops)+matrix(sig.r.naive+dm[ii,3],nrow=nops,ncol=nops)
					V.naive.i2 <- solve(V.naive2)
					for (sub in 1:nsub) {
						w0 <- matrix(1,nrow=nops,ncol=1)
						w1 <- matrix(S[ID==sub],ncol=1)
						Z3 <- V.naive.i2%*%matrix(Q[ID==sub]-dm[ii,1]-w1*(beta.1.naive+dm[ii,2]),ncol=1)
						l.mat3[,sub] <- matrix(c(t(w0)%*%Z3,t(w1)%*%Z3,t(Z3)%*%dvdr%*%Z3-trace(V.naive.i2%*%dvdr),t(Z3)%*%dvde%*%Z3-trace(V.naive.i2%*%dvde)),ncol=1)	
					}
			lmat.deriv[[ii]][,theSim]  <- (rowMeans(l.mat3)-rowMeans(l.mat2))/0.01
			}
			}

			
			
		}
			if (wantDeriv==1) for (ii in 1:4) dm2[,ii] <- rowMeans(lmat.deriv[[ii]])
			if (wantlmat==1)  {dif   <- sim.mat-matrix(beta.0,byrow=F,nrow=4,ncol=ncol(sim.mat))
						 lmat.est  <- -L2%*%dif}
			sim.mat.cent  <- sim.mat-matrix(rowMeans(sim.mat),byrow=F,nrow=nrow(sim.mat),ncol=ncol(sim.mat))
			beta.0.est    <- rowMeans(sim.mat)
	list("SE"=sqrt(diag(sim.mat.cent%*%t(sim.mat.cent)/nsim)), "L2"=dm2)
	}


###########################################################################################################
###############################################################################################################

###Unclear whether this is still a fair check of my functions

fakeComp <- function(){


	##rm(list=ls())
	library(mvtnorm)
	library(stats4)
	library(lme4)
	library(rootSolve)

	sig.r   <- 1
	beta.1  <- 1
	theta1  <- 7
	sig.t   <- 5
	sig.g   <- 1
	sig.e   <- 1
	nops    <- 3
	nsub    <- 10000
	rho.y   <- 0.9
	RHO	  <- "A"
	if (RHO == "I") 	{calcMLE <- calcMLE.I
				 logLik <- logLik.I
				 d = c(1:6)}
	if (RHO != "I") 	{calcMLE <- calcMLE.AT
				 logLik <- logLik.AT
				 d = c(1:7)}

	###all are unbiased so need to worry
	X       <- genData(theta1=theta1,sig.t=sig.t,sig.g=sig.g,beta.1=beta.1,sig.e=sig.e,sig.r=sig.r,nops=nops,RHO=RHO,rho.y=rho.y)
	SE	  <- calcSE(   theta1=theta1,sig.t=sig.t,sig.g=sig.g,beta.1=beta.1,sig.e=sig.e,sig.r=sig.r,d=d,nops=nops,nsub=nsub,RHO=RHO,rho.y=rho.y)
	SE
	calcMLE(tn=max(d))



	X       <- genData(theta1=theta1,sig.t=sig.t,sig.g=sig.g,beta.1=beta.1,sig.e=sig.e,sig.r=sig.r,nops=nops,RHO=RHO,rho.y=rho.y)
	Q  	    <- c(X[,1:nops])
 	S  	    <- c(X[,(nops+1):(2*nops)])
	ID 	    <- rep(c(1:nsub),nops)
	tm 	    <- lmer(Q~S+(1|ID))
	beta.ests   <- c(fixef(tm)[2]-beta.1,sqrt(diag(vcov(tm)))[2])
####	X2 	    <- X[,1:nops]-(fixef(tm)[1] + fixef(tm)[2]*X[,(nops+1):(2*nops)])
	sig.e.ests<- calcMLE()
	empRes    <- rbind(beta.ests,
                   c(sig.e.ests[2,1]-sig.e,sig.e.ests[2,2]), 
                   c(sig.e.ests[5,1]-sig.r,sig.e.ests[5,2]))
	rownames(empRes) <- c("beta.1","sig.e","sig.r")
####	extra.bias  <- c((attr(VarCorr(tm), "sc"))^2-sig.e,VarCorr(tm)$ID[1]-sig.r)
	bias	  <- calcBsSE( theta1=theta1,sig.t=sig.t,sig.g=sig.g,sig.e=sig.e,sig.r=sig.r,nops=nops,nsub=nsub)
	bias							
	empRes

}

###########################################################################################################
###############################################################################################################



SEandBIAS  <- function(sig.t=0,sig.d=1,sig.u=1,beta.1=1,sig.r=0,sig.e=1,N=2,nsub=100,RHO="I",rho.y=0){
outMat <- 
cbind(calcSE2(sig.t,sig.d,sig.u,beta.1,sig.r,sig.e,N,nsub,RHO,rho.y)[1:3],
      calcBsSE2(sig.t,sig.d,sig.u,beta.1,sig.r,sig.e,d2=c(1:6),N,nsub,0,RHO,rho.y))
colnames(outMat) <- c("SE.full","Bias.red","SE.red")
outMat
}

calcSE2  <- function(sig.t=0,sig.d=1,sig.u=1,beta.1=1,sig.r=0,sig.e=1,N=2,nsub=100,RHO="I",rho.y=0,mu=0){
			if (RHO=="I") d=c(1,2,3,4,5,6)
			if (RHO!="I") d=c(1,2,3,4,5,6,7)
			out <- calcSE(theta1=(sig.t+sig.d),sig.t=sig.t,sig.g=sig.u,beta.1=beta.1,sig.r=sig.r,sig.e=sig.e,d=d,nops=N,nsub=nsub,RHO=RHO,rho.y=rho.y,mu=mu)
			row.names(out)[c(1,4)] <- c("sig.d","sig.u")
			new.ord <- c(5,2,6,c(1:length(out))[-c(5,2,6)])
			out2 <- matrix(out[new.ord,],ncol=1)
			row.names(out2) <- row.names(out)[new.ord]
out2
}

calcBsSE2 <- function(sig.t=0,sig.d=1,sig.u=1,beta.1=1,sig.r=0,sig.e=1,d2=c(1:2),N=2,nsub=100,onlyBias=0,RHO="I",rho.y=0){
	calcBsSE(theta1=(sig.t+sig.d),sig.t=sig.t,sig.g=sig.u,beta.1=beta.1,sig.r=sig.r,sig.e=sig.e,d2=d2,nops=N,nsub=nsub,RHO=RHO,rho.y=rho.y)
}



