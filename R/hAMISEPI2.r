hAMISEPI2=function(X, Xbar, Zstar){

	N=length(X)
	
	nj=rep(cj,N/cj)
	
	J=length(Xbar)
	
	#Second moment of a standard normal kernel
	muK2=1
	#\int K^2 for a standard normal kernel
	RK=1/(2*sqrt(pi))

	qantX=quantile(X,c(0.1,0.9))
	
	wofXbarj=dunif(Xbar,qantX[1],qantX[2])*(qantX[2]-qantX[1])
	
	#------------------------------------------------------------------------------
	#Calculation of the bandwidth h2 for nonparam estimation of the bias  term b. 
	#------------------------------------------------------------------------------
	
	#Global polynomial estimation of th24spline coming from the main program
	# globpoly(X,X,ZZ,2,3)*globpoly(X,X,ZZ,4,4)
	th24spline=globpoly(Xbar,Xbar,Zstar,2,3)*globpoly(Xbar,Xbar,Zstar,4,4)
	theta24=sum(th24spline*wofXbarj)/J
	
	#Constants related to the kernel K
	C2K=(3/(8*sqrt(pi)))^(1/7)
	if(theta24>0)
		C2K=(15/(16*sqrt(pi)))^(1/7)


	#-----------------------------------
	#Estimation of the variance term v
	#-----------------------------------
	
	
	#Find the index of the first component of each group	
	intvar = 0
	for(j in 1:(J-1)){
		intvar = intvar + Zstar[j]*(1-Zstar[j+1])*(Xbar[j+1]-Xbar[j])*wofXbarj[j]
	}

	#--------------
	#Bandwidth h2
	#--------------
		
	h2=C2K*(intvar/abs(theta24))^(1/7)*J^(-1/7)

	#-------------------------------------------------------
	#Calculation of the nonparam estimator of the bias term b. 
	#--------------------------------------------------------
	
	gofXipp=locpoldervec(Xbar,Xbar,Zstar,h2,2,3)
	gofXipp2=sum(gofXipp^2*wofXbarj)/J
	
	#-------------------------------------------------------
	#Calculation of the PI bandwidth. 
	#--------------------------------------------------------
	
	
	Num=RK*intvar
	Den=muK2^2*gofXipp2
	h=(Num/Den)^(1/5)*N^(-1/5)
	h
}