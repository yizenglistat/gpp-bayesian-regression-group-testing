
#This routine calculates, at a single point x, a local poly estimator of order ell of the nuth derivatiave g^{(nu)}
fnormder=function(x,Xbar,Zstar,h,nu,ell)
{
	J = length(Xbar)
	a0=dnorm((X-x)/h)
	S0=sum(a0)/(J*h)
	T0=sum(Zstar*a0)/(J*h)

	if(ell>0)
	{	
		for(i in 1:(2*ell))
		{

			eval(parse(text=paste("a",i,"=a",i-1,"*(X-x)/h", sep = "")))
			eval(parse(text=paste("S",i,"=sum(a",i,")/(n*h)", sep = "")))
			eval(parse(text=paste("T",i,"=sum(ZZ*a",i,")/(n*h)", sep = "")))
		}
	}
	
	TT=c()
	SS=matrix(0,nrow=ell+1,ncol=ell+1)
	for(i in 0:ell)
		{
		eval(parse(text=paste("TT=c(TT,T",i,")", sep = "")))
		for(j in 0:ell)
			eval(parse(text=paste("SS[",i+1,",",j+1,"]=S",i+j, sep = "")))
		
		}

	sol=solve(SS)%*%TT
	sol=factorial(nu)*h^(-nu)*sol[nu+1]
	sol

}


#This routine calculates, at a whole vector x, a local poly estimator of order ell of the nuth derivatiave g^{(nu)}
locpoldervec=function(x,X,ZZ,h,nu,ell)
{

	sapply(x,fnormder,X,ZZ,h,nu,ell)

}



hAMISEPI2=function(X,Zstar)
{

	Xbar = filter(X,rep(1/max(nj),max(nj)), sides=1)[seq(4,N,max(nj))]
	N=length(X)
	J=length(Xbar)
	#Second moment of a standard normal kernel
	muK2=1
	#\int K^2 for a standard normal kernel
	RK=1/(2*sqrt(pi))

	qantX=range(X)
	x=seq(qantX[1],qantX[2],(qantX[2]-qantX[1])/100)
	
	qantX=quantile(X,c(0.1,0.9))
	
	wofXbarj=dunif(Xbar,qantX[1],qantX[2])*(qantX[2]-qantX[1])
	
	#------------------------------------------------------------------------------
	#Calculation of the bandwidth h2 for nonparam estimation of the bias  term b. 
	#------------------------------------------------------------------------------
	
	#Global polynomial estimation of th24spline coming from the main program
	# globpoly(X,X,ZZ,2,3)*globpoly(X,X,ZZ,4,4)
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
	gofXipp2=mean(gofXipp^2*wofXbarj)
	
	#-------------------------------------------------------
	#Calculation of the PI bandwidth. 
	#--------------------------------------------------------
	
	
	Num=RK*intvar
	Den=muK2^2*gofXipp2
	h=(Num/Den)^(1/5)*N^(-1/5)
	h
}