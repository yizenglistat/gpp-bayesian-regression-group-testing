#---------------------------------------------------------------------------------------------------------
# R code for calculating the estimator and ROT and PI bandwidths, for the paper by Delaigle and Meister
# "Nonparametric regression analysis for group testing data", 2011.
# HERE ALL GROUPS ARE OF EQUAL SIZE
# Important note: All bandwidth selectors are valid only for standard normal kenel but can be easily adapted to another 
# kernel. The local poly estimators are also all calculated with a standard normal kernel (dnorm function) but can be  
# easily adapted to another kernel.
#---------------------------------------------------------------------------------------------------------




#This routine is needed to calculate the ML estimator of q
fobject=function(qq)
{

	sumj=rep(0,J)
	for(j in 1:J)
		sumj[j]=sum(qq^(0:(nj[j]-1)))

	som=sum(nj/sum(sumj)*(Zstar-qq^nj))
	som
}


#This routine calculates the local linear estimator  of g at a single point x
fnorm=function(x,X,ZZ,h)
{
	n=length(X)
	a0=dnorm((x-X)/h)
	a1=a0*(X-x)/h
	a2=a1*(X-x)/h
	
	S0=sum(a0)/(n*h)
	S1=sum(a1)/(n*h)
	S2=sum(a2)/(n*h)
	T0=sum(ZZ*a0)/(n*h)
	T1=sum(ZZ*a1)/(n*h)
	
	(T0*S2-T1*S1)/(S0*S2-S1^2)
}



#This routine calculates the local linear estimator of g at a vector x
locpolvec=function(x,X,ZZ,h)
{

	sapply(x,fnorm,X,ZZ,h)

}



#This routine calculates, at a single point x, a local poly estimator of order ell of the nuth derivatiave g^{(nu)}
fnormder=function(x,X,ZZ,h,nu,ell)
{
	n=length(X)
	a0=dnorm((X-x)/h)
	S0=sum(a0)/(n*h)
	T0=sum(ZZ*a0)/(n*h)

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


#This routine calculates the weighted PI bandwidth for a standard normal kernel K, with weight {q_0.1,q_0.9}

hAMISEPI2=function(Xbar,Zstar)
{
	J=length(Xbar)
	#Second moment of a standard normal kernel
	muK2=1
	#\int K^2 for a standard normal kernel
	RK=1/(2*sqrt(pi))


	qantX=range(Xbar)
	x=seq(qantX[1],qantX[2],(qantX[2]-qantX[1])/100)
	
	qantX=quantile(Xbar,c(0.1,0.9))
	wofXi=dunif(Xbar,qantX[1],qantX[2])*(qantX[2]-qantX[1])
	
	#------------------------------------------------------------------------------
	#Calculation of the bandwidth h2 for nonparam estimation of the bias  term b. 
	#------------------------------------------------------------------------------
	
	#Global polynomial estimation of th24spline coming from the main program
	theta24=sum(th24spline*wofXi)/J
	
	#Constants related to the kernel K
	C2K=(3/(8*sqrt(pi)))^(1/7)
	if(theta24>0)
		C2K=(15/(16*sqrt(pi)))^(1/7)


	#-----------------------------------
	#Estimation of the variance term v
	#-----------------------------------
	
	
	#Find the index of the first component of each group	
	# cumsom=cumsum(nj)+1
	# cumsom=cumsom[-J]
	# cumsom=c(1,cumsom)
	# multi=qq^(-nj)*mean(ZZ)
	intvar=0
	
	# for(i in 1:max(nj))
	# {	
	
	# 	indices=cumsom+i-1
	# 	indord=order(X[indices])
	# 	Xord=X[indices[indord]]
	# 	Tord=multi*ZZ[indices[indord]]
		
	# 	intvar=intvar+sum(Tord[1:(J-1)]*(1-Tord[2:J])*(Xord[2:J]-Xord[1:(J-1)]))
	# }
	
	for(j in 1:(J-1)){
		intvar = intvar + Zstar[j]*(1-Zstar[j+1])*(Xbar[j+1]-Xbar[j])*wofXi[j]
	}



	# intvar=intvar/max(nj)


	#--------------
	#Bandwidth h2
	#--------------
		
	h2=C2K*(intvar/abs(theta24))^(1/7)*J^(-1/7)

	#-------------------------------------------------------
	#Calculation of the nonparam estimator of the bias term b. 
	#--------------------------------------------------------

	gofXipp2=0
	#ZZ=ZZ*qq^(-nj)*mean(ZZ)

	gofXipp=locpoldervec(Xbar,Xbar,Zstar,h2,2,3)
	gofXipp2=sum(gofXipp^2*wofXi)/J
	
	#-------------------------------------------------------
	#Calculation of the PI bandwidth. 
	#--------------------------------------------------------
	
	Num=RK*intvar
	Den=muK2^2*gofXipp2
	h=(Num/Den)^(1/5)*N^(-1/5)

	return(h)
}


#This routine calculates a global polynomial estimator or order ell of the nuth derivative of g, where nu can be 0, 2 or 4 
#and ell can be nu, nu+1 or nu+2. Can easily be extended to other derivatives nu and other values of ell
globpoly=function(x,Xbar,Zstar,nu,ell)
{
	stX=sqrt(var(Xbar))
	x=x/stX
	Xbar=Xbar/stX

	#find the index of the first component of each group
	cumsom=cumsum(nj)+1
	cumsom=cumsom[-J]
	cumsom=c(1,cumsom)
	hatbeta=rep(0,ell+1)
	
	#These are the T_j*'s
	#ZZ=ZZ*qq^(-Nj)*mean(ZZ)

	#Design matrix
	XD=0*Zstar+1
	for(j in 1:ell)
		XD=cbind(XD,Xbar^j)

	
	#Estimator of beta
	hatbeta=hatbeta+solve(t(XD)%*%XD)%*%t(XD)%*%Zstar
	
	rm(XD)


	#Estimator of g
	gofx=hatbeta[1]
	for(j in 1:ell)
		gofx=gofx+hatbeta[j+1]*x^j


	#Take the nuth derivative
	if((nu==2)&&(ell==4))
		gofxpp=2*hatbeta[3]+6*hatbeta[4]*x+12*hatbeta[5]*x^2
	

	if((nu==2)&&(ell==3))
		gofxpp=2*hatbeta[3]+6*hatbeta[4]*x


	if((nu==2)&&(ell==2))
		gofxpp=2*hatbeta[3]
		

	if((nu==4)&&(ell==6))
		gofxpp=24*hatbeta[5]+120*hatbeta[6]*x+720*hatbeta[7]*x^2
	

	if((nu==4)&&(ell==5))
		gofxpp=24*hatbeta[5]+120*hatbeta[6]*x

	if((nu==4)&&(ell==4))
		gofxpp=24*hatbeta[5]




	if(nu==0)
		gofxpp=gofx
	
	gofxpp


}

#---------------------------------------
#An example of how to use the codes
#---------------------------------------



#Function p(x) to estimate
pofx=function(x)
{
	num=sin(pi*x/2)+1.2
	sgnx=abs(x)/x
	sgnx[x==0]=0
	den=1+2*x^2*(sgnx+1)
	pp=num/(8*den)
	pp/2.5
}

#vector of values where to estimate p(x)
x_low = qunif(0.05,-3,3)
x_up = qunif(0.95,-3,3)

x = seq(x_low,x_up,0.01)


#---------------------------------------
#Generate the data
#---------------------------------------
	
#Total sample size

N=1000
n=N
	
	
#Generate the X_i's and the Y_i's
	
X=runif(n,-3,3)
X=X[order(X)]
pp=pofx(X)



Y=rep(0,n)
U=runif(n,0,1)
Y[U<pp]=Y[U<pp]+1


Xbar=filter(X,rep(1/max(nj),max(nj)), sides=1)[seq(max(nj),N,max(nj))]

#---------------------------------------------
# Group the data in groups of size G (here G=4)
#---------------------------------------------

#----------------------------------------------------------------------------------
# J, nj, Nj, Zstar and th24spline are GLOBAL variables that are used in the routines
#----------------------------------------------------------------------------------

G=20
nj=rep(G,N/G)
J=length(nj)

Ystar=rep(0,J)
for(j in 1:J)
{
	ind1=c(0,cumsum(nj[1:(J-1)]))+1
	ind2=cumsum(nj)

Ystar[j]=max(Y[ind1[j]:ind2[j]])
}



#Calculate Z*
Zstar=1-Ystar

#Estimate q. The resulting estimator qq is used as a GLOBAL variable in the routines
qq=uniroot(f=fobject,lower=0,upper=1)$root






#Calculate global poly estimator of theta24
ZZ=rep(Zstar,nj)
Nj=rep(nj,nj)

th24spline=globpoly(Xbar,Xbar,Zstar,2,3)*globpoly(Xbar,Xbar,Zstar,4,4)


#Calculate bandwidth (here unweighted ROT)
# hROT=hROT(X,ZZ)

#Calculate local poly estimator of g
# fit1=locpolvec(x,X,ZZ,hROT)

#Deduce estimator of p
# hatp1=1-qq*fit1/mean(ZZ)
# hatp1[hatp1<0]=0*hatp1[hatp1<0]
# hatp1[hatp1>1]=1+0*hatp1[hatp1>1]



#plot(x,pofx(x),type='l')
#lines(x,hatp1,col=2)



#Calculate bandwidth (here weighted PI)
hPI2=hAMISEPI2(Xbar,Zstar)

#Calculate local poly estimator of g
fit1=locpolvec(x,Xbar,Zstar,hPI2)
fit1[fit1<0] <- 0*fit1[fit1<0]
fit1[fit1>1] <- 1+0*fit1[fit1>1]
hatp1 = 1-fit1^(1/G)

#Deduce estimator of p
# hatp1=1-qq*fit1/mean(ZZ)
# hatp1[hatp1<0]=0*hatp1[hatp1<0]
# hatp1[hatp1>1]=1+0*hatp1[hatp1>1]

#lines(x,hatp1,col=3)

ise = mean((hatp1-pofx(x))^2)*10^4*(x_up-x_low)