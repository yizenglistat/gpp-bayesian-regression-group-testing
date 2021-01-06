#---------------------------------------------------------------------------------------------------------
# R code for calculating the estimator and ROT and PI bandwidths, for the paper by Delaigle and Meister
# "Nonparametric regression analysis for group testing data", 2011.
# HERE ALL GROUPS CAN BE OF UNEQUAL SIZE
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



#This routine calculates a global polynomial estimator or order ell of the nuth derivative of g, where nu can be 0, 2 or 4 
#and ell can be nu, nu+1 or nu+2. Can easily be extended to other derivatives nu and other values of ell
globpoly=function(x,X,ZZ,nu,ell)
{
	stX=sqrt(var(X))
	x=x/stX
	X=X/stX

	hatbeta=rep(0,ell+1)
	#These are the T_j*'s
	ZZ=ZZ*qq^(-Nj)*mean(ZZ)

	#Design matrix
	XD=0*ZZ+1
	for(j in 1:ell)
		XD=cbind(XD,X^j)

	
	#Estimator of beta
	hatbeta=hatbeta+solve(t(XD)%*%XD)%*%t(XD)%*%ZZ
	
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




#This routine calculates the unweighted ROT

hROT=function(X,ZZ)
{
	N=length(X)
	muK2=1
	RK=1/(2*sqrt(pi))

	#Gobal polynomial estimation of the bias term
	gofXipp=globpoly(X,X,ZZ,2,3)

	#Estimation of the variance term
	
	cumsom=cumsum(nj)+1
	cumsom=cumsom[-J]
	cumsom=c(1,cumsom)
	mult=qq^(-nj)*mean(ZZ)
	intvar=0
	
	sumwi=0
	for(i in 1:max(nj))
	{	
		#keep only groups of size nj>=i
		condnj=(nj>=i)
		Ji=sum(condnj)
		wi=sqrt(Ji)
		sumwi=sumwi+wi
		cumsomi=cumsom[condnj]
		indices=cumsomi+i-1
		indord=order(X[indices])
		Xord=X[indices[indord]]
		Tord=mult[condnj]*ZZ[indices[indord]]
		intvar=intvar+wi*sum(Tord[1:(Ji-1)]*(1-Tord[2:Ji])*(Xord[2:Ji]-Xord[1:(Ji-1)]))
	}
	
	intvar=intvar/sumwi


	#Bandwidth
	Num=RK*intvar
	Den=muK2^2*sum(gofXipp^2)
	h=(Num/Den)^(1/5)
	h
}




#This routine calculates the ROT weighted by weight 1{q_0.1,q_0.9}

hROT2=function(X,ZZ)
{
	N=length(X)
	#Second moment of a standard normal kernel
	muK2=1
	RK=1/(2*sqrt(pi))

	#Gobal polynomial estimation of the bias term
	gofXipp=globpoly(X,X,ZZ,2,3)

	#Estimation of the variance term
	
	cumsom=cumsum(nj)+1
	cumsom=cumsom[-J]
	cumsom=c(1,cumsom)
	mult=qq^(-nj)*mean(ZZ)
	intvar=0
	
	sumwi=0
	for(i in 1:max(nj))
	{	
		#keep only groups of size nj>=i
		condnj=(nj>=i)
		Ji=sum(condnj)
		wi=sqrt(Ji)
		sumwi=sumwi+wi
		cumsomi=cumsom[condnj]
		indices=cumsomi+i-1
		indord=order(X[indices])
		Xord=X[indices[indord]]
		Tord=mult[condnj]*ZZ[indices[indord]]
		intvar=intvar+wi*sum(Tord[1:(Ji-1)]*(1-Tord[2:Ji])*(Xord[2:Ji]-Xord[1:(Ji-1)]))
	}
	
	intvar=intvar/sumwi

	qantX=quantile(X,c(0.1,0.9))
	wofXi=dunif(X,qantX[1],qantX[2])*(qantX[2]-qantX[1])



	#Bandwidth
	Num=RK*intvar
	Den=muK2^2*sum(gofXipp^2*wofXi)
	h=(Num/Den)^(1/5)
	h
}



#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------




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
x=seq(-3,3,6/100)




#---------------------------------------
#Generate the data
#---------------------------------------


#Total sample size
N=10000
n=N
	
	

#Generate the X_i's and the Y_i's
	
X=runif(n,-3,3)
pp=pofx(X)

Y=rep(0,n)
U=runif(n,0,1)
Y[U<pp]=Y[U<pp]+1


#----------------------------------------------------------------
# Group the data in groups of unequal sizes, here n_j ~ U[1,5]
#----------------------------------------------------------------

#----------------------------------------------------------------------------------
# J, nj, Nj, N, Zstar are GLOBAL variables that are used in the routines
#----------------------------------------------------------------------------------


nj=c()

while(sum(nj)<=N-5)
{
	nj=c(nj,round(runif(1,1,5)))
}

if(sum(nj)<N)
	nj=c(nj,N-sum(nj))

J=length(nj)


n=sum(nj)
N=sum(nj)
Nj=rep(nj,nj)


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



ZZ=rep(Zstar,nj)


#Calculate bandwidth (here unweighted ROT)
hROT=hROT(X,ZZ)

#Calculate local poly estimator of g
fit1=locpolvec(x,X,ZZ,hROT)


#Deduce estimator of p
hatp1=1-qq*fit1/mean(ZZ)
hatp1[hatp1<0]=0*hatp1[hatp1<0]
hatp1[hatp1>1]=1+0*hatp1[hatp1>1]

plot(x,pofx(x),type='l')
lines(x,hatp1,col=2)

