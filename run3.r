source('./R/equal.r')

#Function p(x) to estimate
pofx=function(x){
	num=sin(pi*x/2)+1.2
	sgnx=abs(x)/x
	sgnx[x==0]=0
	den=1+2*x^2*(sgnx+1)
	pp=num/(8*den)
	pp/2.5
}

# pofx <- function(x){
# 	A = exp(-4+2*x)
# 	B = 8+8*exp(-4+2*x)
# 	return(A/B)
# }

reps=200
ise = c()
for(rep in 1:reps){
	N=1000#10000
	n=N

	#Generate the X_i's and the Y_i's
		
	X=runif(n,-3,3)
	a = qunif(0.05,-3,3)
	b = qunif(0.95,-3,3)
	x=seq(a,b,0.01)

	X=X[order(X)]

	pp=pofx(X)

	Y=rep(0,n)
	U=runif(n,0,1)
	Y[U<pp]=Y[U<pp]+1


	cj=G=20
	nj=rep(G,N/G)
	J=length(nj)


	Ystar=rep(0,J)
	for(j in 1:J){
		ind1=c(0,cumsum(nj[1:(J-1)]))+1
		ind2=cumsum(nj)
		Ystar[j]=max(Y[ind1[j]:ind2[j]])
	}

	Zstar=1-Ystar

	ZZ=rep(Zstar,nj)
	Nj=rep(nj,nj)

	Xbar = filter(X,rep(1/max(nj),max(nj)), sides=1)[seq(max(nj),N,max(nj))]

	#Calculate bandwidth (here weighted PI)
	hPI2=hAMISEPI2(X,Xbar,Zstar)

	#Calculate local poly estimator of g
	fit1=locpolvec(x,Xbar,Zstar,hPI2)
	fit1[fit1<0] <- 0*fit1[fit1<0]
	fit1[fit1>1] <- 1+0*fit1[fit1>1]
	#Deduce estimator of p
	hatp1 = 1 - fit1^(1/G)

	ise = c(ise, mean((hatp1-pofx(x))^2)*10^4*(b-a))
	plot(x,pofx(x),type='l',ylim=c(0,1))
	lines(x,hatp1,col=3)
}

print(cbind(median(ise),IQR(ise)))
#plot(x,pofx(x),type='l',ylim=c(0,1))
#lines(x,hatp1,col=3)

