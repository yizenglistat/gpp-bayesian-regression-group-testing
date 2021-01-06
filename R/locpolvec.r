fnorm=function(x,Xbar,Zstar,h)
{
	J = length(Xbar)
	S0 = sum(dnorm((Xbar-x)/h))
	S1 = sum(dnorm((Xbar-x)/h)*(Xbar-x))
	S2 = sum(dnorm((Xbar-x)/h)*(Xbar-x)^2)

	T0 = sum(Zstar*dnorm((Xbar-x)/h))
	T1 = sum(Zstar*dnorm((Xbar-x)/h)*(Xbar-x))
	T2 = sum(Zstar*dnorm((Xbar-x)/h)*(Xbar-x)^2)
	
	out = (T0*S2-T1*S1)/(S0*S2-S1^2)
	return(out)
}

#This routine calculates the local linear estimator of g at a vector x
locpolvec=function(x,Xbar,Zstar,h)
{

	sapply(x,fnorm,Xbar,Zstar,h)

}

