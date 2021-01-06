fnormder=function(x,Xbar,Zstar,h,nu,ell){

	S = matrix(0, ell+1, ell+1)
	for(i in 1:(ell+1)){
		for(j in 1:(ell+1)){
			k = i+j-2
			S[i,j] = sum((Xbar - x)^k * dnorm((Xbar-x)/h))
		}
	}

	T = matrix(0, ell+1, 1)
	for(k in 1:(ell+1)){
		T[k] = sum(Zstar * (Xbar - x)^k * dnorm((Xbar-x)/h))
	}

	sol = solve(S)%*%T
	
	out = factorial(nu)*h^(-nu)*sol[nu+1]
	
	return(sol)

}

locpoldervec=function(x,Xbar,Zstar,h,nu,ell){

	sapply(x,fnormder,Xbar,Zstar,h,nu,ell)

}