globpoly=function(x,Xbar,Zstar,nu,ell){

	#Xbar_std = sd(Xbar)
	#x = x/Xbar_std
	#Xbar = Xbar/Xbar_std

	hatbeta = rep(0, ell+1)

	#Design matrix
	XD=0*Zstar+1
	for(j in 1:ell)
		XD=cbind(XD,Xbar^j)
	
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