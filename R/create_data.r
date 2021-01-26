create_data <- function(N=5000, Se_true, Sp_true, model=1){

	rtruncnorm <<- function(N, mean = 0, sd = 1, a = -Inf, b = Inf) {
	  if (a > b) stop('Error: Truncation range is empty');
	  U <- runif(N, pnorm(a, mean, sd), pnorm(b, mean, sd));
	  qnorm(U, mean, sd)
	}

	model_list 	<- list()
	# Define beta functions, age range [-3,3]
	if (model == 1){
		
		u <- runif(N, -3, 3)

		u_a <- quantile(u,0.05)
		u_b <- quantile(u,0.95)

		f0 <- function(x){
			A = sin(pi*x/2)+1.2
			B = 20+40*x^2*(sign(x)+1)
			return(A/B)
		}

		model_list[[1]]		<- list(beta_0=f0)

	}else if (model == 2){
		
		u <- rtruncnorm(N, 0, 1.5, -3, 3)

		u_a <- quantile(u,0.05)
		u_b <- quantile(u,0.95)

		f0 <- function(x){
			A = sin(pi*x/2)+1.2
			B = 20+40*x^2*(sign(x)+1)
			return(A/B)
		}

		model_list[[2]]		<- list(beta_0=f0)

	}else if (model == 3){
		
		u <- runif(N, -1, 4)

		u_a <- quantile(u,0.05)
		u_b <- quantile(u,0.95)

		f0 <- function(x){
			A = exp(-4+2*x)
			B = 8+8*exp(-4+2*x)
			return(A/B)
		}

		model_list[[3]]		<- list(beta_0=f0)

	}else if (model == 4){
		
		u <- rtruncnorm(N, 2, 1.5, -1, 4)

		u_a <- quantile(u,0.05)
		u_b <- quantile(u,0.95)

		f0 <- function(x){
			A = exp(-4+2*x)
			B = 8+8*exp(-4+2*x)
			return(A/B)
		}

		model_list[[4]]		<- list(beta_0=f0)

	}else if (model == 5){

		u <- runif(N, 0, 1)

		u_a <- qunif(0.05,0,1)
		u_b <- qunif(0.95,0,1)

		f0 <- function(x) return(x^2/8)

		model_list[[5]]		<- list(beta_0=f0)

	}else if (model == 6){

		u <- rtruncnorm(N, 0.5, 0.5, 0, 1)

		u_a <- quantile(u,0.05)
		u_b <- quantile(u,0.95)

		f0 <- function(x) return(x^2/8)

		model_list[[6]]		<- list(beta_0=f0)
		
	}else if (model == 7){

		u <- runif(N, -1, 1)

		u_a <- qunif(0.05,-1,1)
		u_b <- qunif(0.95,-1,1)

		f0 <- function(x) return(x^2/8)

		model_list[[7]]		<- list(beta_0=f0)

	}else if (model == 8){

		u <- rtruncnorm(N, 0, 0.75, -1, 1)

		u_a <- quantile(u,0.05)
		u_b <- quantile(u,0.95)

		f0 <- function(x) return(x^2/8)

		model_list[[8]]		<- list(beta_0=f0)
		
	}

	betaTrue <<- function(u, beta_list, center=FALSE){
		u_lower 		<- min(u); u_upper <- max(u)
		n 				<- length(u)
		nbeta 			<- length(beta_list) # total number of beta: beta0 beta1 beta2 beta3
		beta 			<- matrix(NA, n, nbeta)

		for (d in 1:nbeta){
			constant 	<- integrate(beta_list[[d]],u_lower,u_upper)$value/(u_upper-u_lower)
			beta[,d] 	<- logit_fun(beta_list[[d]](u) - constant*center)
		}
		return(beta)
	}
	
	logit_fun 	<<- function(x) log(x/(1-x))
	logit_inv 	<<- function(x) return(1/(1+exp(-x)))
	
	

	u_range 		<- u_b - u_a
	u_pred 			<- seq(u_a, u_b, by=0.01)
	beta_pred_true 	<- betaTrue(u_pred, model_list[[model]], center=FALSE) 

	X    			<- matrix(1,N,1)
	
	beta_true 		<- betaTrue(u, model_list[[model]], center=FALSE) # varying beta

	Y_true 			<- rep(-99, N)
	for(i in 1:N){
			eta 	<- X[i,] %*% beta_true[i,]
	        prob 	<- logit_inv(eta)
	        Y_true[i] <- rbinom(1, 1, prob)
	}

	prevalence 		<- mean(Y_true)
 	
 	g_pred_true 	<-logit_inv(beta_pred_true)

 	############################ homogenous
	
	u_homo			<- u[order(u)]
	Y_homo_true 	<- Y_true[order(u)]
	beta_homo_true 	<- as.matrix(beta_true[order(u),],N,nbeta)

	out_lst 		<- list(Se_true=Se_true, Sp_true=Sp_true,N=N,beta_true=beta_true,beta_pred_true=beta_pred_true,g_pred_true=g_pred_true,
							Y_true=Y_true, X=X, u=u, u_pred=u_pred, u_range=u_range, prevalence=prevalence,model_list=model_list,
							X_homo=X, u_homo=u_homo,Y_homo_true=Y_homo_true, beta_homo_true=beta_homo_true)
	return(out_lst)
}