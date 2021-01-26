# ---------------------------------------------- Preliminary ------------------------------------------------ #
# cat("\014")
rm(list=ls(all=TRUE))
graphics.off()
options(echo=TRUE,warn=-1)
# cluster <- TRUE
cluster <- FALSE
keys <- c('taskid','reps','nburn','nkeep','nthin','nstep',
	'isknown','cjs','Ns','models','random','homo','plot',
	'isIT','isMPT','isDT','isAT',
	'type','degree','penalty','niknots','a_taus','b_taus')
for(key in keys) assign(key, NA)
input_args <- commandArgs(trailingOnly = TRUE)

if(!cluster){
	#setwd("/Users/yizengli/Projects/mixed-vcm-only-one-random-simulation-2")
	# isknown = TRUE for IT and MPT
	# isknown = FALSE for DT and AT
	input_args <- c('taskid=1',
					'reps=1',
					'nburn=0',
					'nkeep=3',
					'nthin=2',
					'nstep=5',
					'cjs=c(1,5,10,20)',
					'Ns=c(1e3)',
					'models=c(1:8)',
					'random=TRUE',
					'homo=FALSE',
					'plot=FALSE',
					'penalty=2',
					'a_taus=c(7e-2)',
					'b_taus=c(9e-4)',
					'niknots=40',
					'isknown=FALSE',
					'isIT=FALSE',
					'isMPT=TRUE',
					'isDT=FALSE',
					'isAT=FALSE')
}

if(length(input_args)!=0){
	input_dict = unlist(strsplit(input_args,'='))
	input_keys = input_dict[seq(1,length(input_dict),2)]
	if(!all(input_keys %in% keys)) stop('[error] incorrect arguments names!')
	eval(parse(text=input_args))
}

if(identical(taskid, NA)) taskid <- 1
if(identical(reps, NA)) reps <- 1
if(identical(nburn, NA)) nburn <- 5000
if(identical(nkeep, NA)) nkeep <- 2000
if(identical(nthin, NA)) nthin <- 10
if(identical(nstep, NA)) nstep <- 20
if(identical(isknown, NA)) isknown <- FALSE
if(identical(cjs, NA)) cjs <- c(5,8)
if(identical(Ns, NA)) Ns <- 5000
if(identical(models, NA)) models=c(4,5,6,7,8,9,10,11)
if(identical(random, NA)) random <- TRUE
if(identical(homo, NA)) homo <- FALSE
if(identical(plot, NA)) plot <- FALSE
if(identical(isIT, NA)) isIT <- FALSE
if(identical(isMPT, NA)) isMPT <- FALSE
if(identical(isDT, NA)) isDT <- FALSE
if(identical(isAT, NA)) isAT <- FALSE
if(identical(type, NA)) type <- 'equal'
if(identical(degree, NA)) degree <- 3
if(identical(penalty, NA)) penalty <- 2
if(identical(niknots, NA)) niknots <- 40
if(identical(a_taus, NA)) a_taus <- 0.005
if(identical(b_taus, NA)) b_taus <- 0.005

source('./R/create_data.r')
source('./R/equal.r')
source('./R/testing_funs.r')
source('./R/save_data.r')
source('./R/bayes_plot.r')

dir.create('./log/', showWarnings = FALSE)
dir.create('./output/', showWarnings = FALSE)

for(rep in 1:reps){

	for(model in models){
		print(model)
		for(N in Ns){
			print(N)
			data_lst <- create_data(N=N, Se_true=c(1,1), Sp_true=c(1,1), model=model)

			for(cj in cjs){
				
				initial_path <- paste0('./output')
				dir.create(paste0(initial_path,'/model=',model), showWarnings = FALSE)
				dir.create(paste0(initial_path,'/model=',model,'/N=',N), showWarnings = FALSE)
				dir.create(paste0(initial_path,'/model=',model,'/N=',N,'/cj=',cj), showWarnings = FALSE)

				prevalence 	<- data_lst$prevalence
				N 			<- data_lst$N
				Se_true		<- data_lst$Se_true
				Sp_true		<- data_lst$Sp_true
				model_list	<- data_lst$model_list
				
				Y_true 		<- data_lst$Y_true
				X 			<- data_lst$X
				u 			<- data_lst$u
				u_pred 		<- data_lst$u_pred
				beta_true 	<- data_lst$beta_true
				g_pred_true <- data_lst$g_pred_true
				u_range 	<- data_lst$u_range

				Y_homo_true <- data_lst$Y_homo_true
				X_homo 		<- data_lst$X_homo
				u_homo 		<- data_lst$u_homo
				ubar_homo 	<- filter(u_homo,rep(1/cj,cj), sides=1)[seq(cj,N,cj)]

				if(random){
					MPT_res						<- try(masterpoolTesting(Y_true, Se_true[1], Sp_true[1], cj),silent=TRUE)
					Z_MPT						<- MPT_res$Z 								
					Y_MPT						<- MPT_res$Y
					Ystar 						<- Z_MPT[,1]							
					Zstar 						<- 1-Ystar
					X 							<- u 
					x 							<- u_pred
					G 							<- cj
					nj 							<- rep(cj,N/cj)
					J 							<- length(nj)
					qq 							<- uniroot(f=fobject,lower=0,upper=1)$root
					ZZ 							<- rep(Zstar,nj)
					Nj 							<- rep(nj,nj)
					th24spline 					<- globpoly(X,X,ZZ,2,3)*globpoly(X,X,ZZ,4,4)
					hPI2 						<- hAMISEPI2(X,ZZ)
					fit1 						<- locpolvec(x,X,ZZ,hPI2)
					hatp1 						<- 1-qq*fit1/mean(ZZ)
					hatp1[hatp1<0] 				<- 0*hatp1[hatp1<0]
					hatp1[hatp1>1] 				<- 1+0*hatp1[hatp1>1]

					g_pred 						<- hatp1

					savepath = paste0(paste0(initial_path,'/model=',model,'/N=',N,'/cj=',cj,'/random/'))
					dir.create(savepath, showWarnings = FALSE)
					save_data(u_pred, g_pred, g_pred_true, taskid=taskid, reps=reps, rep=rep, visual=FALSE, savepath=savepath)
				}
				if(FALSE){
						
					MPT_homo_res				<- try(masterpoolTesting(Y_homo_true, Se_true[1], Sp_true[1], cj),silent=TRUE)
					Z_homo_MPT					<- MPT_homo_res$Z 								
					Y_homo_MPT					<- MPT_homo_res$Y							
					Zstar_homo 					<- 1-Z_homo_MPT[,1]

					#Calculate bandwidth (here weighted PI)
					hPI2 						<- hAMISEPI2(u_homo, ubar_homo, Zstar_homo)

					#Calculate local poly estimator of g
					fit1 						<- locpolvec(u_pred, ubar_homo, Zstar_homo, hPI2)

					fit1[fit1<0]				<- 0*fit1[fit1<0]
					fit1[fit1>1] 				<- 1+0*fit1[fit1>1]
					#Deduce estimator of p
					g_pred 						<- 1 - fit1^(1/cj)


					savepath = paste0(paste0(initial_path,'/model=',model,'/N=',N,'/cj=',cj,'/homo/'))
					dir.create(savepath, showWarnings = FALSE)
					save_data(u_pred, g_pred, g_pred_true, taskid=taskid, reps=reps, rep=rep, visual=FALSE, savepath=savepath)


					#plot(u_pred, g_pred_true,type='l')
					#lines(u_pred, g_pred)

				}

			}

		}

	}
}