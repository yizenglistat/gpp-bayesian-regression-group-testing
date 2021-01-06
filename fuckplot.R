bayes_plot <- function(model, N, cjs, plot=TRUE){

	source('./R/create_data.r')
	library(data.table)
	library(lattice)
	data_lst 	<- createData(N=N, Se_true=c(1,1), Sp_true=c(1,1), model=model)
	model_list 	<- data_lst$model_list
	u_range 	<- data_lst$u_range

	g_df <- c()
	ise_df 	<- c()

	for(homo in c(FALSE, TRUE)){
		for(cj in cjs){
			# configure the path
			if(!homo){
				path <- paste0('./output/a=',a_tau, 'b=',b_tau, '/model=',model,'/N=',N,'/cj=',cj,'/random/')
				method <- 'Random'
			}else{
				path <- paste0('./output/a=',a_tau, 'b=',b_tau, '/model=',model,'/N=',N,'/cj=',cj,'/homo/')
				method <- 'Homogeneous'
			}

			beta_pred_total <- c()
			for(file in list.files(path)){

				if(grepl('beta_new_summ',file)){
					beta_pred 		<- unlist(fread(paste0(path,file)))
					beta_pred_total <- cbind(beta_pred_total, beta_pred)
				}

			}

			alpha_mean 		<- apply(alpha_total, 1, mean)
			alpha_median 	<- apply(alpha_total, 1, median)
			alpha_1st 		<- apply(alpha_total, 1, quantile,probs=c(0.25))
			alpha_3rd 		<- apply(alpha_total, 1, quantile,probs=c(0.75))
			alpha_lower 	<- apply(alpha_total, 1, quantile,probs=c(0.025))
			alpha_upper 	<- apply(alpha_total, 1, quantile,probs=c(0.975))

			beta_pred_mean 		<- apply(beta_pred_total, 1, mean)
			beta_pred_median 	<- apply(beta_pred_total, 1, median)
			beta_pred_1st 		<- apply(beta_pred_total, 1, quantile,probs=c(0.25))
			beta_pred_3rd 		<- apply(beta_pred_total, 1, quantile,probs=c(0.75))
			beta_pred_lower 	<- apply(beta_pred_total, 1, quantile,probs=c(0.025))
			beta_pred_upper 	<- apply(beta_pred_total, 1, quantile,probs=c(0.975))

			g_pred_mean 		<- apply(g_pred_total, 1, mean)
			g_pred_median 		<- apply(g_pred_total, 1, median)
			g_pred_1st 			<- apply(g_pred_total, 1, quantile,probs=c(0.25))
			g_pred_3rd 			<- apply(g_pred_total, 1, quantile,probs=c(0.75))
			g_pred_lower 		<- apply(g_pred_total, 1, quantile,probs=c(0.025))
			g_pred_upper 		<- apply(g_pred_total, 1, quantile,probs=c(0.975))

			theta_mean 		<- apply(theta_total, 1, mean)
			theta_median 	<- apply(theta_total, 1, median)
			theta_1st 		<- apply(theta_total, 1, quantile,probs=c(0.25))
			theta_3rd 		<- apply(theta_total, 1, quantile,probs=c(0.75))
			theta_lower 	<- apply(theta_total, 1, quantile,probs=c(0.025))
			theta_upper 	<- apply(theta_total, 1, quantile,probs=c(0.975))

			u_pred_index 	<- grepl('u_pred', row.names(beta_pred_total))
			mean_index 		<- grepl('mean', row.names(beta_pred_total))
			median_index 	<- grepl('med_beta', row.names(beta_pred_total))
			lower_index 	<- grepl('lower_beta', row.names(beta_pred_total))
			upper_index 	<- grepl('upper_beta', row.names(beta_pred_total))

			npred 			<- sum(u_pred_index)
			nbeta 			<- sum(median_index)/npred
			label 			<- rep(paste0('beta',0:(nbeta-1)),each=npred)
			u_pred 			<- as.numeric(g_pred_median[u_pred_index])
			beta_pred_true 	<- betaTrue(u_pred, model_list[[model]], center=FALSE)

			g_diff 			<- g_pred_total[median_index,] - logit_inv(as.vector(beta_pred_true))
			ise 			<- apply((g_diff^2)[1:npred,],2,mean)*u_range # beta0
			ise_median 		<- round(median(ise)*10^4,3)
			ise_iqr 		<- round(IQR(ise)*10^4,3)

			ise_res 		<- data.frame(a_tau=a_tau,b_tau=b_tau,model = model, N=N, cj=cj, method=method, MED=ise_median, IQR=ise_iqr)

			g_res 			<- data.frame(u_pred=u_pred,
										true = logit_inv(as.vector(beta_pred_true)),
										median = g_pred_median[median_index],
										lower = g_pred_median[lower_index],
										upper = g_pred_median[upper_index],
										label = label,
										size = paste0('Pool Size=',cj),
										method = method,
										row.names=c())
			g_df <- rbind(g_df, g_res)
			ise_df 	<- rbind(ise_df, ise_res)
		}
	}

	print(ise_df)
	write.csv(ise_df, paste0('./output/model=',model,'N=',N,'ise_df.csv'),row.names=FALSE)

	if(plot){

	plot_config = list(line_col=c("#28B463",'#A93226'), band_col='gray', band_opacity=0.3)

	# confidence interval/prediction interval band
	mybands <- function(x, y, upper, lower, col, subscripts, ..., font, fontface){
		upper <- upper[subscripts]
		lower <- lower[subscripts]
		panel.polygon(c(x, rev(x)), c(upper, rev(lower)), col = band_col, alpha=band_opacity, border = FALSE,...)
	}
	# panel setting
	mypanel <- function(x, y, ...){
		panel.superpose(x, y, panel.groups=mybands, type='l',...)
		panel.xyplot(x, y, type='l', cex=0.6, lty=c(1,2,1,1),lwd=c(1,1,0.15,0.15), col=c(line_col,'black','black'),...)
	}
	# strip setting
	mystrip <- function(which.given, ..., factor.levels) {
		levs <- if (which.given == 1) factor.levels
		else   c('Random','Homogeneous')
		strip.default(which.given, ..., factor.levels = levs)
	}

	# configure xyplot setting
	line_col <- plot_config$line_col
	band_col <- plot_config$band_col
	band_opacity <- plot_config$band_opacity

	# generate figure 
	figure <- xyplot(true+median+lower+upper~u_pred |size*method,
			 		data=g_df,
			 		#data = replace(beta_res, "label", structure(beta_res$label, levels = )),
			 		auto.key = list(space='top',
			 					border=FALSE, 
			 					points=F, 
			 					lines=T, 
			 					#rectangles=T,
			 					columns = 3, 
			 					text=c('True','Median Estimate','95% Credible Interval Region')),
			   		lower=g_df$lower,upper=g_df$upper,
			   		panel = mypanel,
			   		par.settings =list(superpose.line = list(col=c(line_col,'#E5E7E9'), lwd = c(1,1,15), lty=c(1,2,1)),
			   							#superpose.polygon = list(col=c(NA,NA,'#D5D8DC'),border=c(0,0,0)),
			   							#strip.background=list(col=c('#9CC3D5FF','#3498DB'))),
			   							#strip.background=list(col=c('#3498DB','#9CC3D5FF'))),
			   							#strip.background=list(col=c('#F1F4FFFF','#A2A2A1FF'))),
			   							strip.background=list(col=c('#E5E7E9','#9CC3D5FF'))),
			 		xlab='u',
			 		ylab=bquote(bar(beta)(u)),
			 		#scales=list(y=list(relation="free")),
			 		#ylim=c(-0.5,1),
			 		strip=mystrip,
			 		layout=c(length(cjs),2))
	
	print(figure)
	
	}

}