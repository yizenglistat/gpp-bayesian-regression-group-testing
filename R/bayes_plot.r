bayes_plot <- function(models, N, cjs, plot=TRUE){
	source('./R/create_data.r')
	library(data.table)
	library(lattice)
	ise_df 	<- c()
	
	for(model in models){

		data_lst 	<- create_data(N=N, Se_true=c(1,1), Sp_true=c(1,1), model=model)
		model_list 	<- data_lst$model_list
		u_range 	<- data_lst$u_range
		
		for(homo in c(TRUE)){
			for(cj in cjs){
				# configure the path
				if(!homo){
					path <- paste0('./output/model=',model,'/N=',N,'/cj=',cj,'/random/')
					method <- 'Random'
				}else{
					path <- paste0('./output/model=',model,'/N=',N,'/cj=',cj,'/homo/')
					method <- 'Homogeneous'
				}

				diffsq_total 	<- c()
				for(file in list.files(path)){
					if(grepl('g_pred_summ',file)){
						diffsq 			<- fread(paste0(path,file))
						diffsq_total 	<- cbind(diffsq_total, diffsq[,'diffsq'])
					}
				}
				
				u_pred 			<- diffsq[,'u_pred']
				npred 			<- sum(u_pred)
				ise 			<- apply(diffsq_total,2,mean)*u_range
				ise_median 		<- round(median(ise)*10^4,3)
				ise_iqr 		<- round(IQR(ise)*10^4,3)

				ise_res 		<- data.frame(model = model, N=N, cj=cj, method=method, MED=ise_median, IQR=ise_iqr)

				ise_df 			<- rbind(ise_df, ise_res)
			}
		}
	}

	print(ise_df)
	write.csv(ise_df, paste0('./output/model=',model,'N=',N,'ise_df.csv'),row.names=FALSE)

}