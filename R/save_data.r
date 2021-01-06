save_data<-function(u_pred, g_pred, g_true, taskid=character(0), reps=reps, rep=1, visual=FALSE, savepath=''){
		
	g_pred_summ <- data.frame(u_pred=u_pred, g_pred=g_pred, g_true=g_true, diffsq=(g_pred-g_true)^2)

	postfix <- reps*(taskid-1)+rep
	write.csv(g_pred_summ, file = paste0(savepath,'g_pred_summ',paste(rep('0',3-nchar(postfix)),collapse = ''),postfix,'.csv'),row.names=FALSE)

}

