hgplot <-
function(x){
	if(length(x[[1]])==2){
		stat=x[[1]][[1]]
		stat2=x[[1]][[2]]
	}else{
		stat=x[[1]]
	}
	conditional=stat$conditional
	mgn <- par(oma=c(0.1,0.1,0.1,0.1))

	legend1="# of the observed where localfdr <"
	legend2="# of the expected in the data."
	legend3="Total obs. count is"
			
	if(conditional==TRUE || length(x[[1]])==2){
		if(length(x[[1]])==2)
		{
			par(ask=TRUE)
			graphs=internalhgplot(stat)
			plotmat=graphs$plotmat
			label=graphs$label
			margin.y=graphs$margin.y
			e.total=graphs$e.total
			o.total=graphs$o.total
			
			legend1=paste(legend1,stat$plot$fdr.cutoff)
			legend3=paste(legend3,o.total); 
			barplot(plotmat,names.arg=label,beside=TRUE,col=c("gold","brown"), ylim=c(0,margin.y),ylab="Counts",
			xlab="COG",legend.text=c(legend1,legend2,legend3),cex.names=0.7,main="ALTERNATIVE : Over-Expression")
			
			graphs=internalhgplot(stat2)
			plotmat=graphs$plotmat
			label=graphs$label
			margin.y=graphs$margin.y
			e.total=graphs$e.total
			o.total=graphs$o.total

			legend3="Total obs. count is"
			legend3=paste(legend3,o.total); 
			barplot(plotmat,names.arg=label,beside=TRUE,col=c("gold","brown"), ylim=c(0,margin.y),ylab="Counts",
			xlab="COG",legend.text=c(legend1,legend2,legend3),cex.names=0.7,main="ALTERNATIVE : Under-Expression")
			par(ask=FALSE)
		}else{
			alternative=stat$alternative
			graphs=internalhgplot(stat)
			plotmat=graphs$plotmat
			label=graphs$label
			margin.y=graphs$margin.y
			e.total=graphs$e.total
			o.total=graphs$o.total
			switch(alternative,
				"over.expr"={
					legend1=paste(legend1,stat$plot$fdr.cutoff)
					legend3=paste(legend3,o.total); 
					barplot(plotmat,names.arg=label,beside=TRUE,col=c("gold","brown"), ylim=c(0,margin.y),ylab="Counts",
					xlab="COG",legend.text=c(legend1,legend2,legend3),cex.names=0.7,main="ALTERNATIVE : Over-Expression")
				},"under.expr"={
					legend1=paste(legend1,-1*stat$plot$fdr.cutoff)
					legend3=paste(legend3,o.total); 
					barplot(plotmat,names.arg=label,beside=TRUE,col=c("gold","brown"), ylim=c(0,margin.y),ylab="Counts",
					xlab="COG",legend.text=c(legend1,legend2,legend3),cex.names=0.7,main="ALTERNATIVE : Under-Expression")
				}
			)
			par(mgn)
		}
	}else{
		graphs=internalhgplot(stat)
		plotmat=graphs$plotmat
		label=graphs$label
		margin.y=graphs$margin.y
		e.total=graphs$e.total
		o.total=graphs$o.total
		
		legend1=paste(legend1,stat$plot$fdr.cutoff)
		legend3=paste(legend3,o.total);
		barplot(plotmat,names.arg=label,beside=TRUE,col=c("gold","brown"), ylim=c(0,margin.y),ylab="Counts",
		xlab="COG",legend.text=c(legend1,legend2,legend3),cex.names=0.7,main="Unconditional Analysis")
		par(mgn)
	}
}
