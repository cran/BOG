gseaplot <-
function(x,cat=NULL){
	if(is.null(cat)) stop("BOG : Category needs to be specifiled.")
	if(length(x[[1]])==2){
		stat=x[[1]][[1]]
		stat2=x[[1]][[2]]
	}else{
		stat=x[[1]]
	}
	if(is.null(stat$gsea)) stop("BOG : gsea was not set to be TRUE in creating BOG object.")
	conditional=stat$conditional
	if(conditional==TRUE || length(x[[1]])==2){
		if(length(x[[1]])==2)
		{
			### Over Expressoin ###########
			par(mfrow=c(2,1))
			map=stat$plot$map
			delta=stat$plot$delta
			S=stat$plot$S_location
			
			index=values(map,keys=cat,USE.NAMES=FALSE)
			delta=stat$plot$delta[,index]
			
			pval=stat$gsea$raw[which(stat$gsea$COG==cat)]
			n.index=length(delta)
			max.y=max(delta)
			max.x=which(delta==max.y)
			
			horiz.y=rep(max.y,n.index)
			horiz.x=c(1:n.index)
			
			ver.y=seq(0,max.y+0.2,length.out=1000)
			ver.x=rep(max.x,1000)
			label=paste("Gene List Rank for over.expr alternative of COG",cat)
			label=paste(label,": p-value=")
			label=paste(label,round(pval,digits=5))
			plot(delta,ylim=c(0,1),xlab=label,ylab="Running Enrichment Score",
			type="l",col="brown")
			points(horiz.x,horiz.y,type="l",col="red",lty=3)
			points(ver.x,ver.y,type="l",col="red",lty=3)
			points(max.x,max.y,col="red",pch=18)
			axis(4,max.y,labels=round(max.y,digits=3),lty=3,las=0)
			axis(3,max.x,labels=max.x,lty=3,las=0)
			
			Spoints.x=which(S[,index]!=1)
			Spoints.ys=rep(0,length(Spoints.x))
			Spoints.ye=rep(0.1,length(Spoints.x))
			segments(Spoints.x,Spoints.ys,Spoints.x,Spoints.ye,col="green")
			
			Spoints.x=which(S[,index]==1)
			Spoints.ys=rep(0,length(Spoints.x))
			Spoints.ye=rep(0.1,length(Spoints.x))
			segments(Spoints.x,Spoints.ys,Spoints.x,Spoints.ye,col="red")
			
			### Under Expressoin ###########	
			map=stat2$plot$map
			delta=stat2$plot$delta
			S=stat2$plot$S_location
			
			index=values(map,keys=cat,USE.NAMES=FALSE)
			delta=stat2$plot$delta[,index]
			pval=stat2$gsea$raw[which(stat2$gsea$COG==cat)]
			n.index=length(delta)
			max.y=max(delta)
			max.x=which(delta==max.y)
			
			horiz.y=rep(max.y,n.index)
			horiz.x=c(1:n.index)
			
			ver.y=seq(0,max.y+0.2,length.out=1000)
			ver.x=rep(max.x,1000)
			
			label=paste("Gene List Rank for under.expr alternative of COG ",cat)
			label=paste(label,": p-value=")
			label=paste(label,round(pval,digits=5))
			plot(delta,ylim=c(0,1),xlab=label,ylab="Running Enrichment Score",type="l",col="brown")
			points(horiz.x,horiz.y,type="l",col="red",lty=3)
			points(ver.x,ver.y,type="l",col="red",lty=3)
			points(max.x,max.y,col="red",pch=18)
			axis(4,max.y,labels=round(max.y,digits=3),lty=3,las=0)
			axis(3,max.x,labels=max.x,lty=3,las=0)
			
			Spoints.x=which(S[,index]!=1)
			Spoints.ys=rep(0,length(Spoints.x))
			Spoints.ye=rep(0.1,length(Spoints.x))
			segments(Spoints.x,Spoints.ys,Spoints.x,Spoints.ye,col="green")
			
			Spoints.x=which(S[,index]==1)
			Spoints.ys=rep(0,length(Spoints.x))
			Spoints.ye=rep(0.1,length(Spoints.x))
			segments(Spoints.x,Spoints.ys,Spoints.x,Spoints.ye,col="red")
		}else{
			par(mfrow=c(1,1))
			map=stat$plot$map
			delta=stat$plot$delta
			S=stat$plot$S_location
			alternative=stat$alternative
			index=values(map,keys=cat,USE.NAMES=FALSE)
			delta=stat$plot$delta[,index]
			pval=stat$gsea$raw[which(stat$gsea$COG==cat)]
			n.index=length(delta)
			max.y=max(delta)
			max.x=which(delta==max.y)
			
			horiz.y=rep(max.y,n.index)
			horiz.x=c(1:n.index)
			
			ver.y=seq(0,max.y+0.2,length.out=1000)
			ver.x=rep(max.x,1000)
			
			switch(alternative,
				"over.expr"={
					label=paste("Gene List Rank for over.expr alternative of COG ",cat)
				},
				"under.expr"={	
					label=paste("Gene List Rank for under.expr alternative of COG ",cat)
				}
			)
			label=paste(label,": p-value=")
			label=paste(label,round(pval,digits=5))
			plot(delta,ylim=c(0,1),xlab=label,ylab="Running Enrichment Score",type="l",col="brown")
			points(horiz.x,horiz.y,type="l",col="red",lty=3)
			points(ver.x,ver.y,type="l",col="red",lty=3)
			points(max.x,max.y,col="red",pch=18)
			axis(4,max.y,labels=round(max.y,digits=3),lty=3,las=0)
			axis(3,max.x,labels=max.x,lty=3,las=0)
			
			Spoints.x=which(S[,index]!=1)
			Spoints.ys=rep(0,length(Spoints.x))
			Spoints.ye=rep(0.1,length(Spoints.x))
			segments(Spoints.x,Spoints.ys,Spoints.x,Spoints.ye,col="green")
			
			Spoints.x=which(S[,index]==1)
			Spoints.ys=rep(0,length(Spoints.x))
			Spoints.ye=rep(0.1,length(Spoints.x))
			segments(Spoints.x,Spoints.ys,Spoints.x,Spoints.ye,col="red")	
		}
	}else{
		par(mfrow=c(1,1))
		map=stat$plot$map
		delta=stat$plot$delta
		S=stat$plot$S_location
		alternative=stat$alternative
		index=values(map,keys=cat,USE.NAMES=FALSE)
		delta=stat$plot$delta[,index]
		pval=stat$gsea$adj[which(stat$gsea$COG==cat)]
		n.index=length(delta)
		max.y=max(delta)
		max.x=which(delta==max.y)
		
		horiz.y=rep(max.y,n.index)
		horiz.x=c(1:n.index)
		
		ver.y=seq(0,max.y+0.2,length.out=1000)
		ver.x=rep(max.x,1000)
		
		label=paste("Gene List Rank for anaylsis of COG ",cat)
		
		label=paste(label,": p-value=")
		label=paste(label,round(pval,digits=5))
		plot(delta,ylim=c(0,1),xlab=label,ylab="Running Enrichment Score",type="l",col="brown")
		points(horiz.x,horiz.y,type="l",col="red",lty=3)
		points(ver.x,ver.y,type="l",col="red",lty=3)
		points(max.x,max.y,col="red",pch=18)
		axis(4,max.y,labels=round(max.y,digits=3),lty=3,las=0)
		axis(3,max.x,labels=max.x,lty=3,las=0)
		
		Spoints.x=which(S[,index]!=1)
		Spoints.ys=rep(0,length(Spoints.x))
		Spoints.ye=rep(0.1,length(Spoints.x))
		segments(Spoints.x,Spoints.ys,Spoints.x,Spoints.ye,col="green")
		
		Spoints.x=which(S[,index]==1)
		Spoints.ys=rep(0,length(Spoints.x))
		Spoints.ye=rep(0.1,length(Spoints.x))
		segments(Spoints.x,Spoints.ys,Spoints.x,Spoints.ye,col="red")	
	}
	
}
