BOGest <-
function(data, data_type, cog_file, conditional, alternative, 
fdr.cutoff, gsea,gsea.fdr.alpha,DIME.K,DIME.iter,DIME.rep)
{   
	anthracis_iron <- NULL
	anthracis <- NULL
	brucella <- NULL
	coxiella <- NULL
	difficile <- NULL
	ecoli <- NULL
	francisella <- NULL
	
	if(is.null(data)){
		cat( "Data file is not specified so that the default anthracis data is being loaded.\n")
		data(anthracis_iron,envir = environment())
		dime_data <- get("anthracis_iron",envir = environment())
		data_type="adj.pval"
	}else if(is.data.frame(data)){
		if(is.numeric(data[[1]])) stop("BOG : The fisrt element of data frame should be character")
		if(!is.numeric(data[[2]])) stop("BOG : The second elemnt of data frame should be numeric")
		dime_data=data	
	}else{
		if(tolower(data) %in% c("anthracis_iron")) {
		data(anthracis_iron,envir = environment())
		dime_data <- get("anthracis_iron",envir = environment())
		data_type="adj.pval"
		}else dime_data <- read.table(file=data,header = TRUE)
	}
	
	if(is.null(cog_file)){
		cat( "COG file is not specified so that the default Anthracis COG data is uploaded.\n")
		data(anthracis,envir = environment())
		db_cogs <- get("anthracis",envir = environment())
	}else if(is.data.frame(cog_file)){
			db_cogs=cog_file
	}else{	
			if(tolower(cog_file) %in% c("anthracis","brucella","coxiella","difficile","ecoli","francisella")){
				switch(tolower(cog_file),
				"anthracis"={data(anthracis,envir = environment());db_cogs <- get("anthracis",envir = environment())},
				"brucella"={data(brucella,envir = environment());db_cogs <- get("brucella",envir = environment())},
				"coxiella"={data(coxiella,envir = environment());db_cogs <- get("coxiella",envir = environment())},
				"difficile"={data(difficile,envir = environment());db_cogs <- get("difficile",envir = environment())},
				"ecoli"={data(ecoli,envir = environment());db_cogs <- get("ecoli",envir = environment())},
				"francisella"={data(francisella,envir = environment());db_cogs <- get("francisella",envir = environment())},
				)
			}else {
				cat("Reading COG information from",cog_file,"\n")	
				db_cogs <- read.table(file=cog_file,header = TRUE)
			}
  	}
		
	all_cog_type=LETTERS[c(1:26)]
	
	if(length(db_cogs)!=2) stop("BOG : COG file should be two columns")
	
	#### Local fdr ################
#############  DATA TYPE
	dime_list=list(dime_data[[2]])
	
	k=DIME.K
	iter=DIME.iter
	rep=DIME.rep
		
	switch(data_type,
		"data"={	
			fit_dime=DIME(dime_list,gng.K=k,inudge.K=k,gng.max.iter=iter,gng.rep=rep,
				inudge.max.iter=iter,inudge.rep=rep,nudge.max.iter=iter,nudge.rep=rep)
			best_fdr=fit_dime$best$fdr
			best_fdr=best_fdr/max(best_fdr)
		},
		"raw.pval"={
			if(min(unlist(dime_list))<0) stop("There is a negative p-value")
			if(max(unlist(dime_list))>1) stop("There is a p-value greater than 1")
			pval_sort=sort(unlist(dime_list),decreasing=FALSE,index.return=TRUE)
			adj_pval=p.adjust(pval_sort$x,method="BH")
			best_fdr=adj_pval[order(pval_sort$ix)]
			fit_dime <- NULL
		},
		"adj.pval"={
			if(min(unlist(dime_list))<0) stop("There is a negative p-value")
			if(max(unlist(dime_list))>1) stop("There is a p-value greater than 1")
			best_fdr=unlist(dime_list)
			fit_dime <- NULL
		}
	)
	
	
########################	
	best_fdr[best_fdr<1e-12]=1e-12;
	signed_fdr=best_fdr*sign(dime_data[[2]])
	
	f.dime=as.data.frame(list(dime_data,best_fdr=best_fdr,signed_fdr=signed_fdr))
	f.cog=as.data.frame(db_cogs)
	f.dc=merge(f.dime,f.cog,by=names(f.cog)[1])
	
	quality.data=f.dc[substr(f.dc[,5],1,3)=="COG",]
	
	###########quality_data = (geneID, log_ratio, best_fdr, signed_fdr, COG)
	###########dime = (geneID, best_fdr, signed_fdr)
	dime=list(quality.data[,1],quality.data[,3],quality.data[,4])
	
	################################
	
	index=c(1:length(dime[[1]]))
	deGene <- dime[[1]] 
	de_cogs <- quality.data[,5]
	tmp <- unlist(deGene);
	expr <- tmp 

	if(conditional==TRUE){
		min.log=min(-log(dime[[2]]))
		s = (-log(dime[[2]])-min.log) /max(-log(dime[[2]])-min.log)*sign(dime[[3]])
		db_Gene <- list(index=index,id=deGene,cog=de_cogs,fdr= dime[[3]],gsea=s)
	
		if(alternative=="both.expr"){
			stat1=BOGstat(db_Gene,conditional,"over.expr",fdr.cutoff,gsea,gsea.fdr.alpha)
			stat2=BOGstat(db_Gene,conditional,"under.expr",fdr.cutoff,gsea,gsea.fdr.alpha)
			stat=list(stat1,stat2)
		}else stat=BOGstat(db_Gene,conditional,alternative,fdr.cutoff,gsea,gsea.fdr.alpha)	
	}else{
			min.log=min(-log(dime[[2]]))
			s = (-log(dime[[2]])-min.log) /max(-log(dime[[2]])-min.log)
			db_Gene <- list(index=index,id=deGene,cog=de_cogs,fdr= dime[[2]],gsea=s)
			stat=BOGstat(db_Gene,conditional,alternative,fdr.cutoff,gsea,gsea.fdr.alpha)
	}
	
	list(stat=stat,dime=fit_dime,dime_data=dime_list)	
}
