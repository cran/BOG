BOG <-
function(data=NULL, data_type=c("data","raw.pval","adj.pval"), cog_file = NULL, conditional=FALSE,
alternative=c("over.expr","under.expr","both.expr"),fdr.cutoff=0.05, 
gsea=FALSE,gsea.fdr.alpha=0.05,DIME.K=5,DIME.iter=50,DIME.rep=5)
{
	set.seed(12345)
	data_type=match.arg(data_type)
	alternative=match.arg(alternative)
	
	result=BOGest(data, data_type, cog_file, conditional, alternative, fdr.cutoff=fdr.cutoff, 
	gsea=gsea, gsea.fdr.alpha=gsea.fdr.alpha, DIME.K, DIME.iter, DIME.rep)
	result$call <- match.call()
	class(result) <- "BOG"
	result
}
