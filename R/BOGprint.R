BOGprint <-
function(x)
{
	cat("Call:\n")
	print(x$call)
	cat("\n")
	
	if(length(x[[1]])==2){
		stat=x[[1]][[1]]
		stat2=x[[1]][[2]]
	}else{
		stat=x[[1]]
	}
	
	op=getOption("scipen")
	options(scipen=7)
	conditional=stat$conditional
	if(conditional==TRUE || length(x[[1]])==2){
		cat("Conditional Analysis: ")
		if(length(x[[1]])==2)
		{
			cat("Both-Side Alternatives are tested.\n\n")
			cat("Over-Expression Alternative\n")
			cat("\n HyperGeometric Test\n")
			print(stat$hyper)
			cat("\n")
			
			cat("Rank Test\n")
			print(stat$rank)
			cat("\n")
			
			if(!is.null(stat$gsea)){
				cat("GSEA Test\n")
				print(stat$gsea)
				cat("\n")
			}
			
			cat("\n")
			cat("Under-Expression Alternative\n")
			cat("HyperGeometric Test\n")
			print(stat2$hyper)
			cat("\n")
			
			cat("Rank Test\n")
			print(stat2$rank)
			cat("\n")
			
			if(!is.null(stat2$gsea)){
				cat("GSEA Test\n")
				print(stat2$gsea)
				cat("\n")
			}
		}else{	
			check=charmatch(stat$alternative,"over.expr",nomatch=0)
			cat("One-Side Alternative is tested.\n\n")
			if(check==1) cat("Over-Expression Alternative\n")
			else cat("Under-Expression Alternative\n")
			
			cat("\n HyperGeometric Test\n")
			print(stat$hyper)
			cat("\n")
			
			cat("Rank Test\n")
			print(stat$rank)
			cat("\n")
			
			if(!is.null(stat$gsea)){
				cat("GSEA Test\n")
				print(stat$gsea)
				cat("\n")
			}
		}
	}else{
			cat("Unconditional Analysis\n\n")
			cat("HyperGeometric Test\n")
			print(stat$hyper)
			cat("\n")
			
			cat("Rank Test\n")
			print(stat$rank)
			cat("\n")
			
			if(!is.null(stat$gsea)){
				cat("GSEA Test\n")
				print(stat$gsea)
				cat("\n")
			}
	}
	options(scipen=op)
}
