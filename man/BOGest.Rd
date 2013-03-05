\name{BOGest}
\alias{BOGest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BOGest
}
\description{
This is an internal function so that it is not expected for a user to use it.
}
\usage{
BOGest(data, data_type, cog_file, conditional, alternative, fdr.cutoff, 
gsea, gsea.fdr.alpha, DIME.K, DIME.iter, DIME.rep)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	The definition of all the arguments are same as the ones described in the BOG() command so that a user may refer to BOG() command for details.
  \item{data}{
  	This input file can be either a dataframe or a text file consisting of two
		columns. The first column is the geneIDs (charaters). The second column provides numerical measures for the corresponding genes, which has
		three possible options controlled by the data_type argument. If data is not specified, BOG load a built-in data, anthracis_adjpval, by default.
}
  \item{data_type}{
1. data_type="data" : normalized ``differences'' of gene expressions between two comparison groups.

2. data_type="raw.pval" : raw p-values for each gene if differential analysis is carried out beforehand

3. data_type="adj.pval" : multiple testing adjusted p-values. 
}
  \item{cog_file}{
This can either be a user specified input file (R dataframe), a raw text file, or simply the 
		specification of the name of one of the six built-in COGs: anthracis, brucella, coxiella, difficile, ecoli, or francisella.
}
  \item{conditional}{
For the detail, refer to description in BOG()
}
  \item{alternative}{
For the detail, refer to description in BOG()
}
  \item{fdr.cutoff}{
In statistical analysis, BOG uses local-fdr as a score for strength of evidence for differences between groups. The smaller absolute value of local fdr, the stronger is the evidence for differences in gene expression. fdr.cutoff is a threshhold used in hypergeometric test. By default, it is set 0.05. 
}
  \item{gsea}{
By default, gsea is set FALSE so that unless user specify it to be TRUE, BOG does not perform GSEA test.
}
  \item{gsea.fdr.alpha}{
It is the alpha used in p-value calculation in GSEA. By default, gsea.fdr.alpha=0.05
}
  \item{DIME.K}{
The number of mixture components in fitting an ensemble of mixture models. The default is 5.
}
  \item{DIME.iter}{
The number of iterations in fitting an ensemble of mixture models. The default is 50.
}
  \item{DIME.rep}{
The number of repitions in fitting an ensemble of mixture models. The default is 5.
}
}