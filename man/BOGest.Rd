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
BOGest(data, data.type, cog.file, hg.thresh, gsea, DIME.K, DIME.iter, DIME.rep)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	The definition of all the arguments are same as the ones described in the BOG() command so that a user may refer to BOG() command for details.
  \item{data}{
  	This input file can be either a dataframe or a text file consisting of two
		columns. The first column is the geneIDs (charaters). The second column provides numerical measures for the corresponding genes, which has
		three possible options controlled by the data_type argument. If data is not specified, BOG load a built-in data, anthracis_adjpval, by default.
}
  \item{data.type}{
1. data.type="data" : normalized ``differences'' of gene expressions between two comparison groups.

2. data.type="pval" : raw p-values or multiple testing adjusted p-values for each gene if differential analysis is carried out beforehand

If the data is specified for option(1), then DIME will be called to perform the differental analysis. 
Under option(2), no preprocessing is needed before carrying out the tests. 

Default data.type is "data".
}
  \item{cog.file}{
		This can either be a user specified input file (R dataframe), a raw text file, or simply the 
		specification of the name of one of the six built-in COGs: anthracis, brucella, coxiella, difficile, ecoli, or francisella.
		If the virus/bateria being analyzed is not one of the six built-in varieties, then a data frame or a text file with two columns is required: the first column provides
		geneIDs as in the input data file; the second column specifies the cluster of orthologous groups to which each gene belongs. 
		BOG will perform statistical tests by first merging the data and cog_file using geneID as the key, hence it is important that geneIDs in both dataframes match. 
		If cog_file is not specified, BOG loads ''anthracis'' by default.
}
  \item{hg.thresh}{
In statistical analysis, BOG uses local-fdr(or p-value) as a score for strength of evidence for differences between groups. The smaller absolute value of the score, the stronger is the evidence for differences in gene expression. hg.thresh is a threshhold used in hypergeometric test. By default, it is set 0.05. 
}
  \item{gsea}{
By default, gsea is set FALSE so that unless user specify it to be TRUE, BOG does not perform GSEA test.
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
\value{
	List with three elements : stat, dime, dime_data.
\item{stat}{
	stat is a list consisting of statistical analysis outputs. 
}
\item{dime}{
	dime is a list consisting of DIME outputs.
}
\item{dime_data}{
	dime_data is a list of the data specified in the data argument.
}
}