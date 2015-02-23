\name{BOG}
\alias{BOG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BOG
}
\description{
This function is the flagship function of BOG. It reads data and COG annotation files with user specified setting for analysis. It performs hypergeometric test, rank test (Mann-Whiteny test), and GSEA test (optional).
}
\usage{
BOG(data = NULL, data.type = c("data", "pval"), cog.file = NULL, 
hg.thresh = 0.05, gsea = FALSE,
DIME.K = 5, DIME.iter = 50, DIME.rep = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
	This input file can either be a dataframe or a text file consisting of two
		columns. The first column is the geneIDs (characters). The second column provides numerical measures for the corresponding genes, which has
		three possible options controlled by the data_type argument. If data is not specified, BOG load a built-in data, anthracis_adjpval, by default, as an example data set.
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
In the statistical analysis, BOG uses local-fdr(or p-value) as a score for strength of evidence for differences between groups. The smaller value of the score, the stronger is the evidence for differences in gene expression. hg.thresh is a threshhold used in hypergeometric test. By default, it is set to be 0.05. 
}
  \item{gsea}{
By default, gsea is set FALSE so that unless user specify it to be TRUE, BOG does not perform GSEA test.
}
   \item{DIME.K}{
The number of mixture components in fitting an ensemble of mixture models, if DIME processing is activated. If user select data.type="pval", this is irrelevant. The default is 5.
}
  \item{DIME.iter}{
The number of iterations in fitting an ensemble of mixture models. This is only relevant if data.type="data". The default is 50. 
}
  \item{DIME.rep}{
The number of repetitions in fitting an ensemble of mixture models. This is only relevant if data.type="data". The default is 5. 
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
\examples{
	bog=BOG(data="anthracis_iron",data.type="pval",cog.file="anthracis",gsea=FALSE)
}

