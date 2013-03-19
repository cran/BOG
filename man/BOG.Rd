\name{BOG}
\alias{BOG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BOG
}
\description{
This function is a flagship function in BOG. It reads data and COG annotation files with user specified setting for analysis. It performs hypergeometric test, rank test (Mann-Whiteny test), and GSEA test (optional).
}
\usage{
BOG(data = NULL, data_type = c("data", "raw.pval", "adj.pval"), cog_file = NULL, 
conditional = FALSE, alternative = c("over.expr", "under.expr", "both.expr"), 
fdr.cutoff = 0.05, gsea = FALSE, gsea.fdr.alpha = 0.05, 
DIME.K = 5, DIME.iter = 50, DIME.rep = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
	This input file can be either a dataframe or a text file consisting of two
		columns. The first column is the geneIDs (characters). The second column provides numerical measures for the corresponding genes, which has
		three possible options controlled by the data_type argument. If data is not specified, BOG load a built-in data, anthracis_adjpval, by default.
}
  \item{data_type}{
1. data_type="data" : normalized ``differences'' of gene expressions between two comparison groups.

2. data_type="raw.pval" : raw p-values for each gene if differential analysis is carried out beforehand

3. data_type="adj.pval" : multiple testing adjusted p-values. 

If the data is specified for option(1), then DIME will be called to perform the differental analysis. 
Under option(2), the raw p-values will be converted to FDRs. If option(3) is specified, no preprocessing is needed before carrying out the tests. 

Default data-type is "data".
}
  \item{cog_file}{
		This can either be a user specified input file (R dataframe), a raw text file, or simply the 
		specification of the name of one of the six built-in COGs: anthracis, brucella, coxiella, difficile, ecoli, or francisella.
		If the virus/bateria being analyzed is not one of the six built-in varieties, then a data frame or a text file with two columns is required: the first column provides
		geneIDs as in the input data file; the second column specifies the cluster of orthologous groups to which each gene belongs. 
		BOG will perform statistical tests by first merging the data and cog_file using geneID as the key, hence it is important that geneIDs in both dataframes match. 
		If cog_file is not specified, BOG loads ''anthracis'' by default.
}
  \item{conditional}{
This argument works with ''alternative'' argument. By default, conditional is FALSE under which BOG use all the data for statistical analysis.
 
1. conditional=FALSE : alternative argument is ignored.

2. conditional=TRUE, alternative="over.expr" : statistical tests are performed only with over-expressed genes.

3. conditional=TRUE and alternative="under.expr" : statistical tests are performed only with under-expressed genes.

4. conditional=TRUE and alternative="both.expr" : statistical tests are performed with over-expressed and under-expressed genes respectively.

Default values are conditional=FALSE and alternative="over.expr".
}
  \item{alternative}{
	This argument works with ''conditional'' argument.
}
  \item{fdr.cutoff}{
In statistical analysis, BOG uses local-fdr as a score for strength of evidence for differences between groups. The smaller absolute value of local fdr, the stronger is the evidence for differences in gene expression. fdr.cutoff is a threshhold used in hypergeometric test. By default, it is set 0.05. 
}
  \item{gsea}{
By default, gsea is set FALSE so that unless user specify it to be TRUE, BOG does not perform GSEA test.
}
  \item{gsea.fdr.alpha}{
It is the alpha used in p-value calculation in GSEA. By default, gsea.fdr.alpha=0.05.
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
	dime_data is a list of the data specofied in the data argument.
}
}
\examples{
	library(hash)
	library(DIME)
	bog=BOG(data="anthracis_iron",data_type="adj.pval",cog_file="anthracis",
	conditional=FALSE,gsea=FALSE)
}

