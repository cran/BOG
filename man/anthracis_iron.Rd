\name{anthracis_iron}
\alias{anthracis_iron}
\docType{data}
\title{
anthracis_iron
}
\description{
This is the anthracis data of (Carson et.al,2009) with adjusted p-value produced by DIME. 
It consists of anthracis geneID and adjusted p-value for gene expression between two conditions. 
If a user does not specify the data argument in BOG() command, this built-in data will be loaded by default as an example.
}
\usage{data(anthracis_iron)}
\format{
  A data frame with 5537 observations on the following 2 variables: geneID, adj.pval
}
\references{
	Carlson, P. et al. (2009) Transcriptional Profiling of Bacillus anthracis Sterne (34F2) during Iron Starvation. PLOS ONE 4(9): e6988. Doi:10.1371/journal.pone.0006988.
}
\examples{
data(anthracis_iron)
}
