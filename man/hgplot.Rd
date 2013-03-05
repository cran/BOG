\name{hgplot}
\alias{hgplot}
\title{
hgplot
}
\description{
A user can visualize hypergeometric test using this command.
It display the most significant COG groups (p-value < 0.1) with observed and expected counts displayed. 
}
\usage{
hgplot(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
This is a BOG object.
}
}
\examples{
	library(hash)
	library(DIME)
	bog=BOG(gsea=FALSE,gsea.fdr.alpha=0.05)
	hgplot(bog)
}
