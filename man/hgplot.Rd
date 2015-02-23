\name{hgplot}
\alias{hgplot}
\title{
hgplot
}
\description{
A user can visualize the hypergeometric test results using this command.
It display the most significant COG groups (adjusted p-value < 0.1) with observed and expected counts displayed. 
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
	bog=BOG(gsea=FALSE)
	hgplot(bog)
}
