\name{BOGprint}
\alias{BOGprint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BOGprint
}
\description{
This command print out all the statistical tests. 
It provides p-values and adjusted p-values (adjusted for multiple comparison by Benjamini and Hochberg method) ordered by significance. If gsea argument is set FALSE in BOG() command, only hypergeometric and Mann-Whitney (rank) tests are printed.
}
\usage{
BOGprint(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
It is a BOG object.
}
}
\examples{
	library(hash)
	library(DIME)
	bog=BOG(conditional=FALSE,gsea=FALSE)
	BOGprint(bog)
}
