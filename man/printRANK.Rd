\name{printRANK}
\alias{printRANK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
printRANK
}
\description{
This command print the outcome of RANK(Mann-Whiteny) test. 
}
\usage{
printRANK(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	It is a BOG object.
}
}
\examples{
	bog=BOG(data="anthracis_iron",cog.file="anthracis")
	printRANK(bog)
}
