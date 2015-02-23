\name{printHG}
\alias{printHG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
printHG
}
\description{
This command print the outcome of HG(Hypergeometric) test. 
}
\usage{
printHG(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	It is a BOG object.
}
}
\examples{
	bog=BOG(data="anthracis_iron",cog.file="anthracis")
	printHG(bog)
}

