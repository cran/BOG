\name{gseaplot}
\alias{gseaplot}
\title{
gseaplot
}
\description{
This is the command to visualize the path of GSEA scores. The GSEA takes all genes into account by constructing a test statistic based on their local-fdr without pre-selection of threshold.  
}
\usage{
gseaplot(x, cat = NULL)
}
\arguments{
  \item{x}{
This is a bog object.
}
  \item{cat}{
To plot an enrichment score behavior, a user needs to specify a cat argument from one of COG groups. 
The default of cat is set NULL to raise BOG specific warning instead of R-system warning when a user fogets to specify COG, which should match one of the COGs specified in the input file or in the built-in files. The warning message is ""BOG : Category needs to be specified.""
}
}
\references{
Subramian, A. et al.(2005) Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. PNAS vol.102 no.43, 15545-15550.	
}

