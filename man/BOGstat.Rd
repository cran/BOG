\name{BOGstat}
\alias{BOGstat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BOGstat
}
\description{
This is an internal function so that it is not expected for a user to use it.
}
\usage{
BOGstat(db_Gene, conditional, alternative, fdr.cutoff, gsea, gsea.fdr.alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db_Gene}{
Internally defined list object.
}
  \item{conditional}{
his argument works with ''alternative'' argument. 
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
}
