\name{print.gpmap}
\alias{print.gpmap}

\title{ Print function for \code{\link{gpmap}} objects }

\description{ Print a summary of a genotype-phenotype (GP) map (an object of class \code{\link{gpmap}}) 
              with 1-3 \eqn{N} biallelic loci.}

\usage{\method{print}{gpmap}(x, ...)}

\arguments{
  \item{x}{ A \code{\link{gpmap}} object }
  \item{...}{ignored}
  }
  
\details{
Prints name(s) of GP map(s) and loci, a summary of genotypic values. 
Monotonicity measures are printed if available.
}

\author{
Arne B. Gjuvsland <arne.gjuvsland@umb.no> and Yunpeng Wang <yunpeng.wng@gmail.com>
}

\examples{
data(GPmaps)
print(A)
}
