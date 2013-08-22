\name{plot.gpmap}
\alias{plot.gpmap}
\alias{plot.gpmap.decomp.mono}

\title{ Functions for creating lineplots of genotype-phenotype (GP) maps }

\description{ Function for creating lineplots for genotype-phenotype (GP) map (an object of class \code{\link{gpmap}}) 
              with 1-3 biallelic loci.}

\usage{ 
\method{plot}{gpmap}(x, show=1, decomposed=FALSE, ...)
}

\arguments{
  \item{x}{ A \code{gpmap} object }
  \item{show}{ Which map (only used if >1 map in \code{gpmap} object) to plot}
  \item{decomposed}{ Decomposition into monotone and non-monotone component plotted if TRUE}
  \item{...}{ignored}
  }

\author{
Arne B. Gjuvsland <arne.gjuvsland@umb.no> and Yunpeng Wang <yunpeng.wng@gmail.com>
}

\examples{
data(GPmaps)

#plot additive gpmap
plot(A)

#plot decomposition of GPmap from Cheverud & Routman (1995)
decomp <- decompose_monotone(mouseweight)
plot(decomp, decomposed=TRUE)
}
