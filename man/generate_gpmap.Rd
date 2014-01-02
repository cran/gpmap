\name{generate_gpmap}
\alias{generate_gpmap}
\alias{gpmap}

\title{ Function for creating genotype-phenotype (GP) maps }

\description{ Function for creating a \code{\link{gpmap}} object representing a genotype-phenotype (GP) map 
              for \eqn{N} biallelic loci or more generally \eqn{K} such maps, from a matrix of genotypic values. }

\usage{ generate_gpmap(y, locinames = NULL, allelenames = NULL, mapnames = NULL) }

\arguments{
  \item{y}{ A \eqn{(3^N\times K)}{(3^NxK)} matrix or numeric with each column specifying the \eqn{3^N} genotypic values for \eqn{K} GP maps }
  \item{locinames}{ An optional character vector with \eqn{N} names of loci }
  \item{allelenames}{ An optional character object specifying allele names }
  \item{mapnames}{ An optional character vector with \eqn{K} names of GP maps / phenotypes }
}

\details{
Arguments \code{locinames} and \code{allelenames} are passed on to \code{\link{enumerate_genotypes}}, and the genotypic values in \code{y} should be given 
in the same sequence as the sequence of genotypes returned by \code{\link{enumerate_genotypes}}. If \code{mapnames} is not specified then the GP maps 
will be named "GPmap_1", "GPmap_2",..,"GPmap_K".
}
\value{
The function returns an object of class \code{\link{gpmap}} containing the following components
\item{values}{ The vector or matrix of genotypic values }
\item{nloci}{ The number of loci in the map }
\item{genotypes}{ Data frame with enumeration of genotypes }
\item{locinames}{ Character vector with names for all loci }
\item{mapname}{ The name(s) of the GP map }
}
\references{
Gjuvsland AB, Vik JO, Woolliams JA, Omholt SW (2011) Order-preserving principles underlying 
	genotype-phenotype maps ensure high additive proportions of genetic variance.
	Journal of Evolutionary Biology 24(10):2269-2279 [\href{http://onlinelibrary.wiley.com/doi/10.1111/j.1420-9101.2011.02358.x/full}{link}]
}
\author{
Arne B. Gjuvsland <arne.gjuvsland@nmbu.no> and Yunpeng Wang <yunpeng.wng@gmail.com>
}

\examples{
#inter- and intra-locus additive GPmap with two loci 
generate_gpmap(c(-2,1,0,-1,0,1,0,1,2),mapnames="Additive") 

#two-locus GP map with AxA epistasis for loci named A and B
generate_gpmap(c(-1,0,1,0,0,0,1,0,-1),locinames=c("A","B"),mapnames="AxA")

#random GP map with 3 loci
set.seed(0)
generate_gpmap(rnorm(27))
}
