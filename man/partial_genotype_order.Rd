\name{partial_genotype_order}
\alias{partial_genotype_order}

\title{Generate partial ordering of genotype space based from allele content of genotypes}

\description{ Generate the strict partial order on genotype space specified in eq. (13) in Gjuvsland \emph{et al.} (2011). 
For a genotype space with \eqn{N} biallelic loci a minimal description of this partial order is given by \eqn{2N*3^{(N-1)}}
inequalities. } 
              
\usage{ partial_genotype_order(plusallele) }

\arguments{
  \item{plusallele}{ A \eqn{N}-vector with indexes (1 or 2) for the +allele for each locus }
}

\details{

In short a partial order on a set is a binary relation defining a parwise ordering of some pairs of elements in the set, for a formal definiton see \url{http://en.wikipedia.org/wiki/Partial_order}. 
In the partial order on the set of genotypes defined in Gjuvsland et al. (2011) the comparable pairs of genotypes are equal at every locus except one, 
while all other pairs of genotypes are incomparable. This partial ordering of genotype space is implicit in the regression on gene content 
(the number of alleles with a given index in each genotype) used for decomposition of tne genotypic value in quantitative genetics (see e.g. Lynch and Walsh page 65).
}
\value{
Returns a \eqn{(2N*3^{(N-1)}) \times 2} matrix of genotype indexes. The genotype indexes refer to row number in the genotype sequence set up in \code{\link{enumerate_genotypes}}.  
Each row vector in the matrix contains the genotype indexes of one comparable pairs, and if the first index is \eqn{k} and the second is \eqn{l} then genotype[\eqn{k}] < genotype[\eqn{l}]. 
}
\references{
Gjuvsland AB, Vik JO, Woolliams JA, Omholt SW (2011) Order-preserving principles underlying 
	genotype-phenotype maps ensure high additive proportions of genetic variance.
	Journal of Evolutionary Biology 24(10):2269-2279 [\href{http://onlinelibrary.wiley.com/doi/10.1111/j.1420-9101.2011.02358.x/full}{link}]

Lynch M & Walsh B (1998) Genetics and Analysis of Quantitative Traits, Sunderland, MA: Sinauer Associates
}

\author{
Arne B. Gjuvsland <arne.gjuvsland@nmbu.no> and Yunpeng Wang <yunpeng.wng@gmail.com>
}





