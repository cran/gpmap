\name{enumerate_genotypes}
\alias{enumerate_genotypes}

\title{Function for enumerating genotypes for \eqn{N} biallelic loci}

\description{ Function for enumerating all \eqn{3^N} genotypes for \eqn{N} biallelic loci. Optional specification of names of loci and alleles. 
              genotypes. Generates a data frame of multilocus genotypes in the sequence used for objects of class \code{\link{gpmap}}.   }

\usage{ enumerate_genotypes(nloci=1, locinames=NULL, allelenames=NULL) }

\arguments{
  \item{nloci}{ The number of loci \eqn{N} }
  \item{locinames}{ An optional character vector with \eqn{N} names of loci }
  \item{allelenames}{ An optional character object specifying allele names }
}

\details{
Unless specified locinames default to "Locus 1", "Locus 2",..,"Locus N". 

If \code{allelenames} is not specified then the alleles will be named "1" and "2".
}
\value{

Returns a data frame with locinames as colnames, and with \eqn{3^N} rows specifying all possible genotypes in the sequence 
used for all GP maps in the package (the same sequence as used in Gjuvsland \emph{et al.}( 2011)), where the genotype at the first locus varies fastest,
then the second locus, and so on:

\preformatted{
  Locus_1 Locus_2
1      11      11
2      12      11
3      22      11
4      11      12
5      12      12
6      22      12
7      11      22
8      12      22
9      22      22
}

}
\references{
Gjuvsland AB, Vik JO, Woolliams JA, Omholt SW (2011) Order-preserving principles 
underlying genotype-phenotype maps ensure high additive proportions of genetic 
variance. Journal of Evolutionary Biology 24(10):2269-2279

}
\author{
Arne B. Gjuvsland <arne.gjuvsland@umb.no> and Yunpeng Wang <yunpeng.wng@gmail.com>
}

\examples{
#genotypes for a single locus
enumerate_genotypes()

#genotypes for two loci "A" and "B", with alleles named "H" and "L"
enumerate_genotypes(2,c("A","B"),c("H","L"))

#genotypes for the two-locus example in Cheverud & Routman (1995)
enumerate_genotypes(2,c("D7Mit17","D1Mit7"),rbind(c('A1','A2'),c('B1','B2')))
}


