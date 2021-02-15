\name{monotone_regression}
\alias{monotone_regression}

\title{ Perform monotone regression on a genotype-phenotype (GP) map }
\description{
The function uses \code{\link{partial_genotype_order}} and \code{\link[isotone]{activeSet}} from the \href{https://cran.r-project.org/package=isotone}{isotone} package to do 
monotone regression (Leeuw \emph{et al.}, 2009) on a GP map.
}

\usage{
monotone_regression(gpmap, plusallele)
}

\arguments{
  \item{gpmap}{ An object of class \code{\link{gpmap}} }
  \item{plusallele}{An \eqn{N} vector of allele indexes (1 or 2)}
}

\details{
Element \emph{i} in \code{plusallele} specifies the ordering of the genotypes at locus \emph{i}, if the element is 1 then 11 < 12 < 22 and conversely if it is 2 then
22 < 12 < 11. \code{monotone_regression} calls \code{\link{partial_genotype_order}} to obtain the partial ordering of genotypic values for the given \code{plusalleles}.
This partial ordering is then used together with the GP map itself as input to the \code{\link[isotone]{activeSet}} function from the package 
\href{https://cran.r-project.org/package=isotone}{isotone}.


}
\value{
\code{monotone_regression} returns the output from \code{\link[isotone]{activeSet}} directly. 
}
\references{
Leeuw J, Hornik K and Mair P (2009) Isotone Optimization in R: Pool-Adjacent-Violators 
Algorithm (PAVA) and Active Set Methods. Journal of Statistical Software 32(5) [\href{https://www.jstatsoft.org/v32/i05}{link}]

Gjuvsland AB, Wang Y, Plahte E and Omholt SW (2013) Monotonicity is a key feature of genotype-phenotype maps. 
Front. Genet. 4:216. doi: 10.3389/fgene.2013.00216 [\href{https://www.frontiersin.org/articles/10.3389/fgene.2013.00216/full}{link}]
}

\author{
Arne B. Gjuvsland <arne.gjuvsland@nmbu.no> and Yunpeng Wang <yunpeng.wng@gmail.com>
}

\examples{
data(GPmaps)

#Additive GP map is monotone 
monotone_regression(A,c(2,2))

#Pure AxA epistasis map
monotone_regression(AA,c(2,2))

#two-locus example in Cheverud & Routman (1995)
monotone_regression(mouseweight,c(1,1))
}
