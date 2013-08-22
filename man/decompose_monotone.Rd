\name{decompose_monotone}
\alias{decompose_monotone}
\alias{decompose_monotone_single}
\title{Decompose genotype-phenotype map(s) using monotone regression }
\description{
Functions for decomposing genotype-phenotype (GP) maps with \eqn{N} biallelic loci using monotone regression from the \href{http://cran.r-project.org/package=isotone}{\code{isotone}} package.
}

\usage{
decompose_monotone(gpmap)
decompose_monotone_single(gpmap, parallel=FALSE)
}

\arguments{
  \item{gpmap}{ A \code{\link{gpmap}} object }
  \item{parallel}{ Indicating whether computations should be run in parallel }
}

\details{
\code{decompose_monotone} works for any \code{\link{gpmap}} object (\code{values} is \eqn{(3^N)xK} matrix of genotypic values) 
and calls the internal function \code{decompose_monotone_single} for each column. \code{decompose_monotone_single} takes a 
\code{\link{gpmap}} object with a single set of genotypic values (\eqn{K=1}), loops through all \eqn{2^N} possible
combinations of plusalleles, calls \code{\link{monotone_regression}} and identifies the best fit. If \code{parallell=TRUE} then
computations are run in parallel on multiple cores.

}

\value{

The input \code{\link{gpmap}} is returned with two added elements
  \item{monoR2}{The coefficient of determination of the monotone regression}
  \item{values.mono}{A matrix of genotypic values for the monotone component of genotype-phenotype map(s)}
}

\references{
Leeuw J, Hornik K and Mair P (2009) Isotone Optimization in R: Pool-Adjacent-Violators Algorithm (PAVA) and Active Set Methods. Journal of Statistical Software 32(5) [\href{http://www.jstatsoft.org/v32/i05}{link}]

Gjuvsland AB, Wang Y, Plahte E, Omholt SW (2013) Monotonicity is a key feature of genotype-phenotype maps. Submitted to Frontiers in Genetics
}

\author{
Arne B. Gjuvsland <arne.gjuvsland@umb.no> and Yunpeng Wang <yunpeng.wng@gmail.com>
}

\examples{
data(GPmaps)

#Additive GP map is monotone so monoR2=1 and values.mono=values
decompose_monotone(A)

#Pure AxA epistasis map
decompose_monotone(AA)

#two-locus example in Cheverud & Routman (1995)
decompose_monotone(mouseweight)

#decompose four random 3-locus GP maps
set.seed(0)
randomGP <- rnorm(3^2*4)
dim(randomGP) <- c(9,4)
decompose_monotone(generate_gpmap(randomGP))

}

