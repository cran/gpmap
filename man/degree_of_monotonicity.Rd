\name{degree_of_monotonicity}
\alias{degree_of_monotonicity}
\alias{degree_of_monotonicity_single}

\title{Degree of monotonicity of GP map}

\description{ Functions for computing degree of monotonicity \eqn{m} for \code{\link{gpmap}} objects. 
}
   
\usage{ 
  degree_of_monotonicity(gpmap) 
  degree_of_monotonicity_single(gpmap) 
}

\arguments{
  \item{gpmap}{ A \code{\link{gpmap}} object }
}
\details{ 
\code{degree_of_monotonicity} works for any \code{\link{gpmap}} object (\code{values} is \eqn{(3^N)xK} matrix of genotypic values) 
and calls the internal function \code{degree_of_monotonicity_single} for each column. \code{degree_of_monotonicity_single} computes
substitution effect, locus weights and per-locus and overall degree of monotonicity as decribed in Gjuvsland \emph{et al.} (2013).
}

\value{ 
\code{degree_of_monotonicity} returns the input \code{gpmap} with the following added fields:
  \item{degree.monotonicity}{Overall degree of monotonicity for the \eqn{K} GP maps}
  \item{degree.monotonicity.locus}{Data frame with per locus degree of monotonicity for the \eqn{K} GP maps}
  \item{locus.weight}{ Data frame with locus weigths}
 }

\references{ 
Gjuvsland AB, Wang Y, Plahte E and Omholt SW (2013) Monotonicity is a key feature of genotype-phenotype maps. 
Front. Genet. 4:216. doi: 10.3389/fgene.2013.00216 [\href{http://www.frontiersin.org/Journal/10.3389/fgene.2013.00216/full}{link}]

}

\author{ 
Arne B. Gjuvsland <arne.gjuvsland@nmbu.no> and Yunpeng Wang <yunpeng.wng@gmail.com>
}

\examples{
data(GPmaps)

#Additive GP map is monotone 
degree_of_monotonicity(A)

#Pure AxA epistasis map
degree_of_monotonicity(AA)

#two-locus example in Cheverud & Routman (1995)
degree_of_monotonicity(mouseweight)
}

