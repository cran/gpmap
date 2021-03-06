\name{gpmap-package}
\alias{gpmap-package}
\concept{genotype-phenotype map}
\concept{monotonicity}
\concept{functional epistasis}
\docType{package}
\title{
  Tools for analysing and plotting genotype-phenotype maps}
\description{
  The \pkg{gpmap} package deals with genotype-phenotype maps for biallelic loci underlying 
  quantitative phenotypes. The package provides a class \code{gpmaps}, analysis functions and basic
  lineplots. The package is designed for studying the properties of GP maps without
  reference to any particular population, i.e. the physiological (Cheverud & Routman, 1995)
  or functional (Hansen, 2001) properties of the GP map. This is opposed to statistical effects underlying 
  most of quantitative genetics, where the GP-map is analysed togehter with genotype frequencies in a given 
  population (e.g. Lynch & Walsh, 1998).

  In version 0.1 which is released as part of the publication of Gjuvsland \emph{et al.} (2013) 
  we have implemented functionality for studying monotonicity Gjuvsland \emph{et al.} (2011) of GP maps. 
  The package utilizes the \href{https://cran.r-project.org/package=isotone}{\pkg{isotone}} package for monotone regression, 
  and the \href{https://cran.r-project.org/package=foreach}{\pkg{foreach}} package for parallel computation.
  
   
  The package consists of the following high-level functions :
  \code{\link{enumerate_genotypes}}, \code{\link{generate_gpmap}},  
  \code{\link{degree_of_monotonicity}}, \code{\link{decompose_monotone}} and
  \code{\link{plot.gpmap}}
}

\author{
    Arne B. Gjuvsland <arne.gjuvsland@nmbu.no> and Yunpeng Wang <yunpeng.wng@gmail.com>
}

\references{

Cheverud JM & Routman EJ (1995) Epistasis and Its Controbution to Genetic Variance Components. 
	Genetics 139:1455-1461 [\href{https://www.genetics.org/content/139/3/1455.abstract}{link}]

Gjuvsland AB, Vik JO, Woolliams JA, Omholt SW (2011) Order-preserving principles underlying 
	genotype-phenotype maps ensure high additive proportions of genetic variance.
	Journal of Evolutionary Biology 24(10):2269-2279 [\href{https://onlinelibrary.wiley.com/doi/full/10.1111/j.1420-9101.2011.02358.x}{link}]

Gjuvsland AB, Wang Y, Plahte E and Omholt SW (2013) Monotonicity is a key feature of genotype-phenotype maps. 
Front. Genet. 4:216. doi: 10.3389/fgene.2013.00216 [\href{https://www.frontiersin.org/articles/10.3389/fgene.2013.00216/full}{link}]

Hansen T & Wagner GP (2001) Modeling genetic Architecture: A Multilinear Theory of gene Interaction. 
	Theoretical Population Biology 59:61-86 [\href{https://www.sciencedirect.com/science/article/pii/S0040580900915089}{link}]

Leeuw J, Hornik K and Mair P (2009) Isotone Optimization in R: Pool-Adjacent-Violators 
Algorithm (PAVA) and Active Set Methods. Journal of Statistical Software 32(5) [\href{https://www.jstatsoft.org/v32/i05}{link}]

Lynch M & Walsh B (1998) Genetics and Analysis of Quantitative Traits, Sunderland, MA: Sinauer Associates
}

\keyword{package}
