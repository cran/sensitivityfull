\name{separable1f}
\alias{separable1f}
\title{
Computes the Separable Approximation.
}
\description{
Of limited interest to most users, separable1f() is called by
the main function, senfm(), in computing the sensitivity analysis.
separable1f() is given scores produced by mscoref() and computes the
separable approximation to the upper bound on the P-value.
}
\usage{
separable1f(ymat, gamma = 1)
}
\arguments{
  \item{ymat}{
A matrix of scores produced by mscoref.
}
  \item{gamma}{
The sensitivity parameter \eqn{\Gamma \ge 1}.
}
}
\details{
See Gastwirth, Krieger and Rosenbaum (2000) and
Rosenbaum (2007, section 4) for discussion of the
separable approximation.
}
\value{
\item{pval }{The upper bound on the 1-sided P-value.}
\item{deviate }{The deviate that was compared to the Normal
distribution to produce pval.}
\item{statistic }{The value of the score-statistic.}
\item{expectation }{The maximum expectation of the score-statistic for the given \eqn{\Gamma}.}
\item{variance }{The maximum variance of the score-statistic among treatment assignments
that achieve the maximum expectation.  Part of the separable approximation.}
}
\references{
Gastwirth, J. L., Krieger, A. M. and Rosenbaum, P. R. (2000).
Asymptotic separability in sensitivity analysis.
J. Roy. Statist. Soc. B. 62 545-555.

Rosenbaum, P. R. (2007). Sensitivity analysis for m-estimates, tests
and confidence intervals in matched observational studies.
Biometrics 63 456-64. (See section 4.) (R package sensitivitymv)
}
\author{
Paul R. Rosenbaum
}


\examples{
# The artificial example that follows has I=9
# matched sets.  The first 3 sets have one treated
# individual and two controls with treated subjects
# in column 1.  The next 3 sets are
# matched pairs, with treated subjects in column 1.
# The next 3 sets have one control and two treated
# subjects, with the control in column 1.  Simulated
# from a Normal distribution with an additive effect
# of tau=1.

y<-c(2.2, 1.4, 1.6, 2.4, 0.7, 1.3, 1.2, 0.6, 0.3,
0.5, -0.1, -1.3, -0.3, 0.1, 0.4, 3.0, 1.1, 1.4, -0.8,
0.1, 0.8, NA, NA, NA, 1.1, 0.5, 1.8)
y<-matrix(y,9,3)
treated1<-c(rep(TRUE,6),rep(FALSE,3))

# The same calculation done in two equivalent ways.
separable1f(mscoref(y,treated1),gamma=2)
senfm(y,treated1,gamma=2)

}


