\name{mscoref}
\alias{mscoref}
\title{
Computes M-scores for a Full Match
}
\description{
Of limited interest to most users, mscoref() computes the
scores that form the basis for the hypothesis test
performed by senfm.
}
\usage{
mscoref(ymat, treated1, inner = 0, trim = 3, qu = 0.5)
}
\arguments{
  \item{ymat}{
If there are I matched sets and the largest matched set contains J individuals,
then y is an I by J matrix with one row for each matched set.  If matched
set i contains one treated individual and k controls, where k is at least 1 and
at most J-1, then y[i,1] is the treated individual's response, y[i,2],...,y[i,k+1]
are the responses of the k controls, and y[i,k+2],...,y[i,J] are equal to NA.
If matched set i contains one control and k>1 treated individuals, then y[i,1]
is the control's response, y[i,2],...,y[i,k+1] are the responses of the k
treated individuals, and y[i,k+2],...,y[i,J] are equal to NA.
}
  \item{treated1}{
The vector treated1 is a logical vector of length I,
where treated1[i]=TRUE if there is one
treated subject in matched set i and treated1[i]=FALSE if there is more than one
treated subject in matched set i.
}
  \item{inner}{
inner and trim together define the \eqn{\psi}-function for the M-statistic.
The default values yield a version of Huber's \eqn{\psi}-function, while
setting inner = 0 and trim = Inf uses the mean within each matched set.
The \eqn{\psi}-function is an odd function, so \eqn{\psi(w) = -\psi(-w)}.  For
w >= 0, the \eqn{\psi}-function is \eqn{\psi(w)=0} for 0 <= w <= inner,
is \eqn{\psi(w)} = trim for w >= trim, and rises linearly from 0 to trim
for inner < w < trim.
}
  \item{trim}{
inner and trim together define the \eqn{\psi}-function for the M-statistic.
See inner.
}
  \item{qu}{
Before applying the \eqn{\psi}-function to treated-minus-control differences,
the differences are scaled by dividing by the qu quantile of all
within set absolute differences.  Typically, qu = 1/2 for the median.
The value of qu has no effect if trim=0 and inner=Inf.  See Maritz (1979)
for the paired case and Rosenbaum (2007) for matched sets.
}
}

\value{
Returns a matrix with the same dimensions as ymat and the same
pattern of NAs.  The returned value in position [i,j] compares
ymat[i,j] to the other observations in row i of ymat, scoring
the differences using \eqn{\psi}-function, totalling them, and applying
a weight.  Matched sets of unequal size are weighted
using weights that would
be efficient in a randomization test under a simple model with
additive set and treatment effects and errors with constant variance.
See Rosenbaum (2007, section 4.2).

When a matched set contains one control and several treated subjects,
this is reflected in the returned scores by a sign reversal.
}
\references{
Huber, P. (1981) Robust Statistics.  New York: John Wiley.

Maritz, J. S. (1979). A note on exact robust confidence
intervals for location. Biometrika 66 163--166.

Rosenbaum, P. R. (2007). Sensitivity analysis for m-estimates, tests
and confidence intervals in matched observational studies.
Biometrics 63 456-64. (R package sensitivitymv)

Rosenbaum, P. R. (2010). Design of Observational Studies.  New York:
Springer.  Table 2.12, page 60, illustrates the calculations for
the simple case of matched pairs.

Rosenbaum, P. R. (2013). Impact of multiple matched controls on
design sensitivity in observational studies. Biometrics 69 118-127.
(Introduces inner trimming to increase design sensitivity.)
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

mscoref(y,treated1) # Huber scores
mscoref(y,treated1,inner=0.5,trim=3) #inner trimmed scores
mscoref(y,treated1,qu=.9,trim=1) #trimming the outer 10 percent

# For an additional example, install and load package sensitivitymv
# The following example is a match with variable controls.
# Both mscorev() (in sensitivitymv) and mscoref() (in sensitivityfull)
# reproduce the example in Rosenbaum (2007, Table 3).
# data(tbmetaphase)
# mscorev(tbmetaphase,trim=1)
# mscoref(tbmetaphase,rep(TRUE,15),trim=1)

}

