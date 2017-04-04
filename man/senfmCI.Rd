\name{senfmCI}
\alias{senfmCI}
\concept{Sensitivity analysis}
\concept{Full matching}
\concept{Observational study}
\concept{Causal inference}
\concept{Confidence interval}
\title{
Sensitivity Analysis for a Confidence Interval in a Full Match.
}
\description{
In a full match, each matched set contains either one treated
individual and one or more controls or one control and one or
more treated individuals.  Uses Huber's M-statistic as the basis for
a confidence interval for an additive constant treatment effect,
\eqn{\tau}.  Performs either a randomization
inference or an analysis of sensitivity to departures from random
assignment.  The confidence interval inverts the test in the function
senfm() in the sensitivityfull package.
}
\usage{
senfmCI(y,treated1,gamma=1,inner=0,trim=3,lambda=1/2,
                  alpha=0.05,twosided=TRUE,upper=TRUE)
}
\arguments{
  \item{y}{
If there are I matched sets and the largest matched set contains J individuals,
then y is an I by J matrix with one row for each matched set.  If matched
set i contains one treated individual and k controls, where k is at least 1 and
at most J-1, then y[i,1] is the treated individual's response, y[i,2],...,y[i,k+1]
are the responses of the k controls, and y[i,k+2],...,y[i,J] are equal to NA.
If matched set i contains one control and k>1 treated individuals, then y[i,1]
is the control's response, y[i,2],...,y[i,k+1] are the responses of the k
treated individuals, and y[i,k+2],...,y[i,J] are equal to NA.  Although y
may, and typically does, contain NA's, y[i,1] and y[i,2] must not be NA for
all i.

If you have matched pairs, not matched sets, use the senmwCI function in the
sensitivitymw package.
}
  \item{treated1}{
The vector treated1 is a logical vector of length I,
where treated1[i]=TRUE if there is one
treated subject in matched set i and treated1[i]=FALSE if there is more than one
treated subject in matched set i.
}
  \item{gamma}{
gamma is the sensitivity parameter \eqn{\Gamma}, where \eqn{\Gamma \ge 1}.  Setting
\eqn{\Gamma = 1} is equivalent to assuming ignorable treatment assignment given
the matched sets, and it returns a randomization-based confidence interval.
}
  \item{inner}{
inner and trim together define the \eqn{\psi}-function for the M-statistic.
The default values yield a version of Huber's \eqn{\psi}-function, while
setting inner = 0 and trim = Inf uses the mean within each matched set.
The \eqn{\psi}-function is an odd function, so \eqn{\psi(w) = -\psi(-w)}.  For
\eqn{w \ge 0}, the \eqn{\psi}-function is \eqn{\psi(w)=0} for \eqn{0 \le w \le } inner,
is \eqn{\psi(w)= } trim for \eqn{w \ge } trim, and rises linearly from 0 to trim
for inner < w < trim.

An error will result unless \eqn{0 \le } inner \eqn{ \le } trim.

Taking trim < Inf limits the influence of outliers; see Huber (1981).
Taking inner > 0 often increases design sensitivity; see Rosenbaum (2013).
}
  \item{trim}{
inner and trim together define the \eqn{\psi}-function for the M-statistic.
See inner.
}
  \item{lambda}{
Before applying the \eqn{\psi}-function to treated-minus-control differences,
the differences are scaled by dividing by the lambda quantile of all
within set absolute differences.  Typically, lambda = 1/2 for the median.
The value of lambda has no effect if trim=0 and inner=Inf.  See Maritz (1979)
for the paired case and Rosenbaum (2007) for matched sets.

An error will result unless 0 < lambda < 1.
}
  \item{alpha}{
The coverage rate of the confidence interval is 1-alpha. If the bias in
treatment assignment is at most \eqn{\Gamma}, then the confidence interval
will cover the true \eqn{\tau} with probability at least \eqn{1-\alpha}.
}
  \item{twosided}{
If twosided==TRUE, then a two-sided \eqn{1-\alpha} confidence interval is
constructed.  If twosided==FALSE, then a one-sided \eqn{1-\alpha}
confidence interval is constructed.  The two sided interval is the
intersection of two one-sided \eqn{1-\alpha/2} intervals.
  }
  \item{upper}{
If twosided==TRUE, then upper is ignored. If twosided==FALSE and
upper=TRUE, then an upper \eqn{1-\alpha} confidence interval
is constructed.  If twosided==FALSE and upper=FALSE, then
a lower \eqn{1-\alpha} confidence interval is constructed.
  }
}
\details{
For the given \eqn{\Gamma}, senfmCI() inverts the test in the function
senfm() to produce the confidence interval.  That is, it tests every
\eqn{\tau} and retains the values not rejected at level \eqn{\alpha}.

The test is a large sample approximation
based on a Normal approximation to the null distribution;
see Rosenbaum (2007).

Matched sets of unequal size are weighted using weights that would
be efficient in a randomization test under a simple model with
additive set and treatment effects and errors with constant variance;
see Rosenbaum (2007).

The upper bound on the P-value is based on the separable approximation
described in Gastwirth, Krieger and Rosenbaum (2000); see also
Rosenbaum (2007).
}
\value{
\item{PointEstimates }{The interval of point estimates of \eqn{\tau}.
If gamma=1, then the interval is a single point estimate.}
\item{ConfidenceInterval }{The confidence interval for \eqn{\tau}.}
\item{description }{Reminder of the coverage rate, gamma, and
type of interval.}
}


\references{

Hansen, B. B. (2007). Optmatch. R News 7 18-24. (R package optmatch)
(Optmatch can create an optimal full match.)

Hansen, B. B. and Klopfer, S. O. (2006). Optimal full matching and
related designs via network flows. J. Comput. Graph.
Statist. 15 609-627. (R package optmatch)

Huber, P. (1981) Robust Statistics.  New York: John Wiley.
(M-estimates based on M-statistics.)

Maritz, J. S. (1979). A note on exact robust confidence
intervals for location. Biometrika 66 163--166. (Introduces
exact permutation tests based on M-statistics by redefining
the scaling parameter.)

Rosenbaum, P. R. (1991). A characterization of optimal designs for
observational studies. J. Roy. Statist. Soc. B 53 597-610.
(Introduces full matching.)

Rosenbaum, P. R. (2007). Sensitivity analysis for m-estimates, tests
and confidence intervals in matched observational studies.
Biometrics 63 456-64. (R package sensitivitymv)

Rosenbaum, P. R. (2013). Impact of multiple matched controls on
design sensitivity in observational studies. Biometrics 69 118-127.
(Introduces inner trimming.)

Rosenbaum, P. R. (2014). Weighted M-statistics with superior design
sensitivity in matched observational studies with multiple controls.
J. Am. Statist. Assoc. 109 1145-1158. (R package sensitivitymw)

}
\author{
Paul R. Rosenbaum.
}
\note{
In a sensitivity analysis, a one-sided confidence interval is
not conservative; however, two-sided intervals formed
as the intersection of two one-sided \eqn{1-\alpha/2}
intervals are somewhat conservative.  See the
discussion of two-sided tests in the documentation
for senfm().
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

# Randomization interval and point estimate, Huber scores:
senfmCI(y,treated1)

# Uses senfm() to show how senfmCI() inverts the test.
senfm(y,treated1,tau=0.6172307) #P-value is 0.025
senfm(y,treated1,tau=2.0612746,alternative = "less") #P-value is 0.025
senfm(y,treated1,tau=1.345436) #Statistic is 0

senfmCI(y,treated1,gamma=1.5) #Sensitivity of two-sided CI
# The next two calculations relate one and two-sided intervals
senfmCI(y,treated1,gamma=1.5,twosided=FALSE,upper=TRUE,alpha=0.025)
senfmCI(y,treated1,gamma=1.5,twosided=FALSE,upper=FALSE,alpha=0.025)

# If an estimator is approximately Normal, then +/- a standard
# error is approximately a 2/3 confidence interval.  Going the
# other way, people sometimes suggest looking at a 2/3
# confidence interval as analogous to +/- a standard error.
senfmCI(y,treated1,gamma=1.5,alpha=1/3)

# For an additional example, install and load package sensitivitymw
# library(sensitivitymw)
# The mercury data is 397 triples, 1 treated, 2 controls.
# It is the example in Rosenbaum (2014).
# data(mercury)
# help(mercury)
# In this balanced design, senmwCI() and senfmCI() give the same CI.
# senmwCI(mercury,gamma=3)
# senfmCI(mercury,rep(TRUE,397),gamma=3,twosided=FALSE)
}

\keyword{ htest }
\keyword{ robust }
