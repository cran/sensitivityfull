\name{senfm}
\alias{senfm}
\concept{Sensitivity analysis}
\concept{Full matching}
\concept{Observational study}
\concept{Causal inference}
\title{
Sensitivity analysis for a full match in an observational study.
}
\description{
In a full match, each matched set contains either one treated
individual and one or more controls or one control and one or
more treated individuals.  Uses Huber's M-statistic as the basis for
the test, for instance, a mean.  Performs either a randomization
test or an analysis of sensitivity to departures from random
assignment.
}
\usage{
senfm(y, treated1, gamma = 1, inner = 0, trim = 3, lambda = 1/2,
     tau = 0, alternative="greater")
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
treated individuals, and y[i,k+2],...,y[i,J] are equal to NA.

Exception: If y is a vector, not a matrix, then every matched set is understood
to be a matched pair, one treated subject and one control, and y is understood
to be a vector of I=length(y) matched pair differences, treated-minus-control.
The same results would be produced by setting y equal to an I-by-2 matrix whose
first column contains treated responses and whose second column contains control
responses such that the difference between the first and second column equals
the vector y.  For matched pairs, the package sensitivitymw does more, but
sensitivityfull can be used.
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
the matched sets, and it performs a within-set randomization test.
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
  \item{tau}{
The null hypothesis asserts that the treatment has an additive effect, tau.
By default, tau=0, so by default the null hypothesis is Fisher's sharp null
hypothesis of no treatment effect.
}
  \item{alternative}{
If alternative="greater", the null hypothesis of a treatment effect of tau
is tested against the alternative of a treatment effect larger than tau.  If
alternative="less", the null hypothesis of a treatment effect of tau
is tested against the alternative of a treatment effect smaller than tau.
In particular, alternative="less" is equivalent to: (i) alternative="greater",
(ii) y replaced by -y, and (iii) tau replaced by -tau.  See the note for
discussion of two-sided sensitivity analyses.
  }
}
\details{
For the given \eqn{\Gamma}, senfm computes the upper bound on the 1-sided
P-value testing the null hypothesis
of an additive treatment effect tau against the alternative hypothesis of
a treatment effect larger than tau.  By default, senfm tests the null hypothesis of
no treatment effect against the alternative of a positive treatment effect.
The P-value is an approximate P-value
based on a Normal approximation to the null distribution; see Rosenbaum (2007).

Matched sets of unequal size are weighted using weights that would
be efficient in a randomization test under a simple model with
additive set and treatment effects and errors with constant variance;
see Rosenbaum (2007).

The upper bound on the P-value is based on the separable approximation
described in Gastwirth, Krieger and Rosenbaum (2000); see also
Rosenbaum (2007).
}
\value{
\item{pval }{The upper bound on the 1-sided P-value.}
\item{deviate }{The deviate that was compared to the Normal
distribution to produce pval.}
\item{statistic }{The value of the M-statistic.}
\item{expectation }{The maximum expectation of the
M-statistic for the given \eqn{\Gamma}.}
\item{variance }{The maximum variance of the
M-statistic among treatment assignments
that achieve the maximum expectation.
Part of the separable approximation.}
}
\references{
Cox, D. R. (1977).  The role of signficance tests (with Discussion).
Scand. J. Statist. 4, 49-70.

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

Rosenbaum, P. R. and Small, D. S. (2016).  An adaptive Mantel-Haenszel
test for sensitivity analysis in observational studies.  Biometrics,
to appear.
}
\author{
Paul R. Rosenbaum.
}
\note{
The function senfm performs 1-sided tests.  One approach
to a 2-sided, \eqn{\alpha}-level test does both 1-sided tests
at level \eqn{\alpha/2}, and rejects the null hypothesis if either
1-sided
test rejects.  Equivalently, a bound on the two sided
P-value is the smaller of 1 and twice the smaller of the two 1-sided
P-values.  This approach views a 2-sided test as two 1-sided tests
with a Bonferroni correction; see Cox (1977, Section 4.2).  In all
cases, this approach is valid large sample test: a true
null hypothesis is falsely
rejected with probability at most \eqn{\alpha} if the bias in
treatment assignment is at most \eqn{\Gamma}; so, this procedure
is entirely safe to use.  For a randomization test, \eqn{\Gamma=1}, this
Bonferroni procedure is not typically conservative. For large \eqn{\Gamma},
this Bonferroni procedure tends to be somewhat conservative.

Related packages are sensitivitymv, sensitivitymw and sensitivity2x2xk.
sensitivitymv is for matched sets with one treated
subject and a variable number of controls.
sensitivitymw is for matched sets with one treated subject and a fixed
number of controls, including matched pairs.  For their special cases,
sensitivitymv and sensitivitymw contain additional features not available
in sensitivityfull.  sensitivitymw is faster and computes confidence intervals
and point estimates.  sensitivitymw also implements methods
from Rosenbaum (2014).  sensitivity2x2xk is for 2x2xk contingency tables,
treatment x outcome x covariates; see Rosenbaum and Small (2016).
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

#Randomization test of no effect, Huber scores:
senfm(y,treated1)

#Sensitivity analysis, Huber scores:
senfm(y,treated1,gamma=2)

#Randomization test of tau=1 vs tau>1
senfm(y,treated1,tau=1)

#Randomization test of tau=1 vs tau<1
senfm(y,treated1,tau=1,alternative="less")

#Same randomization test of tau=1 vs tau<1
senfm(-y,treated1,tau=-1)

#Sensitivity analysis testing tau=1 at gamma=2
senfm(y,treated1,tau=1,gamma=2,alternative="greater")
senfm(y,treated1,tau=1,gamma=2,alternative="less")

# For an additional example, install and load package sensitivitymv
# The following example is a match with variable controls.
# So this example has one treated subject per matched set.
# Both mscorev (in sensitivitymv) and mscoref (in sensitivityfull)
# reproduce parts of the example in Rosenbaum (2007, Section 4).
# data(tbmetaphase)
# senmv(tbmetaphase,gamma=2,trim=1)
# senfm(tbmetaphase,rep(TRUE,15),gamma=2,trim=1)
# senmv(tbmetaphase,gamma=2,trim=1,tau=0.94)
# senfm(tbmetaphase,rep(TRUE,15),gamma=2,trim=1,tau=.94)
# senmv(tbmetaphase,gamma=2,trim=1,tau=0.945)
# senfm(tbmetaphase,rep(TRUE,15),gamma=2,trim=1,tau=.945)
# mscoref(tbmetaphase,rep(TRUE,15),trim=1)
}

\keyword{ htest }
\keyword{ robust }
