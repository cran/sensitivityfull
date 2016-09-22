senfm <-
function (y, treated1, gamma = 1, inner = 0, trim = 3, lambda = 1/2,
              tau = 0, alternative="greater")
    {
        stopifnot(is.logical(treated1))
        stopifnot(inner>=0)
        stopifnot(trim>=inner)
        stopifnot((lambda>0)&(lambda<1))
        stopifnot(gamma>=1)
        stopifnot((alternative=="greater")|(alternative=="less"))
        if (is.vector(y)) {
            y <- y[!is.na(y)]
            treat <- y/2
            cont <- (-y/2)
            y <- cbind(treat, cont)
        }
        stopifnot(min(apply(!is.na(y), 1, sum)) >= 2)
        stopifnot(length(treated1)==(dim(y)[1]))
        if (alternative=="less"){
          y<-(-y)
          tau<-(-tau)
        }
        if (!(tau == 0)){
            if (sum(treated1)>0) y[treated1, 1] <- y[treated1, 1] - tau
            if (sum((1-treated1))>0) y[!treated1,-1] <- y[!treated1,-1]-tau
		}
        ms <- mscoref(y, treated1, inner = inner, trim = trim, qu = lambda)
        separable1f(ms, gamma = gamma)
    }
