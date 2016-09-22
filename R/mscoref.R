mscoref <-
function (ymat, treated1, inner = 0, trim = 3, qu = 0.5)
    {
        stopifnot(is.logical(treated1))
        stopifnot(inner>=0)
        stopifnot(trim>=inner)
        stopifnot((qu>0)&(qu<1))
        ymat <- as.matrix(ymat)
        n <- dim(ymat)[1]
        m <- dim(ymat)[2]
        out <- matrix(NA, n, m)
        one <- rep(1, m - 1)
        difs <- array(NA, c(n, m, m - 1))
        for (j in 1:m) {
            difs[, j, ] <- outer(as.vector(unlist(ymat[, j])), one,
                                 "*") - ymat[, -j]
        }
        ms <- as.vector(difs)
        if ((trim < Inf) | (inner > 0)) {
            hqu <- as.numeric(stats::quantile(abs(ms), qu, na.rm = TRUE))
            if (hqu > 0) {
                ms <- ms/hqu
                if ((trim < Inf) & (inner < trim)) {
                    ab <- pmin(1, pmax(0, (abs(ms) - inner))/(trim -
                                                                  inner))
                }
                else if ((trim < Inf) & (inner == trim)) {
                    ab <- 1 * (abs(ms) > inner)
                }
                else {
                    ab <- pmax(0, abs(ms) - inner)
                }
                ms <- sign(ms) * ab
            }
            else {
                "Error: Scale factor is zero.  Increase lambda."
            }
        }
        ms <- array(ms, c(n, m, m - 1))
        ms <- apply(ms, c(1, 2), sum, na.rm = TRUE)
        ms[is.na(ymat)] <- NA
        colnames(ms) <- colnames(ymat)
        ni <- apply(!is.na(ymat), 1, sum)
        use <- (ni >= 2) & (!is.na(ms[, 1]))
        ms <- ms[use, ]
        ni <- ni[use]
        ms <- ms/outer(ni, rep(1, m), "*")
        ms[!treated1,]<-(-ms[!treated1,])
        ms
    }
