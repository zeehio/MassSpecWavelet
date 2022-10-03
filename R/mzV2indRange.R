"mzV2indRange" <-
    function(mzV, error = 0.003) {
        mzIndR <- NULL
        for (i in 1:length(mzV)) {
            from.i <- round(u2i(mzV[i] * (1 - error)))
            to.i <- round(u2i(mzV[i] * (1 + error)))
            mzIndR <- c(mzIndR, from.i:to.i)
        }
        mzIndR <- sort(unique(mzIndR))
        return(mzIndR)
    }
