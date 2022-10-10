
test01_localMaximum <- function() {
    checkEquals(localMaximum(c(1,2,3,4,2,1), winSize = 5), c(0,0,0,1,0,0), msg = "Check localMaximum on a short vector")
}

test02_checkIdenticalResults_classic_vs_faster <- function() {
    on.exit({
        options(MassSpecWavelet.localMaximum.algorithm = NULL)
    })
    set.seed(5413L)
    winSizes <- c(5, 31, 301)
    xlengths <- c(20, 200, 2000, 20000)
    simulations <- 1:20
    for (winSize in winSizes) {
        for (xlength in xlengths) {
            for (simulation in simulations) {
                x <- round(10*runif(xlength), 1)*10
                options(MassSpecWavelet.localMaximum.algorithm = "classic")
                localmax_classic <- localMaximum(x, winSize = winSize)
                options(MassSpecWavelet.localMaximum.algorithm = "faster")
                localmax_faster <- localMaximum(x, winSize = winSize)
                checkEquals(
                    localmax_classic,
                    localmax_faster,
                    msg = "faster and classic do not give identical results"
                )
            }
        }
    }
}


test03_check_new_does_not_miss_maxima <- function() {
    on.exit({
        options(MassSpecWavelet.localMaximum.algorithm = NULL)
    })
    set.seed(5413L)
    winSizes <- c(5, 31, 301)
    xlengths <- c(20, 200, 2000, 20000)
    simulations <- 1:20
    for (winSize in winSizes) {
        for (xlength in xlengths) {
            for (simulation in simulations) {
                x <- round(10*runif(xlength), 1)*10
                options(MassSpecWavelet.localMaximum.algorithm = "faster")
                localmax_faster <- localMaximum(x, winSize = winSize)
                options(MassSpecWavelet.localMaximum.algorithm = "new")
                localmax_new <- localMaximum(x, winSize = winSize)
                localmax_faster <- which(localmax_faster > 0)
                localmax_new <- which(localmax_new > 0)
                missing <- setdiff(localmax_faster, localmax_new)
                # Remove borders because classic&faster have false positives in them
                missing <- missing[-which(missing < winSize/2 | missing > (xlength - winSize/2))]
                # Remove plateaus because we have another criteria
                missing <- setdiff(missing, missing[x[missing] == x[missing + 1L]])
                checkEquals(0, length(missing), msg = sprintf("new algorithm misses some peaks"))
            }
        }
    }
}