
test01_localMaximum <- function() {
    checkEquals(localMaximum(c(1,2,3,4,2,1), winSize = 5), c(0,0,0,1,0,0), msg = "Check localMaximum on a short vector")
}


test02_checkCorrectness <- function() {
    set.seed(5413L)
    winSizes <- c(5, 31, 301)
    xlengths <- c(20, 200, 2000, 20000)
    simulations <- 1:20
    for (winSize in winSizes) {
        for (xlength in xlengths) {
            for (simulation in simulations) {
                x <- round(10*runif(xlength), 1)*10
                options(MassSpecWavelet.localMaximum.algorithm = "classic")
                localmax_classic <- which(localMaximum(x, winSize = winSize) > 0)
                options(MassSpecWavelet.localMaximum.algorithm = "new")
                localmax_new <- which(localMaximum(x, winSize = winSize) > 0)
                missing <- setdiff(localmax_classic, localmax_new)
                # Remove borders because classic has false positives in them
                missing <- missing[-which(missing < winSize/2 | missing > (xlength - winSize/2))]
                # Remove plateaus because we have another criteria
                missing <- setdiff(missing, missing[x[missing] == x[missing + 1L]])
                checkEquals(0, length(missing), msg = sprintf("new algorithm misses some peaks"))
            }
        }
    }
}