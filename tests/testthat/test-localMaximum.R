test_that("localMaximum works on a short vector", {
    expect_equal(localMaximum(c(1, 2, 3, 4, 2, 1), winSize = 5), c(0, 0, 0, 1, 0, 0))
})

test_that("faster and classic algorithms give identical results", {
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
                x <- round(10 * runif(xlength), 1) * 10
                options(MassSpecWavelet.localMaximum.algorithm = "classic")
                localmax_classic <- localMaximum(x, winSize = winSize)
                options(MassSpecWavelet.localMaximum.algorithm = "faster")
                localmax_faster <- localMaximum(x, winSize = winSize)
                expect_equal(localmax_faster, localmax_classic)
            }
        }
    }
})

test_that("new algorithm does not miss any local maxima found by the faster algorithm", {
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
                x <- round(10 * runif(xlength), 1) * 10
                options(MassSpecWavelet.localMaximum.algorithm = "faster")
                localmax_faster <- localMaximum(x, winSize = winSize)
                options(MassSpecWavelet.localMaximum.algorithm = "new")
                localmax_new <- localMaximum(x, winSize = winSize)
                localmax_faster <- which(localmax_faster > 0)
                localmax_new <- which(localmax_new > 0)
                missing <- setdiff(localmax_faster, localmax_new)
                # Remove borders because classic&faster have false positives in them
                missing <- missing[-which(missing < winSize / 2 | missing > (xlength - winSize / 2))]
                # Remove plateaus because we have another criteria
                missing <- setdiff(missing, missing[x[missing] == x[missing + 1L]])
                expect_length(missing, 0)
            }
        }
    }
})
