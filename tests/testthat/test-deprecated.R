## getRidgeValue(), i2u(), u2i(), mzInd2vRange(), mzV2indRange() and
## smoothDWT() are unused, both internally in MassSpecWavelet and by known
## downstream packages (checked against xcms). They are marked deprecated
## via .Deprecated() as a heads-up before removal; none of them are exported,
## so they are only reachable via `:::`.

test_that("getRidgeValue() warns that it is deprecated", {
    ridgeList <- list("1_2" = c(2L, 2L, 2L))
    attr(ridgeList, "scales") <- 0:2
    wCoefs <- matrix(1:15, nrow = 5, ncol = 3)
    expect_warning(
        MassSpecWavelet:::getRidgeValue(ridgeList, wCoefs),
        class = "deprecatedWarning"
    )
})

test_that("i2u() and u2i() warn that they are deprecated", {
    expect_warning(MassSpecWavelet:::i2u(100), class = "deprecatedWarning")
    expect_warning(MassSpecWavelet:::u2i(500000), class = "deprecatedWarning")
})

count_deprecated_warnings <- function(expr) {
    n <- 0
    withCallingHandlers(
        expr,
        deprecatedWarning = function(w) {
            n <<- n + 1
            invokeRestart("muffleWarning")
        }
    )
    n
}

test_that("mzInd2vRange() and mzV2indRange() warn exactly once, even for multi-element input", {
    # These call i2u()/u2i() internally in a loop (once per element); a naive
    # implementation would emit one deprecation warning per loop iteration
    # instead of once for the outer call.
    expect_equal(count_deprecated_warnings(MassSpecWavelet:::mzInd2vRange(1:10)), 1)
    expect_equal(count_deprecated_warnings(MassSpecWavelet:::mzV2indRange((1:10) * 100000)), 1)
})

test_that("smoothDWT() warns that it is deprecated even when it later errors", {
    skip_if(requireNamespace("waveslim", quietly = TRUE), "waveslim is installed, smoothDWT() would not error")
    expect_warning(
        expect_error(MassSpecWavelet:::smoothDWT(rnorm(100)), "waveslim"),
        class = "deprecatedWarning"
    )
})
