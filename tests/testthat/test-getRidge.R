test_that("getRidge() never uses scientific notation for large m/z indices (#8)", {
    # Regression test for issue #8: ridge names are built by pasting numeric
    # m/z indices together (e.g. "1_653"). R's default numeric-to-character
    # coercion switches to scientific notation for large "round" doubles
    # (e.g. 200000 -> "2e+05"), and getRidge() internally promotes indices
    # from integer to double partway through its tracing loop. Without
    # consistent, representation-independent string keys, the same peak
    # index could end up looked up under two different strings, silently
    # breaking the ridge tracing (dropped/duplicated ridges).
    #
    # High resolution spectra can easily have m/z indices at or above 1e5,
    # so we build a spectrum with peaks at "round" indices, which is exactly
    # what used to trigger scientific notation.
    n <- 120000
    ms <- rep(0, n)
    centers <- c(100000, 110000)
    for (mu in centers) {
        idx <- max(1, mu - 30):min(n, mu + 30)
        ms[idx] <- ms[idx] + 1000 * exp(-0.5 * ((idx - mu) / 8)^2)
    }

    scales <- c(1, seq(2, 30, 2), seq(32, 64, 4))
    wCoefs <- cwt(ms, scales = scales, wavelet = "mexh")
    wCoefs <- cbind(as.vector(ms), wCoefs)
    colnames(wCoefs) <- c(0, scales)
    localMax <- getLocalMaximumCWT(wCoefs)
    colnames(localMax) <- colnames(wCoefs)

    ridgeList <- getRidge(localMax, gapTh = 3, skip = 2)

    ridgeName <- names(ridgeList)
    expect_gt(length(ridgeList), 0)
    expect_false(any(grepl("e[+-]", ridgeName, ignore.case = TRUE)))
    expect_true(all(grepl("^[0-9]+_[0-9]+$", ridgeName)))

    # The m/z index encoded in the ridge name must match the first element of
    # the ridge itself: if a lookup mismatch had split/duplicated a ridge
    # in getRidge(), this invariant would be the first thing to break.
    ridgeInfo <- matrix(as.numeric(unlist(strsplit(ridgeName, "_"))), nrow = 2)
    mzIndFromName <- ridgeInfo[2, ]
    mzIndFromRidge <- sapply(ridgeList, function(x) x[1])
    expect_equal(unname(mzIndFromRidge), unname(mzIndFromName))

    # Each synthetic peak should be found by exactly one ridge, reaching the
    # finest scale (i.e. not truncated/orphaned because of a name mismatch).
    for (mu in centers) {
        nearby <- which(abs(mzIndFromRidge - mu) <= 5)
        expect_length(nearby, 1)
        expect_gte(length(ridgeList[[nearby]]), length(scales) - 3)
    }
})
