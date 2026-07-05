## getLocalMaximumCWT()'s relative/exclude0scaleAmpThresh amplitude
## thresholding has been the source of several past bugfixes (see NEWS.md),
## but was only ever exercised indirectly through getRidge()/peakDetectionCWT()
## tests using their own defaults. These tests build a small, hand-crafted
## wCoefs matrix with a scale-0 (raw signal) column that has a much larger
## baseline than the wavelet-scale columns, which is exactly the scenario
## exclude0scaleAmpThresh is meant to fix.
buildWCoefsFixture <- function() {
    col0 <- c(1, 1, 1, 1, 1, 1, 1, 100, 1, 1, 1, 1, 1, 1, 1) # raw signal: big spike at row 8
    col1 <- c(1, 1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1) # scale=1: peaks 5 (row 5) and 2 (row 12)
    col4 <- c(1, 1, 1, 1, 1, 1, 1, 8, 1, 1, 1, 1, 1, 1, 1) # scale=4: peak 8 (row 8)
    wCoefs <- cbind(col0, col1, col4)
    colnames(wCoefs) <- c("0", "1", "4")
    wCoefs
}

test_that("getLocalMaximumCWT() finds the expected local maxima with the default (no threshold)", {
    wCoefs <- buildWCoefsFixture()
    localMax <- getLocalMaximumCWT(wCoefs)

    expect_equal(dim(localMax), dim(wCoefs))
    expect_equal(which(localMax[, "0"] > 0), 8)
    expect_equal(which(localMax[, "1"] > 0), c(5, 12))
    expect_equal(which(localMax[, "4"] > 0), 8)
})

test_that("an absolute amp.Th zeroes out local maxima below the threshold", {
    wCoefs <- buildWCoefsFixture()
    localMax <- getLocalMaximumCWT(wCoefs, amp.Th = 3)

    # The peak of height 2 at row 12 (scale 1) is now below the threshold...
    expect_equal(which(localMax[, "1"] > 0), 5)
    # ...while the taller peaks are unaffected.
    expect_equal(which(localMax[, "0"] > 0), 8)
    expect_equal(which(localMax[, "4"] > 0), 8)
})

test_that("a relative amp.Th without exclude0scaleAmpThresh can be dominated by the scale-0 baseline", {
    wCoefs <- buildWCoefsFixture()
    # max(wCoefs) == 100 (from the scale-0 baseline spike), threshold == 50.
    localMax <- getLocalMaximumCWT(wCoefs, amp.Th = 0.5, isAmpThreshRelative = TRUE, exclude0scaleAmpThresh = FALSE)

    # Real wavelet-scale peaks (5 and 8) are both wiped out by the inflated threshold.
    expect_equal(which(localMax[, "1"] > 0), integer(0))
    expect_equal(which(localMax[, "4"] > 0), integer(0))
    expect_equal(which(localMax[, "0"] > 0), 8)
})

test_that("exclude0scaleAmpThresh excludes the scale-0 column when computing a relative amp.Th", {
    wCoefs <- buildWCoefsFixture()
    # max(wCoefs[, c("1","4")]) == 8, threshold == 4.
    localMax <- getLocalMaximumCWT(wCoefs, amp.Th = 0.5, isAmpThreshRelative = TRUE, exclude0scaleAmpThresh = TRUE)

    expect_equal(which(localMax[, "1"] > 0), 5) # 5 >= 4, but 2 < 4 is dropped
    expect_equal(which(localMax[, "4"] > 0), 8)
    expect_equal(which(localMax[, "0"] > 0), 8) # scale-0 itself is never thresholded away
})

test_that("exclude0scaleAmpThresh has no effect when there is no scale-0 column", {
    wCoefs <- buildWCoefsFixture()[, c("1", "4")]
    withExclude <- getLocalMaximumCWT(wCoefs, amp.Th = 0.5, isAmpThreshRelative = TRUE, exclude0scaleAmpThresh = TRUE)
    withoutExclude <- getLocalMaximumCWT(wCoefs, amp.Th = 0.5, isAmpThreshRelative = TRUE, exclude0scaleAmpThresh = FALSE)
    expect_equal(withExclude, withoutExclude)
})

test_that("minWinSize overrides the scale-derived window size when it would be smaller", {
    wCoefs <- buildWCoefsFixture()
    # scale 1 -> winSize = 1*2+1 = 3, but minWinSize = 9 should be used instead.
    localMax <- getLocalMaximumCWT(wCoefs, minWinSize = 9)
    expect_equal(localMax[, "1"], localMaximum(wCoefs[, "1"], winSize = 9))
})

test_that("getLocalMaximumCWT() preserves column and row names", {
    wCoefs <- buildWCoefsFixture()
    rownames(wCoefs) <- paste0("r", seq_len(nrow(wCoefs)))
    localMax <- getLocalMaximumCWT(wCoefs)
    expect_equal(colnames(localMax), colnames(wCoefs))
    expect_equal(rownames(localMax), rownames(wCoefs))
})
