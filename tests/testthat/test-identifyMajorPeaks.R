## Fixture used by most tests below: a synthetic spectrum with three peaks
## designed to exercise identifyMajorPeaks()'s three filtering rules directly
## (bypassing peakDetectionCWT()'s own defaults):
##  - P1 (index ~800): tall/wide gaussian -> long ridge, high SNR. Passes every
##    rule and is far from both signal boundaries.
##  - P3 (index ~850): weak/narrow gaussian near P1 -> short ridge (fails the
##    ridgeLength rule alone), but close enough to P1 to be rescued when
##    nearbyPeak = TRUE.
##  - P4 (index ~30): identical shape to P1, but placed near the left edge of
##    the signal -> passes ridgeLength/SNR but fails the boundary-exclusion
##    rule under the default excludeBoundariesSize.
buildIdentifyMajorPeaksFixture <- function() {
    set.seed(123)
    n <- 2000
    ms <- rep(0, n) + rnorm(n, sd = 0.1)

    addGauss <- function(ms, mu, sigma, amp) {
        idx <- max(1, round(mu - 6 * sigma)):min(length(ms), round(mu + 6 * sigma))
        ms[idx] <- ms[idx] + amp * exp(-0.5 * ((idx - mu) / sigma)^2)
        ms
    }
    ms <- addGauss(ms, 800, 10, 5000) # P1
    ms <- addGauss(ms, 850, 1.5, 150) # P3
    ms <- addGauss(ms, 30, 10, 5000) # P4

    scales <- c(1, seq(2, 30, 2), seq(32, 64, 4))
    wCoefs <- cwt(ms, scales = scales, wavelet = "mexh")
    wCoefs <- cbind(as.vector(ms), wCoefs)
    colnames(wCoefs) <- c(0, scales)
    localMax <- getLocalMaximumCWT(wCoefs)
    colnames(localMax) <- colnames(wCoefs)
    ridgeList <- getRidge(localMax, gapTh = 3, skip = 2)

    list(ms = ms, wCoefs = wCoefs, ridgeList = ridgeList)
}

test_that("identifyMajorPeaks() returns a well-formed empty result when ridgeList has no ridges", {
    fx <- buildIdentifyMajorPeaksFixture()
    info <- identifyMajorPeaks(fx$ms, list(), fx$wCoefs)

    expect_named(
        info,
        c(
            "peakIndex", "peakValue", "peakCenterIndex", "peakSNR", "peakScale",
            "potentialPeakIndex", "allPeakIndex", "peakRidgeLengthScale",
            "peakNoise", "selInd"
        )
    )
    for (field in setdiff(names(info), "selInd")) {
        expect_length(info[[field]], 0)
    }
    expect_length(info$selInd$selInd1, 0)
    expect_length(info$selInd$selInd2, 0)
    expect_length(info$selInd$selInd3, 0)
})

test_that("ridgeLength filters short ridges from peakIndex, but they remain in allPeakIndex", {
    fx <- buildIdentifyMajorPeaksFixture()
    info <- identifyMajorPeaks(fx$ms, fx$ridgeList, fx$wCoefs, SNR.Th = 3, ridgeLength = 32, nearbyPeak = FALSE)

    # P1 (long ridge, high SNR, away from boundaries) is the only major peak.
    expect_equal(unname(info$peakIndex), 800)
    # P3 (short ridge) is filtered out of peakIndex...
    expect_false(850 %in% info$peakIndex)
    # ...even though its ridge was found and is reported in allPeakIndex.
    expect_true(850 %in% info$allPeakIndex)
})

test_that("nearbyPeak = TRUE rescues a short ridge located near a qualifying long ridge", {
    fx <- buildIdentifyMajorPeaksFixture()

    withoutNearby <- identifyMajorPeaks(
        fx$ms, fx$ridgeList, fx$wCoefs,
        SNR.Th = 3, ridgeLength = 32, nearbyPeak = FALSE, excludeBoundariesSize = 0
    )
    expect_false(850 %in% withoutNearby$peakIndex)

    withNearby <- identifyMajorPeaks(
        fx$ms, fx$ridgeList, fx$wCoefs,
        SNR.Th = 3, ridgeLength = 32, nearbyPeak = TRUE, nearbyWinSize = 150, excludeBoundariesSize = 0
    )
    expect_true(850 %in% withNearby$peakIndex)
})

test_that("SNR.Th filters low-SNR peaks, and SNR.Th = 0 disables the rule", {
    fx <- buildIdentifyMajorPeaksFixture()

    strict <- identifyMajorPeaks(fx$ms, fx$ridgeList, fx$wCoefs, SNR.Th = 1e6, ridgeLength = 32)
    expect_length(strict$peakIndex, 0)
    expect_false(any(strict$selInd$selInd2))

    disabled <- identifyMajorPeaks(fx$ms, fx$ridgeList, fx$wCoefs, SNR.Th = 0, ridgeLength = 32)
    expect_true(all(disabled$selInd$selInd2))
})

test_that("excludeBoundariesSize filters peaks near the signal edges, and excludeBoundariesSize = 0 disables the rule", {
    fx <- buildIdentifyMajorPeaksFixture()

    withBoundary <- identifyMajorPeaks(fx$ms, fx$ridgeList, fx$wCoefs, SNR.Th = 3, ridgeLength = 32, excludeBoundariesSize = 50)
    expect_false(30 %in% withBoundary$peakIndex)

    withoutBoundary <- identifyMajorPeaks(fx$ms, fx$ridgeList, fx$wCoefs, SNR.Th = 3, ridgeLength = 32, excludeBoundariesSize = 0)
    expect_true(30 %in% withoutBoundary$peakIndex)
})

test_that("an invalid SNR.method raises a descriptive error", {
    fx <- buildIdentifyMajorPeaksFixture()
    expect_error(
        identifyMajorPeaks(fx$ms, fx$ridgeList, fx$wCoefs, SNR.method = "bogus"),
        "Invalid SNR.method"
    )
})

test_that("a named 'fixed' minNoiseLevel is used as-is instead of being scaled by max(wCoefs)", {
    fx <- buildIdentifyMajorPeaksFixture()
    info <- identifyMajorPeaks(fx$ms, fx$ridgeList, fx$wCoefs, minNoiseLevel = c(fixed = 500))

    withScale <- info$peakValue != 0
    expect_true(all(info$peakNoise[withScale] >= 500))
})
