## findLocalMaxWinSize() has zero test coverage of its C implementation
## (find_local_maximum.c). Despite its roxygen comment saying `@export`, it is
## not actually listed in NAMESPACE (see R/findLocalMaxWinSize.R).
##
## It's unexported, but testthat's test environment is a clone of the package
## namespace (see testthat:::test_env()), so it resolves here by its bare
## name, same as exported functions.

test_that("findLocalMaxWinSize() matches its own documented example", {
    expect_equal(findLocalMaxWinSize(c(1, 2, 3, 2, 1)), c(0, 0, 5, 0, 0))
})

test_that("findLocalMaxWinSize() reports 0 everywhere for a monotonic signal", {
    expect_equal(findLocalMaxWinSize(c(1, 2, 3, 4, 5)), rep(0, 5))
    expect_equal(findLocalMaxWinSize(c(5, 4, 3, 2, 1)), rep(0, 5))
})

test_that("findLocalMaxWinSize() reports 0 everywhere for a fully flat signal", {
    # A signal with no rising/falling edge anywhere has no discernible peak,
    # even though every point is technically tied for the maximum.
    expect_equal(findLocalMaxWinSize(c(3, 3, 3, 3, 3)), rep(0, 5))
    expect_equal(findLocalMaxWinSize(c(3, 3, 3, 3)), rep(0, 4))
})

test_that("an equal-height neighboring peak does not block window growth", {
    # Two isolated peaks of the same height: neither is strictly greater than
    # the other, so each gets a window spanning the whole signal.
    expect_equal(findLocalMaxWinSize(c(0, 1, 0, 0, 0, 1, 0)), c(0, 7, 0, 0, 0, 7, 0))
    expect_equal(findLocalMaxWinSize(c(0, 1, 0, 1, 0)), c(0, 5, 0, 5, 0))
})

test_that("capWinSize caps the reported window size, and NA means unlimited", {
    x <- c(1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1) # single peak at position 6
    expect_equal(findLocalMaxWinSize(x, capWinSize = 3L)[6], 3L)
    expect_equal(findLocalMaxWinSize(x, capWinSize = 1L)[6], 1L)
    expect_equal(findLocalMaxWinSize(x, capWinSize = 0L)[6], 0L)
    expect_equal(findLocalMaxWinSize(x, capWinSize = NA)[6], 11L) # unlimited: the whole signal
})

test_that("findLocalMaxWinSize() reports the exact center of an odd-length plateau", {
    x <- c(0, 1, 2, 3, 3, 3, 2, 1, 0) # plateau of 3 points at positions 4:6
    result <- findLocalMaxWinSize(x)
    expect_equal(which(result > 0), 5) # the middle of the plateau
})

test_that("findLocalMaxWinSize() reports the first of the two centers of an even-length plateau", {
    x <- c(0, 1, 2, 3, 3, 2, 1, 0) # plateau of 2 points at positions 4:5
    result <- findLocalMaxWinSize(x)
    expect_equal(which(result > 0), 4) # the first of the two candidate centers
})

test_that("findLocalMaxWinSize() handles length 0, 1 and 2 inputs without crashing", {
    expect_equal(findLocalMaxWinSize(numeric(0)), integer(0))
    expect_equal(findLocalMaxWinSize(5), 0L)
    expect_equal(findLocalMaxWinSize(c(1, 2)), c(0L, 0L))
})
