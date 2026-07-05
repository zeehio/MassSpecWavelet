## extendLength()/extendNBase() are internal (not exported) helpers that pad a
## signal/matrix used by cwt() to extend a spectrum before the wavelet
## transform. They are correctness-critical: an off-by-one in the padding
## silently corrupts every CWT coefficient near the signal edges. Only the
## default combination (method = "reflection", direction = "right") is ever
## exercised indirectly through cwt(); the other 8 method/direction
## combinations had no coverage at all.
extendLength <- MassSpecWavelet:::extendLength
extendNBase <- MassSpecWavelet:::extendNBase

test_that("extendLength() pads to the right with each method", {
    x <- 1:5
    expect_equal(extendLength(x, addLength = 2, method = "reflection", direction = "right"), c(1, 2, 3, 4, 5, 5, 4))
    expect_equal(extendLength(x, addLength = 2, method = "open", direction = "right"), c(1, 2, 3, 4, 5, 5, 5))
    expect_equal(extendLength(x, addLength = 2, method = "circular", direction = "right"), c(1, 2, 3, 4, 5, 1, 2))
})

test_that("extendLength() pads to the left with each method", {
    x <- 1:5
    expect_equal(extendLength(x, addLength = 2, method = "reflection", direction = "left"), c(2, 1, 1, 2, 3, 4, 5))
    expect_equal(extendLength(x, addLength = 2, method = "open", direction = "left"), c(1, 1, 1, 2, 3, 4, 5))
    expect_equal(extendLength(x, addLength = 2, method = "circular", direction = "left"), c(4, 5, 1, 2, 3, 4, 5))
})

test_that("extendLength() pads both sides with each method", {
    x <- 1:5
    # direction = "both" adds `addLength` on *each* side (not `addLength` split
    # in half), so the result grows by 2 * addLength.
    expect_equal(extendLength(x, addLength = 2, method = "reflection", direction = "both"), c(2, 1, 1, 2, 3, 4, 5, 5, 4))
    expect_equal(extendLength(x, addLength = 2, method = "open", direction = "both"), c(1, 1, 1, 2, 3, 4, 5, 5, 5))
    expect_equal(extendLength(x, addLength = 2, method = "circular", direction = "both"), c(4, 5, 1, 2, 3, 4, 5, 1, 2))
})

test_that("extendLength() extends a matrix column-wise", {
    m <- matrix(1:15, nrow = 5, ncol = 3)
    extended <- extendLength(m, addLength = 2, method = "reflection", direction = "right")
    expect_equal(dim(extended), c(7, 3))
    expect_equal(extended[, 1], c(1, 2, 3, 4, 5, 5, 4))
    expect_equal(extended[, 2], c(6, 7, 8, 9, 10, 10, 9))
    expect_equal(extended[, 3], c(11, 12, 13, 14, 15, 15, 14))
})

test_that("extendLength() with addLength = 0 returns the input unchanged", {
    expect_equal(extendLength(1:5, addLength = 0), 1:5)
})

test_that("extendLength() requires addLength and validates method/direction", {
    expect_error(extendLength(1:5), "provide the length to be added")
    expect_error(extendLength(1:5, addLength = 2, method = "bogus"), "should be one of")
    expect_error(extendLength(1:5, addLength = 2, direction = "bogus"), "should be one of")
})

test_that("extendNBase() extends row count to the next multiple of base^nLevel", {
    expect_equal(extendNBase(1:5, nLevel = 1, base = 2), c(1, 2, 3, 4, 5, 5))
    expect_equal(extendNBase(1:5, nLevel = 2, base = 2), c(1, 2, 3, 4, 5, 5, 4, 3))
    expect_equal(extendNBase(1:5, nLevel = NULL, base = 2), c(1, 2, 3, 4, 5, 5, 4, 3))
})

test_that("extendNBase() leaves an already-aligned length untouched", {
    # Note: when no extension is needed, extendNBase() returns a 1-column
    # matrix rather than a plain vector (extendLength() is never called, so
    # the vector-vs-matrix normalization at its end never runs). This is
    # existing, if slightly inconsistent, behavior that this test documents.
    result <- extendNBase(1:8, nLevel = 2, base = 2)
    expect_true(is.matrix(result))
    expect_equal(as.vector(result), 1:8)
})
