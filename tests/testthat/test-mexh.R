## mexh() is a trivial closed-form function, so the point of testing it isn't
## hunting for bugs -- it's pinning the exact wavelet formula as a regression
## guard, since cwt() defaults to wavelet = "mexh" and no existing test
## checks this formula directly (only end-to-end CWT behavior).

test_that("mexh() matches its documented closed form at known points", {
    expect_equal(mexh(0), 2 / sqrt(3) * pi^(-0.25))
    expect_equal(mexh(1), 0) # the (1 - x^2) root
    expect_equal(mexh(-1), 0)
    # A generic (non-root, non-origin) reference value, precise to 15 digits,
    # so a subtly wrong formula (e.g. a wrong exponent denominator) that
    # happens to agree at x = 0 and x = +-1 still gets caught.
    expect_equal(mexh(2), -0.352139052257134, tolerance = 1e-14)
})

test_that("mexh() is an even (symmetric) function", {
    x <- c(-3, -1.5, -0.5, 0.5, 1.5, 3)
    expect_equal(mexh(x), mexh(-x))
})

test_that("mexh() is positive for |x| < 1 and negative for |x| > 1", {
    expect_true(all(mexh(c(-0.9, -0.5, 0, 0.5, 0.9)) > 0))
    expect_true(all(mexh(c(-3, -2, -1.1, 1.1, 2, 3)) < 0))
})

test_that("mexh() decays to (near) zero far from the origin", {
    expect_lt(abs(mexh(20)), 1e-10)
})

test_that("mexh() is vectorized and preserves length, including the empty case", {
    expect_length(mexh(1:10), 10)
    expect_length(mexh(numeric(0)), 0)
})
