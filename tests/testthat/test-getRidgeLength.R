test_that("getRidgeLength() trims a ridge at the point it drops below Th * max", {
    expect_equal(getRidgeLength(list(c(5, 4, 3, 2, 1), c(5, 3, 1))), c(3, 2))
})

test_that("getRidgeLength() respects a custom Th", {
    expect_equal(getRidgeLength(list(c(10, 9, 8, 1)), Th = 0.9), 1)
})

test_that("getRidgeLength() returns the full length when no point drops below the threshold", {
    expect_equal(getRidgeLength(list(c(5, 5, 5, 5)), Th = 0.5), 4)
})

test_that("getRidgeLength() returns 1 for a single-element ridge regardless of Th", {
    expect_equal(getRidgeLength(list(7), Th = 0.99), 1)
})

test_that("getRidgeLength() preserves the names of a named ridgeList", {
    result <- getRidgeLength(list(a = c(5, 4, 3, 2, 1), b = c(5, 3, 1)))
    expect_equal(result, c(a = 3, b = 2))
})
