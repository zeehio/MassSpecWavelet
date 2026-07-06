test_that("tuneInPeakInfo() refines peak center/scale via majorPeakInfo (happy path)", {
    wider_peak <- c(5000, 12043, 15344, 12748, 20730, 20781, 24673, 36956, 44600,
                     48596, 57698, 76937, 89422, 106482, 122977, 143989, 157769, 181563,
                     206296, 226309, 251067, 283592, 307523, 324212, 341520, 368568,
                     375716, 388428, 401694, 408352, 399415, 403964, 394144, 382952,
                     368333, 341668, 330255, 301146, 276234, 254643, 231601, 211038,
                     184239, 155817, 140996, 123284, 100121, 90280, 77303, 58708,
                     52817, 44003, 36068, 24637, 20688, 14162, 14836, 16603, 8341,
                     8307)
    ms <- c(rep(0, 2000), wider_peak, rep(0, 2000))
    peakInfo <- peakDetectionCWT(ms, SNR.Th = 3, excludeBoundariesSize = 0, exclude0scaleAmpThresh = TRUE)
    majorPeakInfo <- peakInfo$majorPeakInfo
    expect_equal(unname(majorPeakInfo$peakIndex), 2030)
    expect_equal(unname(majorPeakInfo$peakScale), 22)

    tuned <- tuneInPeakInfo(ms, majorPeakInfo)
    expect_equal(unname(tuned$peakIndex), 2030) # peakIndex itself is never refined
    expect_equal(unname(tuned$peakCenterIndex), 2031)
    expect_equal(unname(tuned$peakScale), 23)
    expect_length(tuned$unProcessedPeak, 0)
    # peakSNR is rescaled by the ratio of the refined to the original peakValue.
    expect_equal(unname(tuned$peakSNR), unname(majorPeakInfo$peakSNR * tuned$peakValue / majorPeakInfo$peakValue))
})

test_that("tuneInPeakInfo() accepts peakIndex/peakScale directly, without majorPeakInfo", {
    wider_peak <- c(5000, 12043, 15344, 12748, 20730, 20781, 24673, 36956, 44600,
                     48596, 57698, 76937, 89422, 106482, 122977, 143989, 157769, 181563,
                     206296, 226309, 251067, 283592, 307523, 324212, 341520, 368568,
                     375716, 388428, 401694, 408352, 399415, 403964, 394144, 382952,
                     368333, 341668, 330255, 301146, 276234, 254643, 231601, 211038,
                     184239, 155817, 140996, 123284, 100121, 90280, 77303, 58708,
                     52817, 44003, 36068, 24637, 20688, 14162, 14836, 16603, 8341,
                     8307)
    ms <- c(rep(0, 2000), wider_peak, rep(0, 2000))
    tuned <- tuneInPeakInfo(ms, peakIndex = c(a = 2030), peakScale = c(a = 22))

    expect_equal(unname(tuned$peakCenterIndex), 2031)
    expect_equal(unname(tuned$peakScale), 23)
    expect_null(tuned$peakSNR) # no original peakSNR to rescale
    expect_length(tuned$unProcessedPeak, 0)
})

test_that("tuneInPeakInfo() requires either majorPeakInfo or peakIndex+peakScale", {
    ms <- rep(0, 100)
    expect_error(tuneInPeakInfo(ms), "majorPeakInfo or peakIndex and peakScale")
    expect_error(tuneInPeakInfo(ms, peakIndex = c(a = 10)), "majorPeakInfo or peakIndex and peakScale")
})

test_that("tuneInPeakInfo() validates the shape of majorPeakInfo", {
    ms <- rep(0, 100)
    expect_error(tuneInPeakInfo(ms, majorPeakInfo = list(foo = 1)), "Format of majorPeakInfo is incorrect")
})

test_that("tuneInPeakInfo() clips the candidate scale range near the start of the spectrum but still refines the peak", {
    skinny_peak <- c(
        9107, 3326, 9523, 3245, 3429, 9394, 1123, 935, 5128, 8576,
        2711, 3427, 7294, 8109, 9288, 6997, 9756, 8034, 1317, 8866, 13877,
        14854, 28296, 57101, 92209, 151797, 222386, 299402, 365045, 394255,
        402680, 363996, 293985, 222989, 147007, 94947, 52924, 32438,
        11511, 10836, 8046, 601, 889, 5917, 2690, 5381, 9901, 8494, 3349,
        8283, 3410, 5935, 3332, 7041, 3284, 7478, 76, 3739, 2158, 5507
    )
    ms <- c(rep(0, 50), skinny_peak, rep(0, 2000))
    prep_wavelets <- prepareWavelets(length(ms))
    peakInfo <- peakDetectionCWT(ms, prep_wavelets, excludeBoundariesSize = 0, exclude0scaleAmpThresh = TRUE)
    majorPeakInfo <- peakInfo$majorPeakInfo
    expect_equal(unname(majorPeakInfo$peakIndex), 81) # close enough to the start to require clipping

    tuned <- tuneInPeakInfo(ms, majorPeakInfo)
    expect_length(tuned$unProcessedPeak, 0)
    expect_equal(unname(tuned$peakCenterIndex), 81)
    expect_equal(unname(tuned$peakScale), 5)
})

test_that("tuneInPeakInfo() clips the candidate scale range near the end of the spectrum but still refines the peak", {
    skinny_peak <- c(
        9107, 3326, 9523, 3245, 3429, 9394, 1123, 935, 5128, 8576,
        2711, 3427, 7294, 8109, 9288, 6997, 9756, 8034, 1317, 8866, 13877,
        14854, 28296, 57101, 92209, 151797, 222386, 299402, 365045, 394255,
        402680, 363996, 293985, 222989, 147007, 94947, 52924, 32438,
        11511, 10836, 8046, 601, 889, 5917, 2690, 5381, 9901, 8494, 3349,
        8283, 3410, 5935, 3332, 7041, 3284, 7478, 76, 3739, 2158, 5507
    )
    ms <- c(rep(0, 2000), skinny_peak, rep(0, 52))
    prep_wavelets <- prepareWavelets(length(ms))
    peakInfo <- peakDetectionCWT(ms, prep_wavelets, excludeBoundariesSize = 0, exclude0scaleAmpThresh = TRUE)
    majorPeakInfo <- peakInfo$majorPeakInfo
    expect_equal(unname(majorPeakInfo$peakIndex), 2031)

    tuned <- tuneInPeakInfo(ms, majorPeakInfo)
    expect_length(tuned$unProcessedPeak, 0)
    expect_equal(unname(tuned$peakCenterIndex), 2003)
    expect_equal(unname(tuned$peakScale), 4)
})

test_that("tuneInPeakInfo() leaves a peak unprocessed when clipping near a boundary removes almost all candidate scales", {
    wider_peak <- c(5000, 12043, 15344, 12748, 20730, 20781, 24673, 36956, 44600,
                     48596, 57698, 76937, 89422, 106482, 122977, 143989, 157769, 181563,
                     206296, 226309, 251067, 283592, 307523, 324212, 341520, 368568,
                     375716, 388428, 401694, 408352, 399415, 403964, 394144, 382952,
                     368333, 341668, 330255, 301146, 276234, 254643, 231601, 211038,
                     184239, 155817, 140996, 123284, 100121, 90280, 77303, 58708,
                     52817, 44003, 36068, 24637, 20688, 14162, 14836, 16603, 8341,
                     8307)
    ms <- c(rep(0, 100), wider_peak, rep(0, 2000))
    peakInfo <- peakDetectionCWT(ms, SNR.Th = 3, excludeBoundariesSize = 0, exclude0scaleAmpThresh = TRUE)
    majorPeakInfo <- peakInfo$majorPeakInfo
    peakName <- names(majorPeakInfo$peakIndex)
    expect_equal(unname(majorPeakInfo$peakIndex), 130)

    tuned <- tuneInPeakInfo(ms, majorPeakInfo)
    expect_equal(tuned$unProcessedPeak, peakName)
    # Unprocessed peaks keep their original values.
    expect_equal(tuned$peakScale, majorPeakInfo$peakScale)
    expect_equal(tuned$peakValue, majorPeakInfo$peakValue)
    expect_equal(tuned$peakCenterIndex, majorPeakInfo$peakCenterIndex)
})

test_that("tuneInPeakInfo() leaves a peak unprocessed instead of crashing when its local window has no ridges", {
    # Regression test: a flat/quiet neighborhood around the peak makes
    # getRidge() return an empty ridge list within tuneInPeakInfo()'s
    # local refinement window, which used to crash with
    # "Error in strsplit(ridgeName.i, "_") : non-character argument".
    wider_peak <- c(5000, 12043, 15344, 12748, 20730, 20781, 24673, 36956, 44600,
                     48596, 57698, 76937, 89422, 106482, 122977, 143989, 157769, 181563,
                     206296, 226309, 251067, 283592, 307523, 324212, 341520, 368568,
                     375716, 388428, 401694, 408352, 399415, 403964, 394144, 382952,
                     368333, 341668, 330255, 301146, 276234, 254643, 231601, 211038,
                     184239, 155817, 140996, 123284, 100121, 90280, 77303, 58708,
                     52817, 44003, 36068, 24637, 20688, 14162, 14836, 16603, 8341,
                     8307)
    ms <- c(rep(0, 2000), wider_peak, rep(0, 2000))
    majorPeakInfo <- list(
        peakIndex = c(a = 500), peakCenterIndex = c(a = 500), peakScale = c(a = 22),
        peakValue = c(a = 1234), peakSNR = c(a = 10)
    )
    tuned <- tuneInPeakInfo(ms, majorPeakInfo)
    expect_equal(tuned$unProcessedPeak, "a")
    expect_equal(unname(tuned$peakScale), 22)
    expect_equal(unname(tuned$peakValue), 1234)
    expect_equal(unname(tuned$peakCenterIndex), 500)
    expect_equal(unname(tuned$peakSNR), 10)
})

test_that("tuneInPeakInfo() uses NA as a placeholder value when refinement is skipped and there is no original peakValue", {
    # Regression test: with the direct peakIndex/peakScale entry point,
    # peakValue/peakSNR are NULL. An unprocessed peak used to silently drop
    # an element (misaligning the output vectors) and then crash with
    # "attempt to set an attribute on NULL".
    wider_peak <- c(5000, 12043, 15344, 12748, 20730, 20781, 24673, 36956, 44600,
                     48596, 57698, 76937, 89422, 106482, 122977, 143989, 157769, 181563,
                     206296, 226309, 251067, 283592, 307523, 324212, 341520, 368568,
                     375716, 388428, 401694, 408352, 399415, 403964, 394144, 382952,
                     368333, 341668, 330255, 301146, 276234, 254643, 231601, 211038,
                     184239, 155817, 140996, 123284, 100121, 90280, 77303, 58708,
                     52817, 44003, 36068, 24637, 20688, 14162, 14836, 16603, 8341,
                     8307)
    ms <- c(rep(0, 2000), wider_peak, rep(0, 2000))
    tuned <- tuneInPeakInfo(ms, peakIndex = c(a = 500), peakScale = c(a = 22))
    expect_equal(tuned$unProcessedPeak, "a")
    expect_true(is.na(tuned$peakValue))
    expect_null(tuned$peakSNR)
})

test_that("tuneInPeakInfo() returns a well-formed empty result for zero peaks", {
    # Regression test: peakDetectionCWT(flat_signal, tuneIn = TRUE) used to
    # crash with "missing value where TRUE/FALSE needed" because
    # `for (i in 1:length(peakIndex))` iterates i = 1, 0 when peakIndex has
    # length 0 (e.g. the well-formed empty majorPeakInfo of a flat signal).
    flat_signal <- rep(5, 2001)
    peakInfo <- peakDetectionCWT(flat_signal, scales = c(1, seq(2, 30, 2), seq(32, 64, 4)), exclude0scaleAmpThresh = TRUE)
    expect_length(peakInfo$majorPeakInfo$peakIndex, 0)

    tuned <- tuneInPeakInfo(flat_signal, peakInfo$majorPeakInfo)
    for (field in setdiff(names(tuned), "unProcessedPeak")) {
        expect_length(tuned[[field]], 0)
    }
    expect_length(tuned$unProcessedPeak, 0)

    # The full peakDetectionCWT(tuneIn = TRUE) pipeline must not crash either.
    expect_no_error(peakDetectionCWT(flat_signal, scales = c(1, seq(2, 30, 2), seq(32, 64, 4)), exclude0scaleAmpThresh = TRUE, tuneIn = TRUE))
})

test_that("maxScale caps the candidate scale range used for refinement", {
    wider_peak <- c(5000, 12043, 15344, 12748, 20730, 20781, 24673, 36956, 44600,
                     48596, 57698, 76937, 89422, 106482, 122977, 143989, 157769, 181563,
                     206296, 226309, 251067, 283592, 307523, 324212, 341520, 368568,
                     375716, 388428, 401694, 408352, 399415, 403964, 394144, 382952,
                     368333, 341668, 330255, 301146, 276234, 254643, 231601, 211038,
                     184239, 155817, 140996, 123284, 100121, 90280, 77303, 58708,
                     52817, 44003, 36068, 24637, 20688, 14162, 14836, 16603, 8341,
                     8307)
    ms <- c(rep(0, 2000), wider_peak, rep(0, 2000))
    majorPeakInfo <- list(
        peakIndex = c(a = 2030), peakCenterIndex = c(a = 2030), peakScale = c(a = 30),
        peakValue = c(a = 1e6), peakSNR = c(a = 50)
    )

    capped <- tuneInPeakInfo(ms, majorPeakInfo, maxScale = 32)
    uncapped <- tuneInPeakInfo(ms, majorPeakInfo, maxScale = 128)

    expect_length(capped$unProcessedPeak, 0)
    expect_length(uncapped$unProcessedPeak, 0)
    # Capping to maxScale = 32 restricts the candidate scales to seq(22, 32, 0.5)
    # instead of seq(26, 34, 0.5), giving a different refined scale.
    expect_equal(unname(capped$peakScale), 23)
    expect_equal(unname(uncapped$peakScale), 26)
})
