test01_peakDetection <- function() {
    skinny_peak <- c(
        9107, 3326, 9523, 3245, 3429, 9394, 1123, 935, 5128, 8576, 
        2711, 3427, 7294, 8109, 9288, 6997, 9756, 8034, 1317, 8866, 13877, 
        14854, 28296, 57101, 92209, 151797, 222386, 299402, 365045, 394255, 
        402680, 363996, 293985, 222989, 147007, 94947, 52924, 32438, 
        11511, 10836, 8046, 601, 889, 5917, 2690, 5381, 9901, 8494, 3349, 
        8283, 3410, 5935, 3332, 7041, 3284, 7478, 76, 3739, 2158, 5507
    )

    prep_wavelets <- prepareWavelets(length(skinny_peak))
    peak1 <- peakDetectionCWT(
        skinny_peak,
        prep_wavelets,
        excludeBoundariesSize = 0,
        exclude0scaleAmpThresh = TRUE
    )$majorPeakInfo$peakIndex
    checkEquals(1, length(peak1), msg = "Peak found in skinny_peak")

    prep_wavelets <- prepareWavelets(length(skinny_peak) + 512)
    peak2 <- peakDetectionCWT(
        c(rep(0, 256), skinny_peak, rep(0, 256)),
        prep_wavelets,
        excludeBoundariesSize = 0,
        exclude0scaleAmpThresh = TRUE
    )$majorPeakInfo$peakIndex
    checkEquals(1, length(peak2), msg = "Peak found in skinny_peak padded")
    checkEquals(unname(peak2), unname(peak1) + 256, msg = "Peaks do not match")
}

test02_peakDetection <- function() {
    
    # Test a wider peak
    # Values from round((dnorm(seq(-3, 3, length.out = 60))*100+runif(60))*10000)
    wider_peak <- c(5000, 12043, 15344, 12748, 20730, 20781, 24673, 36956, 44600, 
                    48596, 57698, 76937, 89422, 106482, 122977, 143989, 157769, 181563, 
                    206296, 226309, 251067, 283592, 307523, 324212, 341520, 368568, 
                    375716, 388428, 401694, 408352, 399415, 403964, 394144, 382952, 
                    368333, 341668, 330255, 301146, 276234, 254643, 231601, 211038, 
                    184239, 155817, 140996, 123284, 100121, 90280, 77303, 58708, 
                    52817, 44003, 36068, 24637, 20688, 14162, 14836, 16603, 8341, 
                    8307)
    prep_wavelets <- prepareWavelets(length(wider_peak))
    peak1 <- peakDetectionCWT(
        wider_peak,
        prep_wavelets,
        excludeBoundariesSize = 0,
        exclude0scaleAmpThresh = TRUE
    )
    
    checkEquals(30, unname(peak1$majorPeakInfo$peakIndex), msg = "Peak not found at index 30")
    
    prep_wavelets <- prepareWavelets(length(wider_peak) + 512)
    peak2 <- peakDetectionCWT(
        c(rep(0, 256), wider_peak, rep(0, 256)),
        prep_wavelets,
        excludeBoundariesSize = 0,
        exclude0scaleAmpThresh = TRUE
    )
    checkEquals(30 + 256, unname(peak2$majorPeakInfo$peakIndex), msg = "Peak not found at index 30+256")
    
    
}
