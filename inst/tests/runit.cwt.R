test01_cwt <- function() {
    skinny_peak <- c(
        9107, 3326, 9523, 3245, 3429, 9394, 1123, 935, 5128, 8576, 
        2711, 3427, 7294, 8109, 9288, 6997, 9756, 8034, 1317, 8866, 13877, 
        14854, 28296, 57101, 92209, 151797, 222386, 299402, 365045, 394255, 
        402680, 363996, 293985, 222989, 147007, 94947, 52924, 32438, 
        11511, 10836, 8046, 601, 889, 5917, 2690, 5381, 9901, 8494, 3349, 
        8283, 3410, 5935, 3332, 7041, 3284, 7478, 76, 3739, 2158, 5507
    )
    scales <- seq(1, 64, 2)
    wCoefs <- cwt(skinny_peak, scales = scales)
    checkEquals(length(skinny_peak), nrow(wCoefs), msg = "wCoefs should have as many rows as the signal length and not error")
}


test02_cwt <- function() {
    # All changes do not affect the result:
    skinny_peak <- c(
        rep(0, 1024),
        9107, 3326, 9523, 3245, 3429, 9394, 1123, 935, 5128, 8576, 
        2711, 3427, 7294, 8109, 9288, 6997, 9756, 8034, 1317, 8866, 13877, 
        14854, 28296, 57101, 92209, 151797, 222386, 299402, 365045, 394255, 
        402680, 363996, 293985, 222989, 147007, 94947, 52924, 32438, 
        11511, 10836, 8046, 601, 889, 5917, 2690, 5381, 9901, 8494, 3349, 
        8283, 3410, 5935, 3332, 7041, 3284, 7478, 76, 3739, 2158, 5507,
        rep(0, 1024)
    )
    
    scales <- c(1,2,4,8)
    wCoefs_classic <- MassSpecWavelet:::cwt_classic(skinny_peak, scales = scales, wavelet = "mexh")
    wCoefs_new <- cwt(skinny_peak, scales = scales, wavelet = "mexh")
    prep_wav <- prepareWavelets(
        length(skinny_peak),
        scales = scales, 
        wavelet = "mexh",
        wavelet_xlimit = 8,
        wavelet_length = 1024L,
        extendLengthScales = FALSE
    )
    wCoefs_new2 <- cwt(skinny_peak, prep_wav)
    checkEquals(wCoefs_classic, wCoefs_new, msg = "cwt classic is equal to new")
    checkEquals(wCoefs_new, wCoefs_new2, msg = "cwt new is equal to prepared wCoefs")
}