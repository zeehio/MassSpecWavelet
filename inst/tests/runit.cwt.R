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
