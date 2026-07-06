u2i <- function(u, t0 = 9.0E-8, a = 3.36301655051424E8, b = 0, U = 20000, delta = 4.0E-9) {
    .Deprecated(msg = paste(
        "u2i() appears to be unused, both internally in MassSpecWavelet and by",
        "known downstream packages, and is scheduled for removal in a future",
        "release. If you rely on it, please open an issue at",
        "https://github.com/zeehio/MassSpecWavelet/issues within the next year."
    ))
    i <- (((u / U - b) / a)^(1 / 2) + t0) / delta
    return(floor(i)) # just take the integer part
}
