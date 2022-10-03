"u2i" <-
    function(u, t0 = 9.0E-8, a = 3.36301655051424E8, b = 0, U = 20000, delta = 4.0E-9) {
        i <- (((u / U - b) / a)^(1 / 2) + t0) / delta
        return(floor(i)) # just take the integer part
    }
