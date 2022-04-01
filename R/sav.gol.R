# ----------------------------------------------------------------------
## code added by Steffen Neumann 
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Savitzky-Golay Algorithm
# ----------------------------------------------------------------------
# T2 <- sav.gol(T, fl, forder=4, dorder=0);
#
# Polynomial filtering method of Savitzky and Golay
# See Numerical Recipes, 1992, Chapter 14.8, for details.
#
# T = vector of signals to be filtered
# (the derivative is calculated for each ROW)
# fl = filter length (for instance fl = 51..151)
# forder = filter order (2 = quadratic filter, 4= quartic)
# dorder = derivative order (0 = smoothing, 1 = first derivative, etc.)
#
sav.gol <- function(T, fl, forder=4, dorder=0)
{
    m <- length(T)
    dorder <- dorder + 1

    # -- calculate filter coefficients --
    fc <- (fl-1)/2 # index: window left and right
    X <- outer(-fc:fc, 0:forder, FUN="^") # polynomial terms and coefficients
    Y <- pinv(X); # pseudoinverse

    # -- filter via convolution and take care of the end points --
    #T2 <- convolve(T, rev(Y[dorder,]), type="o") # convolve(...) ## fails if length(T) is prime
    T2 <- conv(T, rev(Y[dorder,])) # convolve(...)
    T2 <- T2[(fc+1):(length(T2)-fc)]
}


#-----------------------------------------------------------------------
# *** PseudoInvers of a Matrix ***
# using singular value decomposition
#
pinv <- function (A)
{
    s <- svd(A)
    # D <- diag(s$d); Dinv <- diag(1/s$d)
    # U <- s$u; V <- s$v
    # A = U D V'
    # X = V Dinv U'
    s$v %*% diag(1/s$d) %*% t(s$u)
}
 

#-----------------------------------------------------------------------
# call c function "convolve2"
conv <- function(a, b) {
  .Call("convolve2", a, b, PACKAGE = "MassSpecWavelet")
}  

