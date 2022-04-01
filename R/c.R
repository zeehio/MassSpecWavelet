# ----------------------------------------------------------------------
## code added by Steffen Neumann 
# ----------------------------------------------------------------------

conv <- function(a, b) {
  .Call("convolve2", as.double(a), as.double(b), PACKAGE = "MassSpecWavelet")
}  
