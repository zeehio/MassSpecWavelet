# ----------------------------------------------------------------------
## code added by Steffen Neumann 
# ----------------------------------------------------------------------

conv <- function(a, b) {
  .Call("convolve2", a, b, PACKAGE = "MassSpecWavelet")
}  
