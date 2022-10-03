
test01_localMaximum <- function() {
    checkEquals(localMaximum(c(1,2,3,4,2,1), winSize = 5), c(0,0,0,1,0,0), msg = "Check localMaximum on a short vector")
}