stopifnot(require(RUnit, quietly = TRUE))
stopifnot(require(MassSpecWavelet, quietly = TRUE))
## Define tests
testSuite <- defineTestSuite(name = "MassSpecWavelet Unit Tests",
                             dirs = system.file("tests", package = "MassSpecWavelet"),
                             testFuncRegexp = "^[Tt]est+",
                             rngKind = "Mersenne-Twister",
                             rngNormalKind = "Inversion"
)
tests <- runTestSuite(testSuite) # Run tests
printTextProtocol(tests) # Print results
# Return success or failure to R CMD CHECK
if (getErrors(tests)$nFail > 0) stop("TEST FAILED!")
if (getErrors(tests)$nErr > 0) stop("TEST HAD ERRORS!")
if (getErrors(tests)$nTestFunc < 1) stop("NO TEST FUNCTIONS RUN!")
