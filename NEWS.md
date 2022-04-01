# MassSpecWavelet 1.6.1 (2022-04-01)

- Change Maintainer to Sergio Oller
- Add NEWS.md file
- Fix warning (error on R>=4.2) when `cwt()` has a matrix in the `wavelet` argument
- Fix warning due to partial matching of arguments inside cwt()
- Move waveslim from Depends to Suggests
- Make all calls to recommended packages qualified (e.g. `sd` -> `stats::sd`)
- Remove empty sections in `man/` files. 
- Register C routines

