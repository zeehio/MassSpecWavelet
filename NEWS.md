# MassSpecWavelet development (2022-07-08)

- Drop sav.gol alias from the documentation. Related to #2

# MassSpecWavelet 1.61.3 (2022-04-05)

- Add signal to Suggests. Replace sav.gol() with signal::sgolayfilt(). Closes #2
- Minor documentation improvements

# MassSpecWavelet 1.61.2 (2022-04-04)

- Remove xcms and caTools. Those are great packages but we don't use them directly in MassSpecWavelet.
- Replace Sweave with RMarkdown vignette and update styles with BiocStyle

# MassSpecWavelet 1.61.1 (2022-04-01)

- Change Maintainer to Sergio Oller
- Add NEWS.md file
- Fix warning (error on R>=4.2) when `cwt()` has a matrix in the `wavelet` argument
- Fix warning due to partial matching of arguments inside cwt()
- Move waveslim from Depends to Suggests
- Make all calls to recommended packages qualified (e.g. `sd` -> `stats::sd`)
- Remove empty sections in `man/` files. 
- Register C routines

