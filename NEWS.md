# MassSpecWavelet 1.63.2 (2022-09-29)

- Let getLocalMaximumCWT() have a relative amp.Th. Related to #4.
- Fix bug on identifyMajorPeaks() where nearbyWinSize was forced to 150
  if nearbyPeak was set to TRUE. Related to #4.
- Added excludeBoundariesSize argument to identifyMajorPeaks(). Before,
  nearbyWinSize was used for two different but related criteria: the
  range for including peaks close to a large peak AND the range to
  exclude peaks close to the beginning and end of the signal. Now,
  we have two independent arguments for each setting.
  The current behaviour does not change, but it is now more 
  flexible. Related to #4.

# MassSpecWavelet 1.63.1 (2022-07-08)

- Drop sav.gol alias from the documentation. Related to #2
- Explain how the scales argument is defined in relation to the mother wavelet.
- Add sessionInfo() to the vignette

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

