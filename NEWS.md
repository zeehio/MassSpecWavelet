# MassSpecWavelet 1.65.1 (2023-04-07)

- Fix .Call() for R-4.3. Thanks to Steffen Neumann. Closes #5

# MassSpecWavelet 1.64.1 (2023-01-30)

- Fix undefined variable in MassSpecWavelet.Rmd

# MassSpecWavelet 1.63.6 (2022-10-15)

- Fix regression in `cwt()` when scales has length 1.

# MassSpecWavelet 1.63.5 (2022-10-11)

- `getRidge()` supports a `scaleToWinSize` parameter. This argument controls how
  scales get mapped to window sizes. These windows are used to track the local
  maxima into ridges. MassSpecWavelet had a criteria of `winsize <- 2*scale+1`,
  while xcms modified it to `winsize <- floor(scale/2)`. This new argument enables
  xcms maintainers to call MassSpecWavelet's getRidge (if they want to) using their
  criteria, while it still lets us preserve backwards compatibility in our results.
  See `?getRidge` for further details.
  
- The `getLocalMaximumCWT()` `is_amp_thres_relative` parameter is now `isAmpThreshRelative`,
  for consistency with other parameter capitalization in the package. Since it
  was introduced 10 days ago, I don't think there will be more than one user using
  it.

- `getLocalMaximumCWT()` and `peakDetectionCWT` have a `exclude0scaleAmpThresh`
  parameter. When computing the relative `amp.Th`, if this parameter is set
  to `TRUE`, the `amp.Th` will exclude the zero-th scale from the
  `max(wCoefs)`. The zero-th scale corresponds to the original signal, that may
  have a much larger baseline than the wavelet coefficients and can distort the
  threshold calculation. The default value is `FALSE` to preserve backwards
  compatibility.

- `peakDetectionCWT` lets the user pass custom arguments to `getRidge()`.


# MassSpecWavelet 1.63.4 (2022-10-10)

- The improvements in `localMaxima()` and `cwt()` provide significant speed-ups
  to `peakDetectionCWT()` as well as better scalability.

- A `prepareWavelets()` function lets the user pre-compute the daughter wavelets
  for more efficient `cwt()` calculations when applied on multiple spectra. When
  used transforming 1000 spectra, of 2000 points long each, using 25 different
  scales, `cwt()` is twice as fast as in previous versions. Further improvements
  to avoid some memory allocations are still feasible in future versions.
  
- Through the `prepareWavelets()` function, we provide the `extendLengthScales`
  argument, that provides the same functionality than the `extendLengthMSW` argument
  in `xcms:::MSW.cwt()`.

- The `peakDetectionCWT()` function accepts a `prepared_wavelets` object in the
  scales argument for better efficiency.

# MassSpecWavelet 1.63.3 (2022-10-10)

- `localMaxima()` has a more efficient implementation of the algorithm, now being
  10x faster than before, while giving the same results.

- Experimentally, `localMaxima()` can use a new and different algorithm for
  detecting local maxima. See the new "Finding local maxima" vignette for
  further details.

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

