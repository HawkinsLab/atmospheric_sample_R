## Elyse Pennington
## new baselining routine: Fourier transform to remove light brightening
## 1) perform a Fourier transform on each wavelength over time
## 2) remove low frequency values (caused by slow changes in the lamp brightness
##    over time)
## 3) perform a reverse Fourier transform to convert the frequency data for each
##    wavelength back to time data