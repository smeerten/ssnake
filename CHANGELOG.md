Changelog
=========

All notable changes to this project will be documented in this file.

## [1.2] - 2019-12-06
### Added
- New plot type: 2D Colour Plot
- Added experimental support for Bruker ParaVision imaging data
- File browser tree: added 'open file external' and 'open browser here' option
- Reference Manual: added 'Roll' and Diffusion fitting
- Reference Manual: more on fitting 'AddData' and 'ConnectPar'
- Added many docstring to the source code

### Changed
- Minimum required Matplotlib version is now 1.5.0

### Fixed
- Temperature Calibration Tool: crash on start
- Baseline correction in macro did not work
- Fixed cos2 apodization (was not squared)

## [1.1] - 2019-05-03
### Added
- Tooltips have been added
- A fixed startup directory for ssNake can be set from preferences
- Spectra can be aligned based on the maximum position in a given region
- Added autophasing per trace
- The traces in a contour plot can now also display the values at the diagonal
- Added the option for color ranges in array and stacked plots
- The Czjzek distribution can now be stored in a file
- Added a fitting routine for mixed CSA and quadrupole data
- Added a function to roll data
- Added volume integration from the contour plot
- Ctrl and Shift can be used as modifiers to change the step size for various buttons
- Added support of MestreC data
- Linear prediction using LPSVD

### Changed
- In peak deconvolution, the resonances outside of the spectral windows are now dropped

### Fixed
- Support for Numpy 1.16

## [1.0] - 2019-01-21
### Added
- First official version of ssNake
