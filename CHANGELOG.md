Changelog
=========

All notable changes to this project will be documented in this file.

## [1.3] - 
### Added
- Added warning that PyQt4 and python2 support will be dropped soon
- Extract part function now has step size
- Load ASCII now has x-axis unit input
- Added support for Bruker ParaVision data
- Fitting window now shows the fit error (RMSD)
- Relaxation and diffusion fitting now show the equation in the fitting window
- Baseline correction with sine/cosine functions now possible
- A multiplot for contour plots has been added
- Projections in contour plots can now be zoomed by hovering over them and scrolling.
- Leaving x and y limits empty in save figure will reset those values
- Added an error message when loading data did not work
- Added a method to load DMFIT 1D spectrum data
- 2D Bruker EPR data can now be loaded
- Plot box can be turned off in save figure

### Changed
- WindowsRun.bat now also checks for pyw installer 
- Double list interface (e.g. combine workspace) how has buttons.
- New icon for export to ASCII
- Closing ssNake with no open workspace no longer prompts for confirmation
- Number of x and y ticks are now minimized

### Fixed
- Fixed an error on loading data in python2
- Fixed opening some utility windows
- Error on save figure in fitting window
- Deleting indices did not accept single integer
- Fixed a bug on undoing digital filter correction
- Fitting data with negative spectral widths lead to inverted spectra
- Fixed a bug in extract part
- Adding hypercomplex data sets had a bug
- Fixed a memory leak when working with large 1D data sets

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
