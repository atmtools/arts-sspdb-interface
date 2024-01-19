# arts-sspdb-interface

This repository provides interfaces in MATLAB and Python intended for use with the ARTS Microwave Single Scattering Properties Database, which is available at zenodo: <https://doi.org/10.5281/zenodo.1175572>. The interfaces are also available at zenodo at <https://doi.org/10.5281/zenodo.1175572>.

The interfaces support browsing and importing of the data. There is also functionality for extracting data for usage in RTTOV.

## Interfaces

### MATLAB

MATLAB interface The interface should work with any relatively recent MATLAB version. Only installation required is adding the database folder to
MATLAB's search path. Storing data into files suitable as input to the ARTS forward model requires that the Atmlab3 package is installed and added to
MATLAB's search path.

### Python

Python interface The interface generally works in Python 3 environments, but Python 3 (  3.5) is suggested for full functionality. The interface
(hard-)requires the following Python packages to be installed:

- **netCDF4**
- **numpy**
- **os**

Certain features have additional requirements, namely (given in the form requirement: 
(module) feature"):

- **pyarts**
  specically the arts submodule: (assp, demo ssp4arts) Conversion to ARTS format and writing of ARTS SSP to XML file.

- **scipy**
  specically the interpolate submodule: (assp, sph) Grid interpolations.

- **SHTns** library
  (assp, sph)Conversion of azimuthally random orientation Z from Spherical Harmonics coeffcients to discrete angle grid values.

Install instructions for pyarts see <https://atmtools.github.io/arts-docs-master/installation.html>.
The SHTns library is available from <https://pypi.org/project/shtns/> or <https://bitbucket.org/nschaeff/shtns>
Important: Version 3.6.1 of shtns seems to be broken, use version 3.6.2 instead.

The netCDF4 is available with and installable through conda6, as well as available from PyPI7, or github8 and installable with pip9.
All other required Python packages listed above are part of standard Python distributions.

### RTTOV

RTTOV interface For a fully working RTTOV interface, a RTTOV installation (v11.3 or v12.1) as well as the Python interface and typhon are required.
The respective patch  le needs to be copied into the RTTOV installation's src/mw scatt coef/ subfolder and applied there by11
> patch -p1 < mw_scatt_coef_*.patch
After that, RTTOV needs to be rebuild in the usual way.

## Description of the database

The database contains microwave single scattering data, mainly of ice hydrometeors. Main applications lie in microwave remote sensing, both passive and active. Covered hydrometeors range from pristine ice crystals to large aggregates, graupel and hail, 34 types in total. Furthermore, 34 frequencies from 1 to 886 GHz, and 3 temperatures, 190, 230 and 270 K, are included. The main article orientation is currently totally random (i.e. each orientation is equally probable), but the database is designed to also handle azimuthally orientation, and data for this orientation case will be added in the future. Mie code was used for the two completely spherical habits, while the bulk of data were calculated using the discrete dipole approximation method (DDA).
