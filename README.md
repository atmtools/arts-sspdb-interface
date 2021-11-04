# arts-sspdb-interface
This repository provides interfaces in MATLAB and Python intended for use with the ARTS Microwave Single Scattering Properties Database, which is available at zenodo: https://doi.org/10.5281/zenodo.1175572. The interfaces are also available at zenodo at https://doi.org/10.5281/zenodo.1175572.

The interfaces support browsing and importing of the data. There is also functionality for extracting data for usage in RTTOV.

Further information on requirements and installation is provided in the separate **readme.pdf** file.

Description of the database:

The database contains microwave single scattering data, mainly of ice hydrometeors. Main applications lie in microwave remote sensing, both passive and active. Covered hydrometeors range from pristine ice crystals to large aggregates, graupel and hail, 34 types in total. Furthermore, 34 frequencies from 1 to 886 GHz, and 3 temperatures, 190, 230 and 270 K, are included. The main article orientation is currently totally random (i.e. each orientation is equally probable), but the database is designed to also handle azimuthally orientation, and data for this orientation case will be added in the future. Mie code was used for the two completely spherical habits, while the bulk of data were calculated using the discrete dipole approximation method (DDA).
