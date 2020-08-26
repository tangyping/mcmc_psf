# mcmc_psf

This package requires CFITSIO and CCfits, which are C & C++ libraries for reading and writing data files in FITS.

CFITSIO can be found here:
https://heasarc.gsfc.nasa.gov/fitsio/

CCFITS can be found here:
https://heasarc.gsfc.nasa.gov/fitsio/CCfits/

Once CFITSIO and CCFITS are installed, change LDFLAGS in Makefile to the location where the libraries are installed, type:
make
./main
to compile and execute the program.

To setup initial conditions for the sampler, please check "setup" sections in main.cpp and "initial conditions" in Model.cpp.

This package also uses a slightly modified version of nr3.h from Numerical Recipe, which was obtained from:
http://numerical.recipes/codefile.php?nr3

Use readbest_singT.py to convert p[123]_med.dat to fits file.
