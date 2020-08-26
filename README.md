# mcmc_psf

This package requires CFITSIO and CCfits, which are libraries for C and C++ for reading and writing data files in FITS.

CFITSIO can be found here:
https://heasarc.gsfc.nasa.gov/fitsio/

CCFITS can be found here:
https://heasarc.gsfc.nasa.gov/fitsio/CCfits/

Once CFITSIO and CCFITS are installed, change LDFLAGS in Makefile to the location where the libraries are installed, just type:
make
./main
to compile and excute the program.

To setup initial conditions for the sampler, please check "setup" sections in main.cpp and "initial conditions" in Model.cpp.

This package also uses a slightly modified version of nr3.h from Numerical Receipe, which was obtained from:
http://numerical.recipes/codefile.php?nr3

