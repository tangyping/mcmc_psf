#include "Img.h"

using namespace CCfits;
using namespace std;

// The constructor just read in fits files and store the images into 
// signals and errs
Img::Img (vector<string>& f1, vector<string>& f2, vector<string>& f3, vector<string>& f4, VecDoub &wav) {

  nf = f1.size();

  fn1 = f1; // array of signal file names
  fn2 = f2; // array of noise file names
  fn3 = f3; // array of kernals
  fn4 = f4; // l & b maps
  wave = wav;

  //// start to read fits files ////
  /// get dimensionality of images
  unique_ptr<FITS> pIn(new FITS(fn1[0],Read,true));
  PHDU& data = pIn->pHDU();

  long ax1(data.axis(0));
  long ax2(data.axis(1));

  nx = ax2;
  ny = ax1;

  /////////////// read in dimensionality of kernels //////

  nx_ker = VecInt(nf);
  ny_ker = VecInt(nf);

  for(int i=0; i<nf; i++){
    unique_ptr<FITS> pIn_k(new FITS(fn3[i],Read,true));
    PHDU& data_k = pIn_k->pHDU();

    nx_ker[i] = long(data_k.axis(1));
    ny_ker[i] = long(data_k.axis(0));
  }

  /////////////// read in images //////////////////
  signals = new MatDoub[nf];
  errs = new MatDoub[nf];
  kers = new MatDoub[nf];
  lbs = new MatDoub[2];

  //// There should be nf signal/err/ker files
  for(int i=0; i<nf; i++){
    signals[i] = MatDoub(nx, ny);
    errs[i] = MatDoub(nx, ny);
    kers[i] = MatDoub(nx_ker[i], ny_ker[i]);
  }

  lbs[0] = MatDoub(nx,ny);
  lbs[1] = MatDoub(nx,ny);

  for(int i=0; i<nf; i++){

    unique_ptr<FITS> pIn1(new FITS(fn1[i],Read,true)); //img
    PHDU& data1 = pIn1->pHDU();

    unique_ptr<FITS> pIn2(new FITS(fn2[i],Read,true)); //err
    PHDU& data2 = pIn2->pHDU();

    unique_ptr<FITS> pIn3(new FITS(fn3[i],Read,true)); //ker
    PHDU& data3 = pIn3->pHDU();

    valarray<double>  cont1;
    valarray<double>  cont2;
    valarray<double>  cont3;

    data1.read(cont1);
    data2.read(cont2);
    data3.read(cont3);

    for (long j=0; j < nx; j++){
      for (long k=0; k < ny; k++){    
        signals[i][j][k] = cont1[j*ny+k];
        errs[i][j][k] = cont2[j*ny+k];
      }
    }

    for (long j=0; j < nx_ker[i]; j++){
      for (long k=0; k < ny_ker[i]; k++){    
        kers[i][j][k] = cont3[j*ny_ker[i]+k];
      }
    }
  }

  /////////////// read in l and bs ///////////////

  for(int i=0; i<2; i++){

    unique_ptr<FITS> pIn_lb(new FITS(fn4[i],Read,true));
    PHDU& data4 = pIn_lb->pHDU();
    valarray<double>  cont4;

    data4.read(cont4);

    for (long j=0; j < nx; j++){
      for (long k=0; k < ny; k++)
	lbs[i][j][k] = cont4[j*ny+k];	  	
    }
  }

}
