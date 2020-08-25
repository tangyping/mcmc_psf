#include <stdio.h>
#include <iostream>
#include "nr3.h"
#include "Writer.h"


bool writeMatOut(const char* outFile, MatDoub &outMat)
{
  ofstream out(outFile);
  if(out.bad()){
    cerr << "vector_utilities::writeMatOut():";
    cerr << " Error opening " << outFile << endl;
    return 0;
  }

  //set precision to something reasonable for det values
  out.precision(14);

  //and the output
  for(int i=0;i<outMat.nrows();i++)
  for(int j=0;j<outMat.ncols();j++)
    {
      out << outMat[i][j] << "\n";
      if(out.bad()){
	cerr << "vector_utilities::writeMatOut():";
	cerr << " Error writing to " << outFile << endl;
	return 0;
      }
    }

  //flush the buffer just to be sure
  out.flush();

  //announce what you just did
  cout << "Wrote out " << outFile <<  " with nrows=" << outMat.nrows();
  cout << " and ncols=" << outMat.ncols() << endl;

  return 1;
}



bool writeVecOut(const char* outFile,  VecDoub outData)
{
  ofstream out(outFile);
  if(out.bad()){
    cerr << "vector_utilities::writeVecOut():";
    cerr << " Error opening " << outFile << endl;
    return 0;
  }

  //set precision to something reasonable for det values
  out.precision(14);

  //and the output
  for(size_t i=0;i<outData.size();i++) out << outData[i] << "\n";
  if(out.bad()){
    cerr << "vector_utilities::writeVecOut():";
    cerr << " Error writing to " << outFile << endl;
    return 0;
  }

  //flush the buffer just to be sure
  out.flush();

  //announce what you just did
  cout << "Wrote out " << outFile << endl;

  return 1;
}
