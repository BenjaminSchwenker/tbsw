#include "DEPFET.h"
#include "MatrixDecoder.h"

// system includes <>
#include <iostream>

using namespace std;
using namespace lcio;
using namespace depfet; 

MatrixDecoder::MatrixDecoder(int xNoOfPixel, int yNoOfPixel)  {

  if ( xNoOfPixel <= 0 ) {
    streamlog_out (ERROR3) << "MatrixDecoder: xNoOfPixel has to be positive" << std::endl;
    exit(-1); 
  }
  if ( yNoOfPixel <= 0 ) {
    streamlog_out (ERROR3) << "MatrixDecoder: yNoOfPixel has to be positive" << std::endl;
    exit(-1); 
  }

  _xNoOfPixel = xNoOfPixel;
  _yNoOfPixel = yNoOfPixel;
  _xMin = 0;
  _yMin = 0;

}

MatrixDecoder::MatrixDecoder(int xNoOfPixel, int yNoOfPixel, int xMin, int yMin) {
  
   if ( xNoOfPixel <= 0 ) {
    streamlog_out (ERROR3) << "MatrixDecoder: xNoOfPixel has to be positive" << std::endl;
    exit(-1); 
  }
  if ( yNoOfPixel <= 0 ) {
    streamlog_out (ERROR3) << "MatrixDecoder: yNoOfPixel has to be positive" << std::endl;
    exit(-1); 
  }
  
  _xNoOfPixel = xNoOfPixel;
  _yNoOfPixel = yNoOfPixel;
  _xMin = xMin;
  _yMin = yMin;
}


int MatrixDecoder::getIndexFromXY(int x, int y) const {
  
  int xCor = x - _xMin;
  int yCor = y - _yMin;
  return xCor + yCor * _xNoOfPixel;

}

void MatrixDecoder::getXYFromIndex(int index, int& x, int& y) const {
  
  y = getYFromIndex(index);
  x = getXFromIndex(index);

}

int MatrixDecoder::getXFromIndex(int index) const {
  
  return ( index % _xNoOfPixel ) + _xMin;

}

int MatrixDecoder::getYFromIndex(int index) const {
  
  return ( index / _xNoOfPixel ) + _yMin;

}




