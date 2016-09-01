

#include "AlignmentConstant.h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCGenericObjectImpl.h>

// system includes <>
#include <iostream>
#include <iomanip>

using namespace lcio;
using namespace depfet;
using namespace std;

AlignmentConstant::AlignmentConstant() : IMPL::LCGenericObjectImpl( ALIGN_CONST_MAX_SIZE, 0 , ALIGN_CONST_MAX_SIZE ) {
  _typeName        = "Alignment constant";
  _dataDescription = "sensorID, xOff, yOff, zOff, alpha, beta, gamma + 13 spare fields";
  _isFixedSize     = true;


  for ( size_t iPos = 0; iPos < ALIGN_CONST_MAX_SIZE; ++iPos ) {
    setIntVal   ( iPos, 0   );
    setDoubleVal( iPos, 0.0 );
  }


}

AlignmentConstant::AlignmentConstant( int sensorID,
                                                double xOff,   double yOff,   double zOff,
                                                double alpha, double beta, double gamma,
                                                double xOffErr,double yOffErr,double zOffErr,
                                                double alphaErr, double betaErr, double gammaErr )  :

  IMPL::LCGenericObjectImpl(ALIGN_CONST_MAX_SIZE,0,ALIGN_CONST_MAX_SIZE) {

  _typeName        = "Alignment constant";
  _dataDescription = "sensorID,\n"
    "xOff,    yOff,    zOff,    alpha,    beta,    gamma\n"
    "xOffErr, yOffErr, zOffErr, alphaErr, betaErr, gammaErr";
  _isFixedSize     = true;

  for ( size_t iPos = 0; iPos < ALIGN_CONST_MAX_SIZE; ++iPos ) {
    setIntVal   ( iPos, 0   );
    setDoubleVal( iPos, 0.0 );
  }

  // set the sensor ID
  setIntVal( 0 , sensorID );

  // set the sensor shifts in mm
  setDoubleVal( 0, xOff );
  setDoubleVal( 1, yOff );
  setDoubleVal( 2, zOff );

  // set the sensor angles
  setDoubleVal( 3, alpha );
  setDoubleVal( 4, beta );
  setDoubleVal( 5, gamma );

  // set the sensor shifts in mm
  setDoubleVal( 6, xOffErr );
  setDoubleVal( 7, yOffErr );
  setDoubleVal( 8, zOffErr );

  // set the sensor angles
  setDoubleVal( 9, alphaErr );
  setDoubleVal( 10, betaErr );
  setDoubleVal( 11, gammaErr );

}

void AlignmentConstant::setSensorID( int id ) { setIntVal( 0 , id ) ; }

void AlignmentConstant::setXOffset( double off ) { setDoubleVal( 0, off ); }
void AlignmentConstant::setYOffset( double off ) { setDoubleVal( 1, off ); }
void AlignmentConstant::setZOffset( double off ) { setDoubleVal( 2, off ); }

void AlignmentConstant::setAlpha( double theta ) { setDoubleVal( 3, theta ); }
void AlignmentConstant::setBeta( double theta ) { setDoubleVal( 4, theta ); }
void AlignmentConstant::setGamma( double theta ) { setDoubleVal( 5, theta ); }

void AlignmentConstant::setXOffsetError( double err ) { setDoubleVal( 6, err ); }
void AlignmentConstant::setYOffsetError( double err ) { setDoubleVal( 7, err ); }
void AlignmentConstant::setZOffsetError( double err ) { setDoubleVal( 8, err ); }

void AlignmentConstant::setAlphaError( double err ) { setDoubleVal( 9, err ); }
void AlignmentConstant::setBetaError( double err ) { setDoubleVal( 10, err ); }
void AlignmentConstant::setGammaError( double err ) { setDoubleVal( 11, err ); }

int AlignmentConstant::getSensorID() const { return getIntVal( 0 ) ; }

double AlignmentConstant::getXOffset() const { return getDoubleVal( 0 ) ; }
double AlignmentConstant::getYOffset() const { return getDoubleVal( 1 ) ; }
double AlignmentConstant::getZOffset() const { return getDoubleVal( 2 ) ; }

double AlignmentConstant::getAlpha() const { return getDoubleVal( 3 ) ; }
double AlignmentConstant::getBeta() const { return getDoubleVal( 4 ) ; }
double AlignmentConstant::getGamma() const { return getDoubleVal( 5 ) ; }

double AlignmentConstant::getXOffsetError() const { return getDoubleVal( 6 ) ; }
double AlignmentConstant::getYOffsetError() const { return getDoubleVal( 7 ) ; }
double AlignmentConstant::getZOffsetError() const { return getDoubleVal( 8 ) ; }

double AlignmentConstant::getAlphaError() const { return getDoubleVal( 9 ) ; }
double AlignmentConstant::getBetaError() const { return getDoubleVal( 10 ) ; }
double AlignmentConstant::getGammaError() const { return getDoubleVal( 11 ) ; }

void AlignmentConstant::print(ostream & os ) const {

  const int maxFieldNo = 7;
  const int narrowWidth = 10;
  const int largeWidth = 14;
  const int lineWidth  = largeWidth * maxFieldNo;

  string dashline( lineWidth, '-' );

  os << dashline << endl
     << setw( lineWidth - narrowWidth ) << setiosflags(ios::left) << "Sensor ID " << resetiosflags(ios::left)
     << setw( narrowWidth ) << setiosflags(ios::right) << getSensorID() << resetiosflags(ios::right) << endl
     << dashline << endl
     << setw(largeWidth) << " X off [mm] "
     << setw(largeWidth) << " Y off [mm] "
     << setw(largeWidth) << " Z off [mm] "
     << setw(largeWidth) << " X angle "
     << setw(largeWidth) << " Y angle "
     << setw(largeWidth) << " Z angle "
     << endl
     << setw(largeWidth) << getXOffset()
     << setw(largeWidth) << getYOffset()
     << setw(largeWidth) << getZOffset()
     << setw(largeWidth) << getAlpha()
     << setw(largeWidth) << getBeta()
     << setw(largeWidth) << getGamma()
     << endl
     << setw(largeWidth) << getXOffsetError()
     << setw(largeWidth) << getYOffsetError()
     << setw(largeWidth) << getZOffsetError()
     << setw(largeWidth) << getAlphaError()
     << setw(largeWidth) << getBetaError()
     << setw(largeWidth) << getGammaError()
     << endl
     << dashline
     << endl;


  return;
}
