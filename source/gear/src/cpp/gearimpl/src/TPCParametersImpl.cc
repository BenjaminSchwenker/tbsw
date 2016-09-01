#include "gearimpl/TPCParametersImpl.h"
#include "gear/PadRowLayout2D.h"


namespace gear{
  
  TPCParametersImpl::TPCParametersImpl() : 
    _padLayout(0),
    _maxDriftLength (0.),
    _driftVelocity (0.),
    _readoutFrequency (0.) {
    
  }
  
//   TPCParametersImpl::TPCParametersImpl( GearParameters& gearParameters )  : 
//     _padLayout(0),
//     _maxDriftLength (0.),
//     _driftVelocity (0.),
//     _readoutFrequency (0.) {
    
//     _intMap         = gearParameters._intMap         ;  
//     _doubleMap      = gearParameters._doubleMap      ;
//     _stringMap      = gearParameters._stringMap      ;
//     _intVecMap      = gearParameters._intVecMap      ; 
//     _doubleVecMap   = gearParameters._doubleVecMap   ;
//     _stringVecMap   = gearParameters._stringVecMap   ;
		      		     
//     _intKeys        = gearParameters._intKeys        ;
//     _doubleKeys     = gearParameters._doubleKeys     ;
//     _stringKeys     = gearParameters._stringKeys     ;
//     _intVecKeys     = gearParameters._intVecKeys     ; 
//     _doubleVecKeys  = gearParameters._doubleVecKeys  ;
//     _stringVecKeys  = gearParameters._stringVecKeys  ;
    
//   }
  



  TPCParametersImpl::~TPCParametersImpl() {
    if( _padLayout != 0 )
      delete _padLayout ;
  }
  
  const PadRowLayout2D & TPCParametersImpl::getPadLayout() const {
    return *_padLayout ;
  }
  
  double TPCParametersImpl::getMaxDriftLength() const { 
    return _maxDriftLength ;
  } 
  
  double TPCParametersImpl::getDriftVelocity() const { 
    return _driftVelocity ;
  } 
  
  double TPCParametersImpl::getReadoutFrequency() const { 
    return _readoutFrequency ;
  } 
  
  void TPCParametersImpl::setPadLayout( PadRowLayout2D * padLayout ) { 
    
    if( _padLayout != 0 )
      delete _padLayout ;
    
    _padLayout = padLayout ;
  } 
  
  void TPCParametersImpl::setMaxDriftLength( double maxDriftLength ) { 
    _maxDriftLength = maxDriftLength ;
  } 
  
  void TPCParametersImpl::setDriftVelocity( double driftVelocity ) { 
    _driftVelocity = driftVelocity ;
  } 
  
  void TPCParametersImpl::setReadoutFrequency( double readoutFrequency ) { 
    _readoutFrequency = readoutFrequency ;
  } 
  
}
