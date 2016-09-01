#include "gearimpl/VXDLayerLayoutImpl.h"
#include <math.h>

namespace gear{

  void VXDLayerLayoutImpl::addLayer(int nLadders, double phi0, 
				    // ladder
				    double ladderDistance, double ladderOffset, double ladderThickness,
				    double ladderLength, double ladderWidth, double ladderRadLength,
				    // sensitive
				    double sensitiveDistance, double sensitiveOffset, double sensitiveThickness,
				    double sensitiveLength, double sensitiveWidth, double sensitiveRadLength )
  {

    Layer lL, sL ;
    
    // F. Gaede, 04.04.2008 : phi0 and offset have not been defined properly in the past
    // what we want is:  phi0: azimuthal angle of first ladder's normal and the offset
    //                   going in the direction of increasing phi 
    // however what is used in the code is phi0 as the angle w/ the y-axis (and wrong orientation) 
    //
    // -> we keep the internal representation of phi0 in the variable internalPhi0 and return
    //    M_PI/2. - internalPhi0
    
    lL.NLadders  = nLadders ;
    lL.internalPhi0      = M_PI/2. - phi0 ;   // phi0 given wrt x-axis  - internalPhi0 is wrt yaxis , negative orientation
    lL.Distance  = ladderDistance ;
    lL.Offset    = ladderOffset ;   // offset given in direction of positive rotation (increasing phi)
    lL.Thickness = ladderThickness ;
    lL.Length    = ladderLength ;
    lL.Width     = ladderWidth ;
    lL.RadLength = ladderRadLength ;

    sL.NLadders  = nLadders ;
    sL.internalPhi0      = M_PI/2. - phi0 ;  // phi0 given wrt x-axis  - internalPhi0 is wrt yaxis , negative orientation
    sL.Distance  = sensitiveDistance ;
    sL.Offset    = sensitiveOffset ;  // offset given in direction of positive rotation (increasing phi)
    sL.Thickness = sensitiveThickness ;
    sL.Length    = sensitiveLength ;
    sL.Width     = sensitiveWidth ;
    sL.RadLength = sensitiveRadLength ;

    _lVec.push_back( lL ) ;
    _sVec.push_back( sL ) ;

  }

      
  double VXDLayerLayoutImpl::getMaxRadius( int layerIndex, bool sensitive ) const 
  {
    Layer l ;
    if (!sensitive) {
      l = _lVec.at( layerIndex ) ;
    }
    else {
      l = _sVec.at( layerIndex ) ;
    }

    float d = l.Distance + l.Thickness ;
    float w = l.Width / 2 + fabs( l.Offset ) ;

    double max = sqrt( d*d + w*w ) ;

    return max ;
  }

  // returns starting phi for first ladder in layer (on side to IP)
  double VXDLayerLayoutImpl::getStartInnerPhi( int layerIndex , bool sensitive ) const 
  {
    Layer l ;
    if ( !sensitive ) {
      l = _lVec.at( layerIndex ) ;
    }
    else {
      l = _sVec.at( layerIndex ) ;
    }
    return ( l.internalPhi0 + atan( (-l.Width /2 - l.Offset) / l.Distance) ) ;
  }

  // returns ending phi for first ladder in layer (on side to IP)
  double VXDLayerLayoutImpl::getEndInnerPhi( int layerIndex , bool sensitive ) const 
  {
    Layer l ;
    if ( !sensitive ) {
      l = _lVec.at( layerIndex ) ;
    }
    else {
      l = _sVec.at( layerIndex ) ;
    }
    return ( l.internalPhi0 + atan( (l.Width/2 - l.Offset)  / l.Distance ) ) ;
  }


  // returns starting phi for first ladder in layer (on side away from IP)
  double VXDLayerLayoutImpl::getStartOuterPhi( int layerIndex , bool sensitive ) const 
  {
    Layer l ;
    if (sensitive) {
      l = _lVec.at( layerIndex ) ;
    }
    else {
      l = _sVec.at( layerIndex ) ;
    }
    return ( l.internalPhi0 - atan( (l.Width /2 + l.Offset) / (l.Distance + l.Thickness) ) )  ;
  }

  // returns ending phi for first ladder in layer (on side away from IP)
  double VXDLayerLayoutImpl::getEndOuterPhi( int layerIndex , bool sensitive ) const 
  {
    Layer l ;
    if (sensitive) {
      l = _lVec.at( layerIndex ) ;
    }
    else {
      l = _sVec.at( layerIndex ) ;
    }
    return ( l.internalPhi0 + atan( (l.Width/2 - l.Offset)  / (l.Distance + l.Thickness) ) ) ;
  }

  // returns thickness under a certain angle
  double VXDLayerLayoutImpl::getThicknessForPhi( int layerIndex , double phi, bool sensitive ) const
  {
    Layer l ;
    if (!sensitive) {
      l = _lVec.at( layerIndex ) ;
    }
    else {
      l = _sVec.at( layerIndex ) ;
    }

    double angularThickness ;

    // first check if layer is completely out of ladder
    if( phi < getStartInnerPhi( layerIndex , sensitive ) || phi > getEndInnerPhi( layerIndex, sensitive ) ) {
      return -1 ;
    }

    // check if angle is withhin outer boundaries - easy calculation then
    if( phi >= getStartOuterPhi( layerIndex , sensitive ) && phi <= getEndOuterPhi( layerIndex , sensitive ) ) {
      return ( l.Thickness / cos( phi ) ) ;
    }

    // the angel is in corner area:

    // calculate distance from space point to intersection point with straight line from IP under phi
    double distanceSpacePoint = l.Distance * tan( phi ) ;
    
    // calculate length of edge that is cut away on side facing IP
    double cutAwayLength = l.Width/2 - l.Offset - distanceSpacePoint ;

    // now one can calculate the angular thickness
    angularThickness = cutAwayLength / sin( phi ) ;

    return angularThickness ;
  }    
      



} //namespace
