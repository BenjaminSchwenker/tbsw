#include "gearimpl/VXDParametersImpl.h"
#include <math.h>
#include <iostream>

#define _EPSILON 0.0001

namespace gear{

    // F. Gaede, 04.04.2008 : phi0 and offset have not been defined properly in the past
    // what we want is:  phi0: azimuthal angle of first ladder's normal and the offset
    //                   going in the direction of increasing phi 
    // however what is used in the code is phi0 as the angle w/ the y-axis (and wrong orientation) 
    //
    // -> we keep the internal representation of phi0 in the variable internalPhi0 and use this here


  VXDParametersImpl::VXDParametersImpl
  ( int vxdType, double shellInnerRadius, double shellOuterRadius, double shellHalfLength, double shellGap, double shellRadLength ) :
    _vxdType( vxdType ) ,
    _shellInnerRadius( shellInnerRadius ) ,
    _shellOuterRadius( shellOuterRadius ) ,
    _shellHalfLength( shellHalfLength ) ,
    _shellGap( shellGap ) ,
    _shellRadLength( shellRadLength ) {}


 bool VXDParametersImpl::isPointInVXD( Vector3D p , bool sensitive ) const 
  {
    
    // very first check for quick calculation
    if( fabs( p[2] ) > ( _shellHalfLength + _shellGap/2 ) ) {
      // Outside shell
      return false ;
    }
    if( fabs( p[2] ) < ( _shellGap/2 ) ) {
      // Inside gap
      return false ;
    }

    // get radius of point in projection 2D
    double pRadius = sqrt( p[0]*p[0] + p[1]*p[1] ) ;

    // get phi of point in projection 2D
    double pPhi = getPhiPoint( p ) ;
    
    for( int possibleLayer = 0 ; possibleLayer < _layer.getNLayers() ; possibleLayer++ ) {

      // first some quick checks to reduce cpu time
      // check if radius is smaller then this layer
      if( ( ( !sensitive ) && ( pRadius < ( _layer.getLadderDistance( possibleLayer ) - _EPSILON ) ) )
	  ||( ( sensitive ) && ( pRadius < ( _layer.getSensitiveDistance( possibleLayer ) - _EPSILON ) ) ) ) {
	
//	std::cout << "Layer " << possibleLayer << " skipped - pRadius " << pRadius << " < " << std::endl ;
	continue ;
      }

      // check if radius is bigger than max radius for this layer
      if( pRadius > ( _layer.getMaxRadius( possibleLayer , sensitive ) + _EPSILON ) ) {
	
//	std::cout << "Layer " << possibleLayer << " skipped - pRadius " << pRadius << " > " << std::endl ;
	continue ;
      }

      // check z-plane
      if( ( ( !sensitive ) && ( fabs( p[2] ) > ( _layer.getLadderLength( possibleLayer ) + _shellGap/2 + _EPSILON ) ) )
	  ||( ( sensitive ) && ( fabs( p[2] ) > (_layer.getSensitiveLength( possibleLayer ) + _shellGap/2 + _EPSILON ) ) ) ) {
	
//	std::cout << "Layer " << possibleLayer << " skipped - z " << p[2] << " ooR " << std::endl ;
	continue ;
      }
      
      // get extensions for layer
      double zStart, zEnd, left, right, thick, d ;
      if( !sensitive ) {
	zEnd = _layer.getLadderLength( possibleLayer ) + _shellGap/2 ;
	left = -_layer.getLadderWidth( possibleLayer ) / 2 - _layer.getLadderOffset( possibleLayer ) ; //-
	right = _layer.getLadderWidth( possibleLayer ) / 2 - _layer.getLadderOffset( possibleLayer ) ;
	thick = _layer.getLadderThickness( possibleLayer ) ;
	d = _layer.getLadderDistance( possibleLayer ) ;
      }
      if (sensitive ) {
	zEnd = _layer.getSensitiveLength( possibleLayer ) + _shellGap/2 ;
	left = -_layer.getSensitiveWidth( possibleLayer ) / 2 - _layer.getSensitiveOffset( possibleLayer ) ;
	right = _layer.getSensitiveWidth( possibleLayer ) / 2 - _layer.getSensitiveOffset( possibleLayer ) ;
	thick = _layer.getSensitiveThickness( possibleLayer ) ;
	d = _layer.getSensitiveDistance( possibleLayer ) ;
      }

      zStart = _shellGap/2 ;
      
      // find out delta phi between each ladder
      double deltaPhi = ( 2 * M_PI ) / _layer.getNLadders( possibleLayer ) ;
      double startPhi = _layer.getStartInnerPhi( possibleLayer , sensitive ) ;
      double endPhi = _layer.getEndInnerPhi( possibleLayer , sensitive ) ;

      // go through all ladders/sensitive
      for( int i = 0 ; i < _layer.getNLadders( possibleLayer ) ; i++ ) {

	double phi = _layer.getInternalPhi0( possibleLayer ) + i*deltaPhi ;
	phi = correctPhiRange( phi ) ;
	
	// check start and end phi for this layer
	double sPhi = correctPhiRange( startPhi + i*deltaPhi - _EPSILON ) ;
	double ePhi = correctPhiRange( endPhi + i*deltaPhi + _EPSILON ) ;
			
	// check if point is in range of ladder
	bool isInRange = 
	  ( ( sPhi <= pPhi )  && ( pPhi <= ePhi ) )
	  ||( ( sPhi > ePhi ) && ( pPhi <= ePhi ) )
	  ||( ( sPhi > ePhi ) && ( pPhi >= sPhi ) ) ;

	if ( !isInRange ) {
	  // std::cout << " ladder/layer " << i << "/" << possibleLayer << " skipped ooR" << std::endl ;
	  continue ;
	}
	
	// take normal vector of planes, Base, Side and z
	Vector3D nBase, nSide, nZ ;
	
	nBase[0] = sin( phi ) ;
	nBase[1] = cos( phi ) ;
	nBase[2] = 0;

	nSide[0] = cos( phi ) ;
	nSide[1] = -sin( phi ) ;
	nSide[2] = 0 ;
      
	nZ[0] = 0 ;
	nZ[1] = 0 ;
	nZ[2] = 1 ;

	if( p[2] < 0 ) {
	  nZ[2] = - nZ[2] ;
	}
	
	Vector3D r( d * nBase[0] , d * nBase[1] , d * nBase[2] ) ;
	
	Vector3D thisVector = distanceToPlane( p, r, nBase, nSide, nZ, left, right, zStart, zEnd ) ;

	
	// calculate distance
	double c = sqrt( thisVector[0] * thisVector[0] + thisVector[1] * thisVector[1] + thisVector[2] * thisVector[2] ) ;

	// debug
//  	std::cout << "**n[0] :" << nBase[0] << " n[1] " <<nBase[1] <<" n[2] " <<nBase[2] << std::endl ;
// 	std::cout << "**nv[0]:" << thisVector[0]  << " nv[1]:" << thisVector[1]  << "nv[2]:" << thisVector[2]  << std::endl ;
//  	std::cout << "**c :" << c << "  thick: " << thick << std::endl ; // debug
//  	std::cout << "**sPhi : " << sPhi*180/M_PI <<" pPhi " << pPhi*180/M_PI << " ePhi: " << ePhi*180/M_PI << std::endl ;
//  	std::cout << "**Phi : " << phi*180/M_PI << " atan() : " << atan( p[0]/p[1] ) << std::endl ;
// 	std::cout << "**isEqual[0] " << isEqual( thisVector[0] , -c*nBase[0] ) 
// 		  << " , " << isEqual( thisVector[1] , -c*nBase[1] ) 
// 		  << " , " << isEqual( thisVector[2] , -c*nBase[2] ) 
// 		  << std::endl ;
		
	// check if vector is same direction as nBase and smaller than thick
	// skip vectors if Basevector close to zero
	if ( ( differsLess( nBase[0] , thisVector[0] ) || isEqual( thisVector[0] , -c*nBase[0] ) )
	     && ( differsLess( nBase[1] , thisVector[1] ) || isEqual( thisVector[1] , -c*nBase[1] ) )
	     && ( differsLess( nBase[2] , thisVector[2] ) || isEqual( thisVector[2] , -c*nBase[2] ) )
	     && ( ( c <= thick ) || isEqual( c , thick ) ) ) {
	  return true ;
	}
	
      } // all ladders

    } // all layers
	  
    return false ;

  }  // function isPointInVXD




  Vector3D VXDParametersImpl::distanceToNearestVXD( Vector3D p , bool sensitive ) const
  {

    Vector3D origin( 0,0,0 ) ;
    
    // check if point is in already
    if( isPointInVXD( p , sensitive ) ) {
      
      return origin ;
    }
    
    // check if point is origin

    if( origin.x() == p.x() &&  origin.y() == p.y() && origin.z() == p.z()  ) {

       if( _layer.getNLayers() > 0 ){

 	// simply use first ladder - distance d , phi0 , gap/2  :
	 //fg: fixme: review definition of phi...
 	return Vector3D( _layer.getLadderDistance(0) ,  -_layer.getInternalPhi0(0)  , _shellGap/2  , Vector3D::cylindrical ) ;
       }

       return origin ;
    }

    //fg:  FIXME: this code needs to be cleaned up and checked 
    //  it can also probably be made a bit faster ...


    double pPhi = getPhiPoint( p ) ;

    // this Point is solution
    Vector3D vPointPlane ;
    double minDistance = 99999.99 ;

    // go through all layers
    for ( int nearestLayer=0 ; nearestLayer < _layer.getNLayers() ; nearestLayer++ ) {

      // get deltaPhi between ladders
      double deltaPhi = ( 2 * M_PI ) / _layer.getNLadders( nearestLayer ) ;
      
      // get extensions for layer
      double dist, zStart, zEnd, left, right, thick ;
      if( !sensitive ) {
	zEnd = _layer.getLadderLength( nearestLayer ) + _shellGap/2 ;
	left = -_layer.getLadderWidth( nearestLayer ) / 2 - _layer.getLadderOffset( nearestLayer ) ;
	right = _layer.getLadderWidth( nearestLayer ) / 2 - _layer.getLadderOffset( nearestLayer ) ;
	thick = _layer.getLadderThickness( nearestLayer ) ;
	dist = _layer.getLadderDistance( nearestLayer ) ;
      }
      if (sensitive ) {
	zEnd = _layer.getSensitiveLength( nearestLayer ) + _shellGap/2 ;
	left = -_layer.getSensitiveWidth( nearestLayer ) / 2 - _layer.getSensitiveOffset( nearestLayer ) ;
	right = _layer.getSensitiveWidth( nearestLayer ) / 2 - _layer.getSensitiveOffset( nearestLayer ) ;
	thick = _layer.getSensitiveThickness( nearestLayer ) ;
	dist = _layer.getSensitiveDistance( nearestLayer ) ;
      }

      zStart = _shellGap/2 ;

      // go through all ladders/sensitive
      for( int i = 0 ; i < _layer.getNLadders( nearestLayer ) ; i++ ) {
	
	double phi = _layer.getInternalPhi0( nearestLayer ) + i*deltaPhi ;
	phi = correctPhiRange( phi ) ;

	// check start and end phi for this layer
	double sPhi = correctPhiRange( phi - deltaPhi ) ;
	double ePhi = correctPhiRange( phi + deltaPhi ) ;
			
	// check if point is in range of ladder
	bool isInRange = 
	  ( ( sPhi <= pPhi )  && ( pPhi <= ePhi ) )
	  ||( ( sPhi > ePhi ) && ( pPhi <= ePhi ) )
	  ||( ( sPhi > ePhi ) && ( pPhi >= sPhi ) ) ;
	
	if ( !isInRange ) {
	  
	  continue ;
	}
	
	// take normal vector of planes, Base, Side, Z and spacePoint
//fg: ------- test vector3d ---------------------------
// 	Vector3D nBase, nSide, nZ ;
	
// 	nBase[0] = sin( phi ) ;
// 	nBase[1] = cos( phi ) ;
// 	nBase[2] = 0;
      
// 	nSide[0] = cos( phi ) ;
// 	nSide[1] = -sin( phi ) ;
// 	nSide[2] = 0 ;
      
// 	nZ[0] = 0 ;
// 	nZ[1] = 0 ;
// 	nZ[2] = 1 ;

// 	if( p[2] < 0 ) {
// 	  nZ[2] = - nZ[2] ;
// 	}
	
// 	Vector3D spacePoint ;
// 	spacePoint[0] = dist * nBase[0] ;
// 	spacePoint[1] = dist * nBase[1] ;
// 	spacePoint[2] = dist * nBase[2] ;
	

	Vector3D nBase( sin( phi ),  cos( phi ), 0.0 ) ;
	Vector3D nSide( cos( phi ), -sin( phi ), 0.0 ) ;

	Vector3D nZ( 0. , 0. , 1. ) ;
 	if( p[2] < 0 ) {
 	  nZ[2] = -nZ[2] ;
 	}
	Vector3D spacePoint = dist * nBase ;

//fg: ------- end test vector3d ---------------------------



	// upper and lower base
	// inner (1) and outer (2) boundary
	for ( int j = 1 ; j<=2 ; j ++ ) {
	  
	  double d = 0 ;
	  if( j == 2 ) {
	    d = thick ;
	  }

//fg: ------- test vector3d ---------------------------
// 	  Vector3D r ;
// 	  r[0] = spacePoint[0] + d * nBase[0] ;
// 	  r[1] = spacePoint[1] + d * nBase[1] ;
// 	  r[2] = spacePoint[2] + d * nBase[2] ;
	  Vector3D r = spacePoint + d * nBase ;
//fg: ------- end test vector3d ---------------------------
  

	  Vector3D myVec = distanceToPlane( p, r, nBase, nSide, nZ, left, right, zStart, zEnd ) ;

	  // debug
// 	  Vector3D iP ( p[0]+myVec[0] , p[1]+myVec[1] , p[2]+myVec[2] ) ;
// 	  bool isCorrect = isPointInVXD( iP ) ;
	  
// 	  if ( !isCorrect ) {
// 	    std::cout << "distanceToPlane returns sth wrong. j =" <<j << std::endl ;
// 	  }
	  
	  double thisDistance = sqrt( myVec[0] * myVec[0] + myVec[1] * myVec[1] + myVec[2] * myVec[2] ) ;
	  if ( thisDistance <= minDistance ) {
	    minDistance = thisDistance ;
	    vPointPlane = myVec ;
	    
// 	    // debug
//  	    std::cout << "\nmyVec set " <<"base " << myVec[0] << " , " << myVec[1] << " , " << myVec[2] << std::endl ; //debug
//  	    std::cout << "      lay " << nearestLayer << " d " << d << " phi " << phi * 180/M_PI << " j " << j <<std::endl ;
// 	    std::cout << "        r " << r[0] << " , " << r[1] << " , " << r[2] 
// 		      << " nSide " << nSide[0] << " , " << nSide[1] << " , " << nSide[2] 
// 		      << " <-->  " <
// 	              << left << " - " << right << std::endl;
	  }
	  
	} // for j -- lower & upper boundary
	

	// left and right side
	// 1 left, 2 right
	
	for ( int j = 1 ; j<=2 ; j ++ ) {
	  
	  double d = 0;
	  if ( (!sensitive) and ( j == 1 ) ) {
	    d = - _layer.getLadderWidth( nearestLayer ) / 2 - _layer.getLadderOffset( nearestLayer ) ;
	  }
	  if( (!sensitive) and ( j == 2 ) ) {
	    d = _layer.getLadderWidth( nearestLayer ) / 2 - _layer.getLadderOffset( nearestLayer ) ;
	  }
	  if( (sensitive) and ( j == 1 ) ) {
	    d = - _layer.getSensitiveWidth( nearestLayer ) / 2 - _layer.getSensitiveOffset( nearestLayer ) ;
	  }
	  if( (sensitive) and ( j == 2 ) ) {
	    d = _layer.getSensitiveWidth( nearestLayer ) / 2 - _layer.getSensitiveOffset( nearestLayer ) ;
	  }
	  
	  // lower left corner as r
	  Vector3D r ;
	  r[0] = spacePoint[0] + d * nSide[0] ; 
	  r[1] = spacePoint[1] + d * nSide[1] ;
	  r[2] = spacePoint[2] + d * nSide[2] ;
	  
	  Vector3D myVec = distanceToPlane( p, r, nSide, nBase, nZ, 0, thick, zStart, zEnd ) ;
	  
	  double thisDistance = sqrt( myVec[0] * myVec[0] + myVec[1] * myVec[1] + myVec[2] * myVec[2] ) ;
	  if ( thisDistance < minDistance ) {
	    minDistance = thisDistance ;
	    vPointPlane = myVec ;
// 	    std::cout << "myVec set " <<"left/right " << myVec[0] << " , " << myVec[1] << " , " << myVec[2] << std::endl ; //debug
	  }
	  
	} // for j -- left and right corner
	
	
	      
	// z Front and Back
	for( int j = 1 ; j <= 2 ; j++ ) {
	  
	  double d = 0 ;
	  
	  if( j == 1 ) {
	    d = zEnd ;
	  }
	  else {
	    d = zStart ;
	  }
	  
	  // lower middle corner as r
	  Vector3D r;
	  r[0] = spacePoint[0] + d * nZ[0] ; 
	  r[1] = spacePoint[1] + d * nZ[1] ;
	  r[2] = spacePoint[2] + d * nZ[2] ;

	  Vector3D myVec = distanceToPlane( p, r, nZ, nBase, nSide, 0, thick, left, right ) ;
	  
	  double thisDistance = sqrt( myVec[0] * myVec[0] + myVec[1] * myVec[1] + myVec[2] * myVec[2] ) ;

	  if ( thisDistance < minDistance ) {
	    minDistance = thisDistance ;
	    vPointPlane = myVec ;
// 	    std::cout << "myVec set " <<" front/back " << myVec[0] << " , " << myVec[1] << " , " << myVec[2] << std::endl ; //debug
// 	    std::cout << "    thick " << thick << " left " << left << " right " << right << std::endl ;
// 	    std::cout << "        r " << r[0] << " , " << r[1] << " , " << r[2] << std::endl ;
	  }
	  
	} // for j -- front and back

      } // for i -- all ladders in layer

    } // for nearestLayer -- all layers in vxd


    return vPointPlane ;

  } // function distanceToVXD





  Vector3D VXDParametersImpl::intersectionVXD( Vector3D p , Vector3D v , bool sensitive ) const {
    
    // return values
    double currentLength = 99999 ;
    Vector3D currentPoint (0. , 0. , 0. ) ;
    
    for ( int takeLayer=0 ; takeLayer < _layer.getNLayers() ; takeLayer++ ) {

      // get deltaPhi between ladders
      double deltaPhi = ( 2 * M_PI ) / _layer.getNLadders( takeLayer ) ;

      // get extensions for layer
      double zStart, zEnd, left, right, thick ;
      if( !sensitive ) {
	zEnd = _layer.getLadderLength( takeLayer ) + _shellGap/2 ;
	left = -_layer.getLadderWidth( takeLayer ) / 2 - _layer.getLadderOffset( takeLayer ) ;
	right = _layer.getLadderWidth( takeLayer ) / 2 - _layer.getLadderOffset( takeLayer ) ;
	thick = _layer.getLadderThickness( takeLayer ) ;
      }
      if (sensitive ) {
	zEnd = _layer.getSensitiveLength( takeLayer ) + _shellGap/2 ;
	left = -_layer.getSensitiveWidth( takeLayer ) / 2 - _layer.getSensitiveOffset( takeLayer ) ;
	right = _layer.getSensitiveWidth( takeLayer ) / 2 - _layer.getSensitiveOffset( takeLayer ) ;
	thick = _layer.getSensitiveThickness( takeLayer ) ;
      }

      // go through both sides
      for( int side = 0 ; side <2 ; side++ ) {
	if( side == 0 ) {
	  zStart = _shellGap /2 ;
	}
	else {
	  zEnd = - zEnd ;
	  zStart = - _shellGap/2 ;
	}

	// go through all ladders/sensitive
	for( int i = 0 ; i < _layer.getNLadders( takeLayer ) ; i++ ) {
	  
	  double phi = _layer.getInternalPhi0( takeLayer ) + i*deltaPhi ;
	  phi = correctPhiRange( phi ) ;
	  
	  // take normal vector of planes, Base, Side, Z and spacePoint
	  Vector3D nBase, nSide, nZ, spacePoint;
	  
	  nBase[0] = sin( phi ) ;
	  nBase[1] = cos( phi ) ;
	  nBase[2] = 0;
	  
	  nSide[0] = cos( phi ) ;
	  nSide[1] = -sin( phi ) ;
	  nSide[2] = 0 ;
      
	  nZ[0] = 0 ;
	  nZ[1] = 0 ;
	  nZ[2] = 1 ;
	
	  // upper and lower base
	  // inner (1) and outer (2) boundary
	  for ( int j = 1 ; j<=2 ; j ++ ) {
	    
	    double d = 0 ;
	    if ( (!sensitive) and ( j == 1 ) ) {
	      d = _layer.getLadderDistance( takeLayer ) ;
	    }
	    if( (!sensitive) and ( j == 2 ) ) {
	      d = _layer.getLadderDistance( takeLayer ) + thick ;
	    }
	    if( (sensitive) and ( j == 1 ) ) {
	      d = _layer.getSensitiveDistance( takeLayer ) ;
	    }
	    if( (sensitive) and ( j == 2 ) ) {
	      d = _layer.getSensitiveDistance( takeLayer ) + thick ;
	    }
	    Vector3D r( d * nBase[0] , d * nBase[1] , d * nBase[2] ) ;
	    
	    // this vector is equivalent to spacepoint -- needed later on
	    if ( j == 1 ) spacePoint = r ;

	    // get intersection point plane/line	    
	    Vector3D interSection = planeLineIntersection( r, nBase, p, v ) ;

	    // get vector from point to intersection point
	    Vector3D vip ;
	    vip[0] = interSection[0] - p[0] ;
	    vip[1] = interSection[1] - p[1] ;
	    vip[2] = interSection[2] - p[2] ;

	    // get vector from spacepoint to intersection point
	    Vector3D vn ;
	    vn[0] = interSection[0] - spacePoint[0] ;
	    vn[1] = interSection[1] - spacePoint[1] ;
	    vn[2] = interSection[2] - spacePoint[2] ;
  
	    // correct in Plane
	    Vector3D inBorder = correctToBorderPoint( vn, nSide, nZ, left, right, zStart, zEnd ) ;
	    
	    // check if point has not changed -> intersects with vxd
	    if ( isEqual( inBorder , vn ) ) {
	      
	      // check if this is closest intersection
	      double thisLength = sqrt( vip[0]*vip[0] + vip[1]*vip[1] +vip[2]*vip[2] ) ;
	      if( thisLength < currentLength ) {
		currentPoint = interSection ;
		currentLength = thisLength ;
	      }
	    }
	    
	  } // for j -- lower & upper boundary
	  

	  // left and right side
	  // 1 left, 2 right
	  
	  for ( int j = 1 ; j<=2 ; j ++ ) {
	    
	    double d = 0;
	    if ( (!sensitive) and ( j == 1 ) ) {
	      d = -( _layer.getLadderWidth( takeLayer ) / 2 + _layer.getLadderOffset( takeLayer ) );
	    }
	    if( (!sensitive) and ( j == 2 ) ) {
	      d = _layer.getLadderWidth( takeLayer ) / 2 - _layer.getLadderOffset( takeLayer ) ;
	    }
	    if( (sensitive) and ( j == 1 ) ) {
	      d = -( _layer.getSensitiveWidth( takeLayer ) / 2 + _layer.getSensitiveOffset( takeLayer ) ) ;
	    }
	    if( (sensitive) and ( j == 2 ) ) {
	      d = _layer.getSensitiveWidth( takeLayer ) / 2 - _layer.getSensitiveOffset( takeLayer ) ;
	    }
	    
	    // lower left corner as r
	    Vector3D r ;
	    r[0] = spacePoint[0] + d * nSide[0] ; 
	    r[1] = spacePoint[1] + d * nSide[1] ;
	    r[2] = spacePoint[2] + d * nSide[2] ;

	    // get intersection point plane/line
	    Vector3D interSection = planeLineIntersection( r, nSide, p, v ) ;
	     
	    // get vector from point to intersection point
	    Vector3D vip ;
	    vip[0] = interSection[0] - p[0] ;
	    vip[1] = interSection[1] - p[1] ;
	    vip[2] = interSection[2] - p[2] ;
	    
	    // get vector from spacepoint to intersection point
	    Vector3D vn ;
	    vn[0] = interSection[0] - spacePoint[0] ;
	    vn[1] = interSection[1] - spacePoint[1] ;
	    vn[2] = interSection[2] - spacePoint[2] ;
  
	    // correct in Plane
	    Vector3D inBorder = correctToBorderPoint( vn, nBase, nZ, 0, thick, zStart, zEnd ) ;

	    if ( isEqual( inBorder , vn )) {
	      // check if this is closest intersection
	      double thisLength = sqrt( vip[0]*vip[0] + vip[1]*vip[1] + vip[2]*vip[2] ) ;
	      if( thisLength < currentLength ) {
		currentPoint = interSection ;
		currentLength = thisLength ;
	      }
	      
	    }

	    
	  } // for j -- left and right corner
	
	  
	  // front and back side in z
	  
	  for( int j = 1 ; j <= 2 ; j++ ) {
	    
	    double d = 0 ;
	    
	    if( j == 1 ) {
	      d = zEnd ;
	    }
	    else {
	      d = zStart ;
	    }
	    
	    // lower left corner as r
	    Vector3D r;
	    r[0] = spacePoint[0] + d * nZ[0] ; 
	    r[1] = spacePoint[1] + d * nZ[1] ;
	    r[2] = spacePoint[2] + d * nZ[2] ;
	    
	    // get intersection point plane/line
	    Vector3D interSection = planeLineIntersection( r, nZ, p, v ) ;
	    
	    // get vector from point to intersection point
	    Vector3D vip ;
	    vip[0] = interSection[0] - p[0] ;
	    vip[1] = interSection[1] - p[1] ;
	    vip[2] = interSection[2] - p[2] ;
	    
	    // get vector from spacepoint to intersection point
	    Vector3D vn ;
	    vn[0] = interSection[0] - spacePoint[0] ;
	    vn[1] = interSection[1] - spacePoint[1] ;
	    vn[2] = interSection[2] - spacePoint[2] ;
  
	    // correct in Plane
	    Vector3D inBorder = correctToBorderPoint( vn, nBase, nSide, 0, thick, left, right ) ;

	    if ( isEqual( inBorder , vn ) ) {
	      // check if this is closest intersection
	      double thisLength = sqrt( vip[0]*vip[0] + vip[1]*vip[1] + vip[2]*vip[2] ) ;
	      if( thisLength < currentLength ) {
		currentPoint = interSection ;
		currentLength = thisLength ;
	      }
	    }
	    
	  } // for j -- front and back
	  

	} // for i -- all ladders in layer
	
      } // for all zSides

    } // for nearestLayer -- all layers in vxd

      return currentPoint ;

  } // function isIntersectionVXD



  Vector3D VXDParametersImpl::distanceToPlane( Vector3D p, Vector3D r, Vector3D n, Vector3D u, Vector3D v, float minU, float maxU, float minV, float maxV ) const {

    // check zeros

    if ( n[0]==0 && n[1]==0 && n[2]==0 ) {
      
      std::cout << "\ndistanceToPlane - Fatal error!" << std::endl ;
      std::cout << "normal vector (0,0,0)" << std::endl ;
      return Vector3D (0,0,0) ;
    }
    
    // check u and v
    if ( !isEqual( n[0]*( u[0]+v[0] ) + n[1]*( u[1]+v[1] ) + n[2]*( u[2]+v[2] ) , 0 ) ){
      std::cout << "\ndistanceToPlan - Fatal error! " 
		<< "n not normal to u and v " << std::endl ;
      return Vector3D (0,0,0) ;
    }
  
    //  The normal n to the plane is n(n0,n1,n2).
    //  The nearest point in the plane is pn(pn0, pn1, pn2)
    //
    //  The line defined by (pn0-p0)/n0 = (pn1-p1)/n1 = (pn2-p2)/n2 = t
    //  goes through p(p0,p1,p2) and is parallel to n.
    //
    //  Solving for the point pn we get
    //
    //    pn0 = n0*t+p0
    //    pn1 = n1*t+p1
    //    pn2 = n2*t+p2
    //
    //  Now place these values in the equation for the plane:
    //
    //    n0 * (n0*t+p0) + n1*(n1*t+p1) + n2*(n2*t+p2) + r = 0
    //
    //  and solve for t:
    //
    //    t = (-n0*p0-n1*p1-n2*p2-d) / ( n0 * n0 + n1 * n1 + n2 * n2 )
 
 
   //  the distance to plane is Vector to spacepoint scalar with n

    double d = ( r[0] * n[0] + r[1] * n[1] + r[2] * n[2] ) ;

    double t = - ( n[0] * p[0] + n[1] * p[1] + n[2] * p[2] - d ) ;

    // nearest point in plane
    Vector3D pn ;
    pn[0] = p[0] + n[0] * t;
    pn[1] = p[1] + n[1] * t;
    pn[2] = p[2] + n[2] * t;

    // vector from spacespoint to nearest point in plane
    Vector3D vn ;
    vn[0] = pn[0] - r[0] ;
    vn[1] = pn[1] - r[1] ;
    vn[2] = pn[2] - r[2] ;

    // correct to borderPoint
    Vector3D vCor = correctToBorderPoint( vn , u, v, minU, maxU, minV, maxV ) ;

    // final vector from point to plane
    Vector3D final ;
    final[0] = r[0] + vCor[0] - p[0] ;
    final[1] = r[1] + vCor[1] - p[1] ;
    final[2] = r[2] + vCor[2] - p[2] ;

    return final ;

  } // function distanceToPlane



  Vector3D VXDParametersImpl::planeLineIntersection( Vector3D r, Vector3D n, Vector3D linePoint, Vector3D lineDir) const {

    // calculate distance of plane
    double distance = ( r[0]*n[0] + r[1]*n[1] + r[2]*n[2] ) / ( sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] ) ) ;
    
    // equations for a point defined by the line with distance to plane = 0
    // n * ( linePoint + t * lineDir ) - distance = 0
    // solve for t

    double numerator   = distance - n[0]*linePoint[0] - n[1]*linePoint[1] - n[2]*linePoint[2] ;
    double denominator = n[0]*lineDir[0] + n[1]*lineDir[1] + n[2]*lineDir[2] ;
    
    if ( denominator == 0 ) {
      // line paralell to plane -> no intersection
      return Vector3D ( 0. , 0. , 0. ) ;
    }

    double t = numerator/denominator ;

    if( t < 0 ) {
      // intersection not in direction of line
      return Vector3D ( 0. , 0. , 0. ) ;
    }

    Vector3D final ;
    final[0] = linePoint[0] + t * lineDir[0] ;
    final[1] = linePoint[1] + t * lineDir[1] ;
    final[2] = linePoint[2] + t * lineDir[2] ;

    return final ;

  }// function planLineIntersectin



  Vector3D VXDParametersImpl::correctToBorderPoint( Vector3D vPlane , Vector3D u, Vector3D v, float minU, float maxU, float minV, float maxV ) const {

    // check if the given vector is within the boundaries
    double neededUs = vPlane[0] * u[0] + vPlane[1] * u[1] + vPlane[2] * u[2] ;
    double neededVs = vPlane[0] * v[0] + vPlane[1] * v[1] + vPlane[2] * v[2] ;

   
    double goUs = neededUs, goVs = neededVs ;
    bool changed = false ;

    if( neededUs < minU ) {
      goUs = minU ;
      changed = true ;
    }
    if( neededUs > maxU ) {
      goUs = maxU ;
      changed = true ;
    }
    if( neededVs < minV ) {
      goVs = minV ;
      changed = true ;
    }
    if( neededVs > maxV ) {
      goVs = maxV ;
      changed = true ;
    }

    if( !changed ) return vPlane ;
    
    Vector3D final ;
    final[0] = goUs * u[0] + goVs * v[0] ;
    final[1] = goUs * u[1] + goVs * v[1] ;
    final[2] = goUs * u[2] + goVs * v[2] ;

    return final ;

  } // correctToBorderPoint

  
  double VXDParametersImpl::correctPhiRange( double Phi ) const {
    
    if( Phi > M_PI ) {
      return ( Phi - 2 * M_PI ) ;
    }
    if( Phi <= -M_PI ) {
      return ( Phi + 2 * M_PI ) ;
    } 
    
    return Phi ;

  } // function correctPhiRange
  

  double VXDParametersImpl::getPhiPoint( Vector3D p ) const {

    //FG: this is the the angle with the negative y-axis 
    // see comment at the top !

    // get phi of point in projection 2D
    double pPhi = 0. ;
    if( ( p[0] >= 0 ) && ( p[1] == 0 ) )
      pPhi = -M_PI/2 ;
    
    if( ( p[0] < 0 ) && ( p[1] == 0 ) )
      pPhi = M_PI/2 ;
    
    if( ( p[0] == 0 ) && ( p[1] < 0 ) ) 
      pPhi = 0 ;
    
    if( ( p[0] != 0 ) && ( p[1] < 0 ) )
      pPhi = atan( p[0] / p[1] ) + M_PI ;
    
    else
      pPhi = atan( p[0] / p[1] ) ;
    
    pPhi = correctPhiRange( pPhi ) ;  
    
    return pPhi ;
    
  } // function getPhiPoint



  bool VXDParametersImpl::isEqual( double valueOne, double valueTwo ) const
  {

    // save calculating time if equal
    if ( valueOne == valueTwo ) return true ;


    // if one value is 0.0 we use the absolute distance 
    if( valueOne * valueTwo == 0.0 )
	       
      return fabs( valueOne - valueTwo ) < _EPSILON ;
    
    
    // if both values are smaller than epsilon they are equal for our purpose 
    if ( fabs(valueOne)  < _EPSILON && fabs(valueTwo)  < _EPSILON )
      
      return true ;
    

    // otherwise we require the relative distance of the two values to be smaller than epsilon
    return fabs( 2 * (valueOne - valueTwo) ) /  (fabs(valueOne) + fabs(valueTwo) ) <  _EPSILON  ;
    
  } 


  bool VXDParametersImpl::isEqual( Vector3D p1 , Vector3D p2 ) const
  {
    for( int i=0 ; i<3 ; i++ ) {
      if( !( p1[i] == p2[i] ) ) {
	return false ;
      }
    }
    return true ;
  } // function isEqual
  

  bool VXDParametersImpl::differsLess( double valueOne, double valueTwo ) const
  {

    return ( ( fabs( valueOne ) + fabs( valueTwo ) ) < _EPSILON ) ;

  } // function differsLess
  
 
} // namespace



