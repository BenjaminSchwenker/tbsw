// -*- C++ -*-
#ifndef GEAR_VXDParametersImpl_H
#define GEAR_VXDParametersImpl_H 1

#include "gear/VXDParameters.h"
#include "gearimpl/GearParametersImpl.h"
#include "gearimpl/VXDLayerLayoutImpl.h"


namespace gear {

class VXDLayerLayout;

  /** Geometry properties of a vertex detector needed for reconstruction code. <br>
   *  <p>The vertex is assumed to consist of a number of layers. Each layer consists of
   *  a number of rectangular ladders that are uniformly distributed in a circle around the IP.<br>
   *  The sensitive volumes can be placed relative to the (insensitive) ladders.
   *  @see addLayer .</p>
   *  <p>The shell is described by outer and inner radius and it's half length in z.</p>
   *  <p>The gap is assumed to be symetrical around z==0 and goes through all layers and ladders.</p>
   *  
   *
   *  @author R. Lippe, DESY
   *  @version $Id: 
   */
class VXDParametersImpl : public GearParametersImpl, public VXDParameters {

public: 
  /** C'tor  
   *  @param vxdType           the type of the vertex detector: CCD, CMOS or HYBRID
   *  @param shellInnerRadius  the inner Radius of the vertex' shell (in mm)
   *  @param shellOuterRadius  the outer Radius of the vertex' shell (in mm)
   *  @param shellHalfLength   the half length in z (in mm)
   *  @param shellGap          the total width of the gap at z==0 - symetrical around z==0
   *  @param shellRadLength    the material property information about the shell's material radiation length
   */
  VXDParametersImpl( int vxdType, double shellInnerRadius, double shellOuterRadius, double shellHalfLength, double shellGap, double shellRadLength ) ;

  // Destructor.
  virtual ~VXDParametersImpl() { /* nop */; }
  
  // static constants enums( for the type of layout )
  
  static const int CCD    = 1 ;
  static const int CMOS   = 2 ;
  static const int HYBRID = 3 ;

  /** adding a Layer to the vertex
   *  one layer consiste of a number of ladders and sensitive areas (nLadders) that are
   *  each equidistant to IP. They are cubic with dimensions length (z), width (perpendicular to
   *  connection IP-spacepoint) and thickness (extension to connection IP-spacepoint).
   *
   * @param nLadders   the number of ladders and sensitive areas inside the layer        
   * @param phi0       azimuthal angle of normal to first ladder
   * @param Distance   the length (mm) of the straigth line between IP and spacepoint to 
   *                   ladder (ladderDistance) or sensitive Area (sensitiveDistance)     
   * @param Offset     the offset of ladder (ladderOffset) or sensitive area             
   *                   (sensitiveOffset) in direction of increasing phi in mm between spacepoint 
   *                   and middle of ladder/sensitive                                    
   * @param Thickness  the thickness in mm of the ladder (ladderThickness) or the        
   *                   sensitive area (sensitiveThickness)                               
   * @param Length     the length (z-direction) in mm of the ladder (ladderLength)       
   *                   or the sensitive area (sensitiveLength)                           
   * @param Width      the width of the ladder (ladderWidth) or the                      
   *                   sensitive area (sensitiveWidth) in mm in clockwise direction      
   * @param RadLenght  the radiation lenght of the ladder (ladderRadLength) or the       
   *                   sensitive area (sensitiveRadLength) in mm                         
   *  		     
   */
  virtual void addLayer(int nLadders, double phi0,
			double ladderDistance, double ladderOffset, double ladderThickness,
			double ladderLength, double ladderWidth, double ladderRadLength,
			double sensitiveDistance, double sensitiveOffset, double sensitiveThickness,
			double sensitiveLength, double sensitiveWidth, double sensitiveRadLength )
  {
    _layer.addLayer( nLadders, phi0,
		     ladderDistance, ladderOffset, ladderThickness, ladderLength, ladderWidth, ladderRadLength,
		     sensitiveDistance, sensitiveOffset, sensitiveThickness, sensitiveLength, sensitiveWidth, sensitiveRadLength ) ;
    return ;
  }
    
      
  /** Returns the layer layout in the Vertex */
  virtual const VXDLayerLayout & getVXDLayerLayout() const { return _layer ; }
  
  /** The type of Vertex detector: VXDParametersImpl.CCD, VXDParametersImpl.CMOS or
   *  VXDParametersImpl.HYBRID.
   */
  virtual int getVXDType() const { return _vxdType ; }
  
  /** The half length (z) of the support shell in mm.
   */
  virtual double getShellHalfLength() const { return _shellHalfLength ; }

  /** The gap length (z) of the support shell in mm.
   */
  virtual double getShellGap() const { return _shellGap ; }
  
  /** The inner radius of the support shell in mm.
   */
  virtual double getShellInnerRadius() const { return _shellInnerRadius ; }
  
  /** The outer radius of the support shell in mm.
   */
  virtual double getShellOuterRadius() const { return _shellOuterRadius ; }
  
  /** The radiation length in the support shell
   */
  virtual double getShellRadLength() const { return _shellRadLength ; }
  
  /** returns whether a point is inside a ladder
   */
  virtual bool isPointInLadder(Vector3D p) const {
    return isPointInVXD( p, false );
  }
  
  /** returns wheter a point is inside a sensitive volume
   */
  virtual bool isPointInSensitive(Vector3D p) const {
    return isPointInVXD( p, true ) ;
  }
  
  /** returns vector from given point p to nearest ladder
   */
  virtual Vector3D distanceToNearestLadder(Vector3D p) const {
    return distanceToNearestVXD( p, false ) ;
  }

  /** returns vector from given point p to nearest sensitive volume
   */
  virtual Vector3D distanceToNearestSensitive(Vector3D p) const {
    return distanceToNearestVXD( p, true ) ;
  }

  /** returns the first point where a given strainght line
   *  (parameters point p and direction v)  crosses a ladder
   */
  virtual Vector3D intersectionLadder( Vector3D p, Vector3D v ) const {
    return intersectionVXD( p, v, false ) ;
  }

  /** returns the first point where a given strainght line
   *  (parameters point p and direction v)  crosses a sensitive volume
   */
  virtual Vector3D intersectionSensitive( Vector3D p, Vector3D v ) const {
    return intersectionVXD( p, v, true ) ;
  }

protected:
  
  VXDLayerLayoutImpl _layer ;
  
  int _vxdType ;
  
  double _shellInnerRadius ;
  
  double _shellOuterRadius ;

  double _shellHalfLength ;

  double _shellGap ;

  double _shellRadLength ;

private:

  /** returns if a point is in ladder (sensitive == false) or in sensitive (sensitive == true)
   */
  bool isPointInVXD(Vector3D p , bool sensitive = false) const ;

  /** returns distance to nearest ladder (sensitive == false) or nearest sensitiv (sensitive == true)
   */
  Vector3D distanceToNearestVXD(Vector3D p, bool sensitive = false) const ;

  /** returns vector to from point to closest point in described plane
   *  r as spacepoint vector
   *  n as normal vector
   *  u and v as in-plane vectors
   *  minU and maxU as min/max sizes of vector u still in plane
   *  minV and maxV as min/max sizes of vector v still in plane
   */
  Vector3D distanceToPlane(Vector3D p, Vector3D r, Vector3D n, Vector3D u, Vector3D v, float minU, float maxU, float minV, float maxV) const ;

  /** returns the first point on the vxd, where it intersects with a given straight line (parameters point p and direction v)
   */     
  Vector3D intersectionVXD( Vector3D p, Vector3D v, bool sensitive = false) const ;

  /** returns the intersection point of a plane and a straight line
   */
  Vector3D planeLineIntersection( Vector3D r, Vector3D n, Vector3D linePoint, Vector3D lineDir) const ;

  /** returns the confiningRatio of a point p in a plane when the plane (r,n) is reduced to a finite square
   *  extensions in direction u (minU to maxU) and in direction v (minV to maxV)
   */
  double confiningRatio( Vector3D p , Vector3D r, Vector3D n, Vector3D u, Vector3D v, float minU, float maxU, float minV, float maxV ) const ;
  
  /** corrects the vector vPlane into the given borders, given by
   *  direction u (minU to maxU) and in direction v (minV to maxV)
   *  the corrected Vector is returned
   */
  Vector3D correctToBorderPoint( Vector3D vPlane , Vector3D u, Vector3D v, float minU, float maxU, float minV, float maxV ) const ;

  bool isEqual( double valueOne , double valueTwo ) const ;

  bool isEqual( Vector3D p1 , Vector3D p2 ) const ;

  bool differsLess( double valueOne , double valueTwo ) const ;

  
  /** returns Phi corrected to Range -M_PI/2 to +M_PI/2
   */
  double correctPhiRange( double Phi ) const ;

  /** returns Phi for a point
   */
  double getPhiPoint( Vector3D p ) const ;

}; // class

} // namespace gear

#endif /* ifndef GEAR_VXDPARAMETERS_H */
