#ifndef REFERENCEFFRAME_H
#define REFERENCEFFRAME_H 1

#include <Eigen/Core>
#include <Eigen/Geometry>

#include<Eigen/StdVector>
using Eigen::Vector3d;
using Eigen::Matrix3d;
//the following is just to be save. It may slow down compilation time and increase binary size...
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Vector3d)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Matrix3d)

namespace depfet {

/** Class ReferenceFrame
 *   
 *  Handles coordinate transformations between a pair of cartesian coordinate
 *  systems. The coordinate systems will be called GLOBAL and LOCAL frames. 
 *  
 *  A ReferenceFrame object consists of a 3D vector and a 3D rotation matrix. The vector 
 *  is pointing to origin of the local frame. The matrix is the 3x3 rotation 
 *  matrix into the local frame. 
 */


class ReferenceFrame {
  
  Vector3d fPosition;
  Matrix3d fRotation;
  	
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  // Constructors - Identity Transformation 
  ReferenceFrame();

  // Comparison operator
  bool operator==(const ReferenceFrame& other) const; 

  
  ReferenceFrame(Vector3d& position, Matrix3d& rotation);
  	
  // Methods for access to member variables
  Vector3d GetPosition() const { return fPosition; }
  Matrix3d GetRotation() const { return fRotation; }
  void SetPosition(Vector3d position) { fPosition = position; }
  void SetRotation(Matrix3d rotation) { fRotation = rotation; }
  const Vector3d& GetPosRef()  { return fPosition; }
  const Matrix3d& GetRotRef()  { return fRotation; }

  // Method transforming given point from global ref. system to local ref. system
  Vector3d TransformPointToLocal(Vector3d& global) const
    { return fRotation * (global - fPosition); }
  
  // Method transforming given point from local ref. system to global ref. system
  Vector3d TransformPointToGlobal(Vector3d& local) const
    { return fRotation.transpose() * local + fPosition; }
   
  // Method transforming given vector from global ref. system to local ref. system
  Vector3d TransformVecToLocal(Vector3d& globalVec)  const
    { return fRotation * globalVec; }
  
  // Method transforming given vector from global ref. system to local ref. system
  Vector3d TransformVecToGlobal(Vector3d& localVec) const
    { return fRotation.transpose() * localVec;  }
   
  // Get direction cosines of local unit vectors UVW wrt. global XYZ  
  Vector3d GetU() ;
  Vector3d GetV() ;
  Vector3d GetW() ;

  // Get Z position of surface 
  double GetZPosition() const
    { return fPosition[2]; }
  
  // Static member to combine two reference frames. Note that the order is
  // important. Note that also it does *NOT* correspond to two consecutive
  // transformations from global to first frame, and then to the delta
  // frame. Instead, the two Rotations are multiplied and the two positions
  // added.  This is done how Veikko Karimaki defined the "misalignment". It
  // has the advantage that once you are satisfied with the combined values,
  // you can make them to the nominal ones...
  static ReferenceFrame combine_karimaki(ReferenceFrame & first, ReferenceFrame & delta);
  
  // Static member to create a Karimaki delta frame. A delta frame parametrizes
  // a small alignment corrections to sensor position and rotation. A delta frame 
  // may be merged with a 'regular' sensor frame using function combine_karimaki. 
  static ReferenceFrame create_karimaki_delta(double dx, double dy, double dz, double dalpha, double dbeta, double dgamma); 
  
  void PrintHepMatrix() const; 
  void PrintParams() ; 

};  // End class ReferenceFrame
 	

} // Namespace

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(depfet::ReferenceFrame)

#endif // REFERENCEFFRAME_H
