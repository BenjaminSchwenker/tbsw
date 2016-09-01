#ifndef REFERENCEFFRAME_H
#define REFERENCEFFRAME_H 1

// Include basic C

// Include CLHEP classes
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/Vector.h>


namespace depfet {

/** Class ReferenceFrame
 *   
 *  Handles coordinate transformations between a pair of cartesian coordinate
 *  systems. The coordinate systems will be called GLOBAL and LOCAL frames. 
 *  
 *  A ReferenceFrame object consists of a HepVector and a HepMatrix. The vector 
 *  is pointing to origin of the local frame. The HepMatrix is the 3x3 rotation 
 *  matrix into the local frame. 
 */
 	
class ReferenceFrame {
  
  CLHEP::HepVector fPosition;
  CLHEP::HepMatrix fRotation;
  	
 public:
  
  // Constructors - Identity Transformation 
  ReferenceFrame();

  // Comparison operator
  bool operator==(const ReferenceFrame& other) const; 

  
  ReferenceFrame(CLHEP::HepVector& position, CLHEP::HepMatrix& rotation);
  	
  // Methods for access to member variables
  CLHEP::HepVector GetPosition() const { return fPosition; };
  CLHEP::HepMatrix GetRotation() const { return fRotation; };
  void SetPosition(CLHEP::HepVector position) { fPosition = position; };
  void SetRotation(CLHEP::HepMatrix rotation) { fRotation = rotation; };
  CLHEP::HepVector& GetPosRef()  { return fPosition; };
  CLHEP::HepMatrix& GetRotRef()  { return fRotation; };
  
  // Method transforming given point from global ref. system to local ref. system
  CLHEP::HepVector TransformPointToLocal(CLHEP::HepVector& global) const  
    { return fRotation * (global - fPosition); };
  
  // Method transforming given point from local ref. system to global ref. system
  CLHEP::HepVector TransformPointToGlobal(CLHEP::HepVector& local) const 
    { return fRotation.T() * local + fPosition; };
   
  // Method transforming given vector from global ref. system to local ref. system
  CLHEP::HepVector TransformVecToLocal(CLHEP::HepVector& globalVec)  const
    { return fRotation * globalVec;   }
  
  // Method transforming given vector from global ref. system to local ref. system
  CLHEP::HepVector TransformVecToGlobal(CLHEP::HepVector& localVec) const
    { return fRotation.T() * localVec;  }
   
  // Get direction cosines of local unit vectors UVW wrt. global XYZ  
  CLHEP::HepVector GetU() ;
  CLHEP::HepVector GetV() ;
  CLHEP::HepVector GetW() ;

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
  
  void PrintHepMatrix() ; 
  void PrintParams() ; 

};  // End class ReferenceFrame
 	

} // Namespace

#endif // REFERENCEFFRAME_H
