#ifndef TBVertexFitter_H
#define TBVertexFitter_H 1

// DEPFETTrackTools includes
#include "TBTrack.h"
#include "TBVertex.h"
#include "TBDetector.h"

//other includes includes
#include <TMath.h>

// CLHEP includes 
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>

namespace depfet {

 
//! Class TBVertexFitter
/* 
 * The TBVertexFitter class performs the fitting of a scattering vertex
 * using trackstates of charged particles in a tracking telescope at 
 * the plane of the device under testing.
 * Inside this class the vertex position r is given in (x,y,z) for readability,
 * however the vertex position is always in local coordinates (u,v,w) of the DUT.
 *
 * The function FitVertex uses a Kalman Filter to reconstruct the
 * Vertex position using the trackstates at the DUT.
 * The relevant measurement equation h used in the Kalman filter reads
 * h = ( a, b, x + a*z, y + b*z)
 */

class TBVertexFitter {
  
 public:
  
  /** Default Constructor  
   */
  TBVertexFitter(int ipl, TBDetector det);

  /** Performs vertex fitting. Returns boolean of success of fit
   */
  bool FitVertex(TBVertex& Vertex);
  
 // Public Members -----------------

  CLHEP::HepMatrix B; //Jacobian B = dh/dq for slope states q = (a,b)

  CLHEP::HepMatrix A; //Jacobian A = dh/dr for vertex state r = (x,y,z)  

 // Private Members -----------------
 private:
  
  int plnr; //plane number of dut
  TBDetector detector; // Detector object

};
 

} // Namespace

#endif 
