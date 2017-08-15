#ifndef TBVERTEX_H
#define TBVERTEX_H 1

// DEPFETTrackTools includes
#include "TBTrack.h"

// CLHEP includes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>

namespace depfet { 


/** Class TBVertex
 *  
 *  The class TBVertex represents a vertex position
 *  relative to the scattering plane used for the 
 *  calculation of said vertex.
 *  It consists of the position r = (x,y,z), the
 *  corresponding covariance matrix C and 
 *  the Chi²-Value of the vertex fit.
 *  
 *  @Author E. Heinrich, University of Göttingen
 *  <mailto:erik.heinrich@stud.uni-goettingen.de>
 */

class TBVertex {

 public: //Members

  //Vertex Position vector
  CLHEP::HepMatrix Pos;
  //Corresponding covariance
  CLHEP::HepMatrix Cov;
  //Vertex fit chi²
  double chi2;
  //Vertex fit number degrees of freedom
  int ndf;
  //filter residual
  CLHEP::HepMatrix Res;

  //Vector containing all states for given measurement
  std::vector<TBTrackState> States;

 public: //Functions

  // Constructors
  TBVertex(); 
  TBVertex(CLHEP::HepMatrix aPos, CLHEP::HepMatrix aCov, double achi2);  
  
  // Get/Set vertex position  
  void SetPos(const CLHEP::HepMatrix& aPos) { Pos= aPos; }; 
  CLHEP::HepMatrix&  GetPos() { return Pos; };
   
  // Get/Set position covariance 
  void SetCov(const CLHEP::HepMatrix& aCov ) { Cov = aCov; }; 
  CLHEP::HepMatrix&  GetCov() { return Cov; };

  // Get/Set chi²-value
  void SetChi2(double achi2 ) { chi2 = achi2; }; 
  double GetChi2() { return chi2; };

  // Get/Set filter residual
  void SetRes(const CLHEP::HepMatrix& aRes ) { Res = aRes; }; 
  CLHEP::HepMatrix& GetRes() { return Res; };

  // Get/Set ndf-value
  void SetNdf(int andf ) { ndf = andf; }; 
  int GetNdf() { return ndf; };

  // add state
  void AddTrackState(TBTrackState& state) { States.push_back(state); };

  // get States vector
  std::vector<TBTrackState>& GetStates() { return States; };

};

} //Namespace

#endif
