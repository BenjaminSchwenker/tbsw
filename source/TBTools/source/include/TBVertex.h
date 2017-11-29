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

  //Vertex Position vector (in local and global coordinates)
  CLHEP::HepMatrix Pos;
  CLHEP::HepMatrix GlobalPos;
  //Corresponding covariance (in local and global coordinates)
  CLHEP::HepMatrix Cov;
  CLHEP::HepMatrix GlobalCov;
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
  TBVertex(CLHEP::HepMatrix aPos, CLHEP::HepMatrix aGlobalPos, CLHEP::HepMatrix aCov, CLHEP::HepMatrix aGlobalCov, double achi2);  
  
  // Get/Set local vertex position  
  void SetPos(const CLHEP::HepMatrix& aPos) { Pos= aPos; }; 
  CLHEP::HepMatrix&  GetPos() { return Pos; };

  // Get/Set global vertex position  
  void SetGlobalPos(const CLHEP::HepMatrix& aPos) { GlobalPos= aPos; }; 
  CLHEP::HepMatrix&  GetGlobalPos() { return GlobalPos; };
   
  // Get/Set local position covariance 
  void SetCov(const CLHEP::HepMatrix& aCov ) { Cov = aCov; }; 
  CLHEP::HepMatrix&  GetCov() { return Cov; };

  // Get/Set global position covariance 
  void SetGlobalCov(const CLHEP::HepMatrix& aCov ) { GlobalCov = aCov; }; 
  CLHEP::HepMatrix&  GetGlobalCov() { return GlobalCov; };

  // Get/Set chi²-value
  void SetChi2(double achi2 ) { chi2 = achi2; }; 
  double GetChi2() { return chi2; };
  double GetChi2Ndof() { return chi2/double(ndf); };

  // Get/Set filter residual
  void SetRes(const CLHEP::HepMatrix& aRes ) { Res = aRes; }; 
  CLHEP::HepMatrix& GetRes() { return Res; };

  // Get/Set ndf-value
  void SetNdf(int andf ) { ndf = andf; }; 
  int GetNdf() { return ndf; };

  // add state
  void AddTrackState(TBTrackState& state) { States.push_back(state); };

  // remove last state
  void RemoveLastTrackState() { States.pop_back(); };

  // get States vector
  std::vector<TBTrackState>& GetStates() { return States; };

};

} //Namespace

#endif
