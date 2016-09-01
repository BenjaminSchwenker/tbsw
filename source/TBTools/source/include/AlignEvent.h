#ifndef _AlignEvent_h
#define _AlignEvent_h
 	
#include "TObject.h"
#include "TClonesArray.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TArrayI.h"

namespace depfet {
 

/** Class AlignEvent
 *  
 *  The AlignEvent class provides root persistency for TBTrack objects. 
 *  This is needed for track based alignment algorithms.  
 *  
 *  The track data is written to and read from AlignEvent objects using 
 *  the class KalmanAlignmentInputProvider. 
 *  
 *  The alignment program basically needs the following data: 
 *  - track quality indicators 
 *  - list of hit detectors  
 *  - measured hit coordinates (+cov)
 *  - mass, charge, momentum
 *  - track parameters at reference plane
 *  
 *  Author: Benjamin Schwenker, Göttingen University 
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
 	
class AlignEvent {
  
 protected:
  double fMass;      // Mass (GeV) 
  double fCharge;    // Charge (e) 
  double fMomentum;  // Momentum (GeV) 
  double fChi2;      // chi^2 of track fit
  TArrayI     * fIndex;                // daqid of hit planes
  TVectorD    * fMeasurements;         // measured hits
  TMatrixDSym * fMeasuredCovariance;   // measured covariance matrix
  TVectorD    * fRefTrackParameters;   // reference parameters 
  int fRefPlane;                       // plane number of reference frame

 public:
  AlignEvent();
  virtual ~AlignEvent();
  void SetMass(double mass) { fMass = mass; };
  void SetCharge(double charge) { fCharge = charge; };
  void SetMomentum(double momentum) { fMomentum = momentum; };
  void SetRefPlaneNumber(int plane) { fRefPlane = plane; };
  void SetChi2(double chi2) { fChi2 = chi2; };
  double GetChi2() { return fChi2; };
  double GetMass() { return fMass; };
  double GetCharge() { return fCharge; };
  double GetMomentum() { return fMomentum; };
  int GetRefPlaneNumber() { return fRefPlane; };
     
  // non-const getters for index, vectors, matrices to allow setting
  // without constructing
  TArrayI     * GetIndex() { return fIndex; };
  TVectorD    * GetMeasurements() { return fMeasurements; };
  TMatrixDSym * GetMeasuredCovariance() { return fMeasuredCovariance; };
  TVectorD    * GetRefTrackParameters() { return fRefTrackParameters; };
  
  
  ClassDef(AlignEvent, 1); // Alignment event
};

} // Namespace
#endif

