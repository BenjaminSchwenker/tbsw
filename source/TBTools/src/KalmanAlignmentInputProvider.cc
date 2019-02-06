//KalmanAlignmentInputProvider.cc
//
// Implementation of KalmanAlignmentInputProvider class.
//
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// C++ includes
//#include <iostream>

// TBTools includes
#include "KalmanAlignmentInputProvider.h"
#include "Utilities.h"

using namespace std;

namespace depfet {


// Fill event structure
void KalmanAlignmentInputProvider::FillEvent(TBTrack& trk, AlignEvent& event)
{
   
  // Copy particle hypothesis
  event.SetMass(trk.GetMass());  
  event.SetCharge(trk.GetCharge());
  event.SetMomentum(trk.GetMomentum());
  
  // Copy quality indicators
  event.SetChi2(trk.GetChiSqu()); 

  // Reference plane number 
  event.SetRefPlaneNumber(trk.GetReferenceState().GetPlane());  
  
  // We assume two measured coord per detector. The track
  // fit needs track parameters at z=0 as seed. 
  
  // Number of detectors in track
  const int nDet = trk.GetNumTEs();
  
  // Track has nHits hits  
  const int nHits = trk.GetNumHits();
  
  // Maps daqid of hit to index
  int * ids = new int[nHits];
  
  ////////////////////////////////////////////////////////
  // Create and initialize vectors and matrices to be used 
   
  TVectorD * newMeasurement =  event.GetMeasurements();
  newMeasurement->ResizeTo(2*nHits);
        
  TMatrixD * newMeasuredCovariance = event.GetMeasuredCovariance();
  newMeasuredCovariance->ResizeTo(2*nHits, 2);
  
  TVectorD * newRefTrackParameters = event.GetRefTrackParameters();
  newRefTrackParameters->ResizeTo(5);   
  
  ////////////////////////////////////////////////////////
  // Fill vectors and matrices  
  
  // Reference track parameters 
  auto refstate = trk.GetReferenceState().GetPars(); 
  (*newRefTrackParameters)[0] = refstate(0);     
  (*newRefTrackParameters)[1] = refstate(1);  
  (*newRefTrackParameters)[2] = refstate(2);  
  (*newRefTrackParameters)[3] = refstate(3);  
  (*newRefTrackParameters)[4] = refstate(4);    
  
  // Fill index, measurments and measured cov      
  int index = 0;  
       
  // Loop over all track elements 
  for (int ipl= 0; ipl<nDet; ++ipl) {  
         
    TBTrackElement& TE = trk.GetTE(ipl);  
         
    // Skip track elements w/o measurment
    if ( !TE.HasHit() ) continue;
            
    // Ok, register hit on this detector  
    ids[index] = TE.GetDet().GetDAQID(); 
    
    // Fill measurments + cov 
    auto hitCoord = TE.GetHit().GetCoord();
    auto hitCov = TE.GetHit().GetCov();
    for (int k=0; k<=1; ++k) {
      (*newMeasurement)[index*2+k] = hitCoord(k); 
      for (int l=0; l<=1; ++l) {     
        (*newMeasuredCovariance)[index*2+k][l] = hitCov(k,l);
      }
    }
       	
    ++index;     
  }
      
  ////////////////////////////////////////////////////////
  // Save all vectors/matrices in event
  event.GetIndex()->Adopt(nHits, ids);

  return; 
}

// Create TBTrack structure
TBTrack KalmanAlignmentInputProvider::MakeTBTrack( AlignEvent& event, TBDetector& detector )
{
  
  ////////////////////////////////////////////////////////
  // Read vectors/matrices from event
  
  TVectorD & RefTrackParameters = *event.GetRefTrackParameters();     
  TArrayI & ids = *event.GetIndex();
  TVectorD & Measurements = *event.GetMeasurements();   
  TMatrixD & MeasuredCovariance = *event.GetMeasuredCovariance();    
   
  ////////////////////////////////////////////////////////
  // Create and fill TBTrack object  
  
  // Create track  
  TBTrack trk(detector);

  // Set particle hypothesis   
  trk.SetMass(event.GetMass());
  trk.SetCharge(event.GetCharge());
  trk.SetMomentum(event.GetMomentum());  
  
  // Set track chisqu
  trk.SetChiSqu(event.GetChi2());

  // Create reference track state (seed)
  TrackState Pars; 
  Pars(0) = RefTrackParameters[0];    // This is u'
  Pars(1) = RefTrackParameters[1];    // This is v' 
  Pars(2) = RefTrackParameters[2];    // This is u 
  Pars(3) = RefTrackParameters[3];    // This is v  
  Pars(4) = RefTrackParameters[4];    // This is q/p 

        
  TBTrackState Seed;
  Seed.Pars = Pars;
  Seed.SetPlane(event.GetRefPlaneNumber()); 
  trk.SetReferenceState(Seed);
                             
  // Add hits to track     
  int nHit = ids.GetSize();
  
  for (int ihit=0; ihit<nHit; ihit++) {
    // Create a TBHit 
    int daqid = ids[ihit]; 
    int ipl = detector.GetPlaneNumber(daqid); 
    Vector2d mCoord;
    mCoord << Measurements[ihit*2+0] , Measurements[ihit*2+1]; 
    Matrix2d mCov;
    mCov(0,0) = MeasuredCovariance[ihit*2+0][0];
    mCov(0,1) = MeasuredCovariance[ihit*2+0][1];  
    mCov(1,0) = MeasuredCovariance[ihit*2+1][0]; 
    mCov(1,1) = MeasuredCovariance[ihit*2+1][1];  
    
    // Add hit to track 
    TBHit Hit(daqid, mCoord, mCov);    
    trk.GetTE(ipl).SetHit(Hit);     
  }
  
  return trk; 
}

} // Namespace;

