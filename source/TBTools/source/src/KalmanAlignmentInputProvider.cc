//KalmanAlignmentInputProvider.cc
//
// Implementation of KalmanAlignmentInputProvider class.
//
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// C++ includes
//#include <iostream>

// DEPFETTrackTools includes
#include "KalmanAlignmentInputProvider.h"
#include "Utilities.h"

using namespace std;
using namespace CLHEP;

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
  
  // Number of track parameters
  const int nTrackParam = 5;
  
  // Maps daqid of hit to index
  int * ids = new int[nHits];
  
  ////////////////////////////////////////////////////////
  // Create and initialize vectors and matrices to be used 
   
  // This is m, the actual detector measurements
  HepVector Measurement(2*nHits, 0);
        
  // This is V, the measurement covariance matrix
  HepSymMatrix MeasuredCovariance(2*nHits, 0);
  
  // This is p0, the reference track parameters 
  HepVector RefTrackParameters(nTrackParam, 0);
  
  ////////////////////////////////////////////////////////
  // Fill vectors and matrices  
  
  // Reference track parameters 
  HepMatrix refstate = trk.GetReferenceState().GetPars(); 
  RefTrackParameters[0] = refstate[0][0];  
  RefTrackParameters[1] = refstate[1][0];    
  RefTrackParameters[2] = refstate[2][0];          
  RefTrackParameters[3] = refstate[3][0];   
  RefTrackParameters[4] = refstate[4][0];           
  
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
    TBHit& Hit = TE.GetHit();
    for (int k=1; k<=2; ++k) {
      Measurement(index*2+k) = Hit.GetCoord()(k,1); 
      for (int l=k; l<=2; ++l) {     
        MeasuredCovariance(index*2+k,index*2+l) = Hit.GetCov()(k,l);
      }
    }
       	
    ++index;     
  }
      
  ////////////////////////////////////////////////////////
  // Save all vectors/matrices in event
  event.GetIndex()->Adopt(nHits, ids);
  CLHEPtoROOT(Measurement, event.GetMeasurements());
  CLHEPtoROOT(MeasuredCovariance, event.GetMeasuredCovariance());
  CLHEPtoROOT(RefTrackParameters, event.GetRefTrackParameters());
  
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
  TMatrixDSym & MeasuredCovariance = *event.GetMeasuredCovariance();    
   
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
  HepMatrix Pars(5,1,0); 
  Pars[0][0] = RefTrackParameters[0];    // This is u'
  Pars[1][0] = RefTrackParameters[1];    // This is v' 
  Pars[2][0] = RefTrackParameters[2];    // This is u 
  Pars[3][0] = RefTrackParameters[3];    // This is v  
  Pars[4][0] = RefTrackParameters[4];    // This is q/p 

        
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
    HepMatrix mCoord(2, 1);
    mCoord[0][0] = Measurements[ihit*2+0];
    mCoord[1][0] = Measurements[ihit*2+1]; 
    HepSymMatrix mCov(2,1);
    mCov[0][0] = MeasuredCovariance[ihit*2+0][ihit*2+0];
    mCov[1][1] = MeasuredCovariance[ihit*2+1][ihit*2+1]; 
    mCov[0][1] = MeasuredCovariance[ihit*2+0][ihit*2+1];  

    // Add hit to track 
    TBHit Hit(daqid, mCoord, mCov);    
    trk.GetTE(ipl).SetHit(Hit);     
  }
  
  return trk; 
}

} // Namespace;

