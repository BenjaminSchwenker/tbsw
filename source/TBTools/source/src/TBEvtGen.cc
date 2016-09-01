// TBEvtGen implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "TBEvtGen.h"

// C++ includes
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <limits>


using namespace std;
using namespace CLHEP;

namespace depfet {

// Constructor
TBEvtGen::TBEvtGen( ) 
{
  // Provide default values for all 
  // settings 
 
  // Particles per sec  
  intensity = 10000; 
  // M26 integration time 100us
  readouttime = 0.0001; 
  // Trigger on first and last hit
  trgmode = 2;  

  // Start at t=0 sec 
  time = 0; 
  // Last trigger in sec
  lasttrg = -numeric_limits< double >::max();
  
  myRng = new TRandom3();
  
}
         
// Destructor
TBEvtGen::~TBEvtGen( ) {
  delete  myRng; 
} 

// Set random number generator seed
void TBEvtGen::SetSeed(int seed) {
  myRng->SetSeed(seed);
}

// ACCEPT
bool TBEvtGen::ACCEPT(TBTrack& TruthTrack) {
  
  // Increment simulation time with each track
  time += 1.0/intensity; 
  
  // All tracks passing telescope within readout  
  // time get accepted
  if ( time < lasttrg + readouttime ) {   
    return true; 
  }
  
  // Otherwise, track is required to generate a new
  // trigger coincidence
  if ( HasCoincidence(TruthTrack) ) {
    lasttrg = time; 
    return true;   
  } 
  
  return false;  
}
  
// HasCoincidence 
bool  TBEvtGen::HasCoincidence(TBTrack& TruthTrack) {
  
  // Return value
  bool val = false; 
  
  // Number of detectors
  int nSensor = TruthTrack.GetNumTEs(); 
  
  if ( trgmode == 0 ) {
     val = HasHit(TruthTrack, 0); 
  } else if ( trgmode == 1 ) {
     val = HasHit(TruthTrack, nSensor-1); 
  } else { 
     val = HasHit(TruthTrack, 0) &&  HasHit(TruthTrack, nSensor-1); 
  }
  
  return val;   
}
 
/** HasHit returns true if layer ipl fires
 */
bool  TBEvtGen::HasHit(TBTrack& TruthTrack, int ipl) {

  // Check sensor fires 
  if ( !TruthTrack.GetTE(ipl).IsCrossed() ) return false; 
             
  // Get truth track intersect 
  double u = TruthTrack.GetTE(ipl).GetState().GetPars()[2][0]; 
  double v = TruthTrack.GetTE(ipl).GetState().GetPars()[3][0]; 
  
  double CenterU = 0; 
  double SizeU = TruthTrack.GetTE(ipl).GetDet().GetSensitiveSizeU();
     
  double CenterV = 0; 
  double SizeV = TruthTrack.GetTE(ipl).GetDet().GetSensitiveSizeV();
  
  if (u < CenterU-SizeU/2.  || u > CenterU+SizeU/2.) {
   return false;
  }
  if (v < CenterV-SizeV/2. || v > CenterV+SizeV/2.) {
    return false;
  }
  
  return true; 
}
 
} // Namespace;
