// TBTrack implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// C++ includes
//#include <iostream>
//#include <cassert>
//#include <cstdlib>

// TBTools include
#include "TBTrack.h"

// Marlin includes 
//#include <marlin/Global.h>
//#include <streamlog/streamlog.h>

// Namespaces
//using namespace marlin;
using namespace std; 

namespace depfet {

/** Constructor    
 */
TBTrack::TBTrack(TBDetector& detector)
{
  // The skeleton track model is defined by the 
  // detector setup
  int nDet = detector.GetNSensors(); 
  TEVec.reserve(nDet);
  for(int ipl=0;ipl<nDet;ipl++)  {  
    // Create track element 
    // Add track element  
    TEVec.push_back(detector.GetDet(ipl));
  }
  
  // Particle hypothesis  
  Mass = 0;
  Charge = 0;
  Mom = 0; 
  // No track fit done yet   
  ChiSqu = -1;
  Ndof = 0;
} 


/** Get track element by plane number 
 */  
TBTrackElement& TBTrack::GetTE(int ipl)
{
  return TEVec.at(ipl);  
}

/** Get track element by plane number 
 */  
const TBTrackElement& TBTrack::GetTE(int ipl) const
{
  return TEVec.at(ipl);  
}


/** Get all track elements  
 */  
std::vector<TBTrackElement>& TBTrack::GetTEs()
{
  return TEVec;
}

/** Get all track elements  
 */  
const std::vector<TBTrackElement>& TBTrack::GetTEs() const
{
  return TEVec;
}

/** Get number of hits in track   
 */  
int TBTrack::GetNumHits() const
{
  // Get all track elements (TEs) 
  int nTE = GetNumTEs();
  int NumHits = 0; 
  // Count hits in track 
  for(int iTE=0;iTE<nTE;++iTE) {
    if ( TEVec[iTE].HasHit() ) ++NumHits;  
  }   
  return NumHits; 
}

/** Get number of track elements    
 */  
int TBTrack::GetNumTEs() const
{
  return (int) TEVec.size();
}

    
} // Namespace;

