// C++ includes
#include <iostream>
#include <cassert>
#include <cstdlib>

// Local includes
#include "HitFactory.h"

using namespace std;

namespace depfet {

/**
 * File HitFactory.cc
 * 
 * Definition of class HitFactory
 */


// Constructor - Hit Factory for nSensors
HitFactory::HitFactory(TBDetector& Detector, double SectorPitch)  
{

  // Deep copy of detector
  _detector = Detector;
  
  // Size of sectors in all sector maps
  _SectorPitch = SectorPitch;
  
  // Manage hits for n subdectors       
  _HitStore.resize(_detector.GetNSensors());
  _HitSubStore.resize(_detector.GetNSensors());  
}
 

// Add hit to factory 
void HitFactory::AddRecoHit(TBHit & hit)
{
  // To which detector does hit belong 
  int daqid = hit.GetDAQID();
  int ipl = _detector.GetPlaneNumber(daqid);  

  // Check that daqid is registered in the gear file
  if ( ipl < 0 ) return;  
  
  // Compute unique id of this hit 
  int hitID = (int) _HitStore[ipl].size(); 
  
  // Compute unique sector id
  int secID = floor( hit.GetCoord()[0][0] / _SectorPitch );
  
  // Now we can proceed to register the hit 
  _HitStore[ipl].push_back(hit);     
  
  // Create a new sector if necessary 
  if ( _HitSubStore[ipl].find(secID) == _HitSubStore[ipl].end() ) {
    _HitSubStore[ipl][secID] = std::vector<int>();
  }
  
  _HitSubStore[ipl][secID].push_back(hitID);
}

// Get all hits for plane ipl 
vector<TBHit>& HitFactory::GetHits(int ipl)
{
  // Check if ipl is valid!!
  return _HitStore[ipl];  
}

// Get reco hit at position ihit from sensor ipl
TBHit& HitFactory::GetRecoHitFromID(int ihit, int ipl)
{
  // Check if ipl is valid!!
  // Check if ihit is valid!!
  return _HitStore[ipl][ihit];  
}

// Get number of hits for sensor ipl 
int HitFactory::GetNHits(int ipl) 
{
  // Check if ipl is valid!! 
  return _HitStore[ipl].size();
}

// Get total number of hits  
int  HitFactory::GetNHits( )
{
  int TotHits = 0; 
  // Loop over detectors  
  for (int ipl=0; ipl<(int)_HitStore.size(); ++ipl) { 
    TotHits += (int)_HitStore[ipl].size();  
  }  
  return TotHits;
} 

/** Get Ids of hits compatible with track at (u,v)  
 *  
 *  A hit with coordinates (um,vm) is compatible iff lying in 
 *  a compatible sector. 
 */  
vector<int> HitFactory::GetCompatibleHitIds(int ipl, double u, double v, double distMaxU, double distMaxV)
{
  
  // List of id for compatible hits 
  vector<int> CompatibleHitIds;
  
  // Compute range of compatible sectors 
  int minSecID = floor( (u-distMaxU) / _SectorPitch );
  int maxSecID = floor( (u+distMaxU) / _SectorPitch );
  
  // Add all hits in all compatible sectors 
  for (int secID = minSecID; secID <= maxSecID; ++secID) {
    
    if ( _HitSubStore[ipl].find(secID) == _HitSubStore[ipl].end() ) {
      continue;   
    } else {
      int nHits = _HitSubStore[ipl][secID].size();
      for (int hitID=0; hitID<nHits; ++hitID) { 
        CompatibleHitIds.push_back( _HitSubStore[ipl][secID][hitID] );  
      }
    }  
    
  }
    
  // Return compatible hits
  return CompatibleHitIds;
}

} // Namespace;
