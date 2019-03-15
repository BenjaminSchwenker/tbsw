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
HitFactory::HitFactory(const TBDetector& detector, double SectorPitch)  : m_detector(detector)
{

  // Size of sectors in all sector maps
  m_sectorPitch = SectorPitch;
  
  // Manage hits for n subdectors       
  m_hitStore.resize(m_detector.GetNSensors());
  m_hitSubStore.resize(m_detector.GetNSensors());  
}
 

// Add hit to factory 
void HitFactory::AddRecoHit(const TBHit & hit)
{
  // To which detector does hit belong 
  int sensorid = hit.GetSensorID();
  int ipl = m_detector.GetPlaneNumber(sensorid);  

  // Check that plane number is valid  
  if ( ipl < 0 ) return;  
  
  // Compute unique id of this hit 
  int hitID = (int) m_hitStore[ipl].size(); 
  
  // Compute unique sector id
  int secID = floor( hit.GetCoord()[0] / m_sectorPitch );
  
  // Now we can proceed to register the hit 
  m_hitStore[ipl].push_back(hit);     
  
  // Create a new sector if necessary 
  if ( m_hitSubStore[ipl].find(secID) == m_hitSubStore[ipl].end() ) {
    m_hitSubStore[ipl][secID] = std::vector<int>();
  }
  
  m_hitSubStore[ipl][secID].push_back(hitID);
}

// Get all hits for plane ipl 
vector<TBHit>& HitFactory::GetHits(int ipl)
{
  // Check if ipl is valid!!
  return m_hitStore[ipl];  
}

// Get reco hit at position ihit from sensor ipl
const TBHit& HitFactory::GetRecoHitFromID(int ihit, int ipl) const
{
  // Check if ipl is valid!!
  // Check if ihit is valid!!
  return m_hitStore[ipl][ihit];  
}

// Get number of hits for sensor ipl 
int HitFactory::GetNHits(int ipl) const 
{
  // Check if ipl is valid!! 
  return m_hitStore[ipl].size();
}

// Get total number of hits  
int  HitFactory::GetNHits( ) const 
{
  int TotHits = 0; 
  // Loop over detectors  
  for (int ipl=0; ipl<(int)m_hitStore.size(); ++ipl) { 
    TotHits += (int)m_hitStore[ipl].size();  
  }  
  return TotHits;
} 

/** Get Ids of hits compatible with track at (u,v)  
 *  
 *  A hit with coordinates (um,vm) is compatible iff lying in 
 *  a compatible sector. 
 */  
vector<int> HitFactory::GetCompatibleHitIds(int ipl, double u, double distMaxU) const
{
  
  // List of id for compatible hits 
  vector<int> CompatibleHitIds;
  CompatibleHitIds.reserve(8);
  
  // Compute range of compatible sectors 
  int minSecID = floor( (u-distMaxU) / m_sectorPitch );
  int maxSecID = floor( (u+distMaxU) / m_sectorPitch );
  
  // Add all hits in all compatible sectors 
  for (int secID = minSecID; secID <= maxSecID; ++secID) {
    
    if ( m_hitSubStore[ipl].find(secID) == m_hitSubStore[ipl].end() ) {
      continue;   
    } else {
      int nHits = m_hitSubStore[ipl].at(secID).size();
      for (int hitID=0; hitID<nHits; ++hitID) { 
        CompatibleHitIds.push_back( m_hitSubStore[ipl].at(secID)[hitID] );  
      }
    }  
    
  }
    
  // Return compatible hits
  return CompatibleHitIds;
}

} // Namespace;
