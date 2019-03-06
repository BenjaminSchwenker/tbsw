#ifndef HITFACTORY_H
#define HITFACTORY_H 1

// Include TrackTools header files
#include "TBHit.h"
#include "TBDetector.h"

// Include STL
#include <vector>
#include <map>

namespace depfet { 
 
  

/** Class HitFactory 
 *  
 *  The HitFactory class is a container for reconstructed hits in a
 *  triggered event. Reco hits are sorted by their sub detector ID. 
 *  The sub detectors are numbered starting from 0 counting in beam 
 *  direction. A unique Id is assigned to each hit. At the moment, 
 *  the hit Id is defined by the order of insertions into the HitFactory.
 *     
 *  The HitFactory class is able to preselect all hits compatible to a
 *  given particle intersection point. The preselection is based on a 
 *  subdivision of the sensitive detector area in strip-like sectors. 
 *  The sector pitch is typically aroung 1mm.
 *  
 *  Author: Benjamin Schwenker, Göttingen University 
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
  	
class HitFactory {
  
 typedef std::map< int , std::vector<int> > SectorMap;
  
 public:
   
  // Constructor
  HitFactory(const TBDetector& Detector, double SectorPitch=1); 
  
  // Add hit to pattern factory 
  void AddRecoHit(const TBHit& hit); 
  
  // Get all hits for plane ipl 
  std::vector<TBHit>& GetHits(int ipl);
  
  // Get reco hit at position ihit from sensor ipl
  const TBHit& GetRecoHitFromID(int ihit, int ipl) const; 
  
  // Get Ids of hits compatible with track at (u,v) 
  std::vector<int> GetCompatibleHitIds(int ipl, double u, double distMaxU) const;
  
  // Get total number of hits  
  int GetNHits() const;
   
  // Get number of hits for detector ipl 
  int GetNHits(int ipl) const; 

 private: 
  
  // Reference to TBDetector
  const TBDetector& m_detector;
  
  // Subdivide detector along u-axis 
  double m_sectorPitch; 
  
  // Vector of hits for each detector 
  std::vector< std::vector<TBHit> > m_hitStore; 
    
  // Vector of sector maps for each detector 
  std::vector< SectorMap > m_hitSubStore; 
  
};

} // Namespace

#endif // HITFACTORY_H

