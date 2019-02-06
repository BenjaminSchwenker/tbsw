#ifndef HITFACTORY_H
#define HITFACTORY_H 1

// Include TrackTools header files
#include "TBHit.h"
#include "TBDetector.h"

// Include STL
#include <vector>
#include <map>

namespace depfet { 
 
  typedef std::map< int , std::vector<int> > SectorMap;

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
  
 private: 
  
  // Total number of pixel sensors 
  TBDetector _detector;
  
  // Subdivide detector along u-axis 
  double _SectorPitch; 
  
  // Vector of hits for each detector 
  std::vector< std::vector<TBHit> > _HitStore; 
    
  // Vector of sector maps for each detector 
  std::vector< SectorMap > _HitSubStore; 
  
 public:
   
  // Constructor
  HitFactory(TBDetector& Detector, double SectorPitch=1); 
  
  // Add hit to pattern factory 
  void AddRecoHit(TBHit& hit); 
  
  // Get all hits for plane ipl 
  std::vector<TBHit>& GetHits(int ipl);
  
  // Get reco hit at position ihit from sensor ipl
  TBHit& GetRecoHitFromID(int ihit, int ipl); 
  
  // Get Ids of hits compatible with track at (u,v) 
  std::vector<int> GetCompatibleHitIds(int ipl, double u, double v, double distMaxU, double distMaxV=0);
  
  // Get total number of hits  
  int GetNHits();
   
  // Get number of hits for detector ipl 
  int GetNHits(int ipl); 
  
};

} // Namespace

#endif // HITFACTORY_H

