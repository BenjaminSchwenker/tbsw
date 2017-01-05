#ifndef TBHIT_H
#define TBHIT_H 1

// CLHEP includes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>

// Include LCIO header files 
#include "lcio.h"
#include <IMPL/TrackerHitImpl.h>

#include "PixelCluster.h"
#include "StripCluster.h"
	
namespace depfet { 


/** Class TBHit 
 *  
 *  A reco hit in a pixel detector. The following data 
 *  is persistently stored in the LCIO file
 *  
 *  - Measured hit coordinates (u,v) as HepMatrix
 *  - Measurement covariance matrix as HepSymMatrix
 *  - Hit quality flag 
 *  - DAQ ID of detector owning the hit
 * 
 *  In addition, the following info may available 
 *   
 *  - Pointer to cluster raw data (default NULL)
 *  - Unique hitid from pattern reco (default -1)
 *  
 *  @Author B. Schwenker, University of Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class TBHit {
   
 private: // members 
 
  // DAQ ID 
  int DAQID;
  // Local hit coordinates (2x1 matrix)
  CLHEP::HepMatrix Coord;  
  // Hit covariance (2x2 matrix)
  CLHEP::HepSymMatrix Cov;   
  // Unique ID (default -1)
  int UniqueID; 
  // Cluster Quality 
  int Quality; 
  // Pointer to raw TrackerHit (default NULL)
  lcio::TrackerHit * RawHitPtr;
    
 public: // functions 
   
  // Constructors 
  TBHit();
  TBHit(int newdaqid, double u, double v, double cov_u, double cov_v, double cov_uv, int quality);
  TBHit(int newdaqid, CLHEP::HepMatrix coord, CLHEP::HepSymMatrix cov);
  TBHit(lcio::TrackerHit* lciohit);
   
  // Build LCIO TrackerHit
  lcio::TrackerHitImpl * MakeLCIOHit();  
  
  // Get/Set plane number 
  void SetDAQ(int newdaqid) { DAQID = newdaqid; };
  int GetDAQID()  { return DAQID; };
  
  // Dimension of a pixel hit
  int GetDim() { return 2;};
  
  // Get/Set quality of raw cluster
  void SetQuality(int flag) { Quality=flag; }; 
  int GetQuality() { return Quality; };
  
  // Get/Set measured hit coord
  void SetCoord(const CLHEP::HepMatrix& aCoord) { Coord = aCoord; }; 
  CLHEP::HepMatrix&  GetCoord() { return Coord; };
   
  // Get/Set measurment covariance 
  void SetCov(const CLHEP::HepSymMatrix& aCov ) { Cov = aCov; }; 
  CLHEP::HepSymMatrix&  GetCov() { return Cov; };
  
  // Get/Set original raw cluster data   
  void SetRawHit(lcio::TrackerHit* lciohit ) { RawHitPtr = lciohit; }; 
  lcio::TrackerHit* GetRawHit()  { return RawHitPtr; };

  
  PixelCluster GetCluster();

  StripCluster GetStripCluster(); 
  
  // Get/Set unique ID for hit  
  void SetUniqueID(int ID ) { UniqueID = ID; }; 
  int GetUniqueID() { return UniqueID; }; 
   
  // Get 3D point in sensor coordinates 
  CLHEP::HepVector GetLocalSpacePoint(); 
  
};

} // Namespace

#endif // TBHIT_H

