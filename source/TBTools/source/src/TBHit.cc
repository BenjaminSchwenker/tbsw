// C++ includes
#include <iostream>
#include <cstdlib>

// DEPFETTrackTools includes
#include "TBHit.h"


using namespace lcio;
using namespace std;
using namespace CLHEP;

namespace depfet {

/**
 * File TBHit.cc
 * 
 * Definition of class TBHit.
 * 
 *  @Author B. Schwenker, University of Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */


// Constructors 

TBHit::TBHit() : DAQID(-1), Coord(2, 1), Cov(2,1)
{
  // Set hit coordinates 
  Coord[0][0] = 0;
  Coord[1][0] = 0;   
  // Set hit covariance
  Cov[0][0] = 0;
  Cov[1][1] = 0;   
  Cov[0][1] = 0; 
  // Set hit id 
  UniqueID = -1; 
  // Set raw data pointer
  RawHitPtr = NULL;
  // Set quality 
  Quality = 0; 
}
 
TBHit::TBHit(int newdaqid, double u, double v, double cov_u, double cov_v, double cov_uv, int q) : DAQID(newdaqid), Coord(2, 1), Cov(2,1)
{
  // Set hit coordinates 
  Coord[0][0] = u;
  Coord[1][0] = v;   
  // Set hit covariance
  Cov[0][0] = cov_u;
  Cov[1][1] = cov_v;   
  Cov[0][1] = cov_uv; 
  // Set hit id 
  UniqueID = -1; 
  // Set raw data pointer
  RawHitPtr = NULL;
  // Set quality 
  Quality = q; 
}


TBHit::TBHit(int newdaqid, HepMatrix mcoord, HepSymMatrix mcov) : DAQID(newdaqid), Coord(2, 1), Cov(2,1) 
{
  // Set hit coordinates 
  Coord = mcoord;
  // Set hit covariance
  Cov = mcov;
  // Set hit id 
  UniqueID = -1; 
  // Set raw data pointer
  RawHitPtr = NULL;
  // Set quality 
  Quality = 0; 
}

TBHit::TBHit(lcio::TrackerHit* lciohit) : DAQID(-1), Coord(2, 1), Cov(2,0) 
{
  // Set DAQID
  DAQID =  (int) lciohit->getPosition()[2];
  // Set hit coordinates 
  Coord[0][0] = lciohit->getPosition()[0];
  Coord[1][0] = lciohit->getPosition()[1]; 
  // Set covariance
  Cov[0][0] = lciohit->getCovMatrix ()[0]; 
  Cov[1][1] = lciohit->getCovMatrix ()[1];     
  Cov[0][1] = lciohit->getCovMatrix ()[2];  
  // Set hit id 
  UniqueID = -1; 
  // Set raw data pointer
  RawHitPtr = lciohit;
  // Set quality 
  Quality = lciohit->getType();
}


/** Build LCIO TrackerHit
 * 
 * The lcio::TrackerHit class is used to store reco TBHits persistently in 
 * LCIO files.  
 */
TrackerHitImpl * TBHit::MakeLCIOHit(  ) 
{
  
  TrackerHitImpl * trackerhit = new TrackerHitImpl;
  
  // Set hit position 
  double localPos[3] = {0.,0.,0.};
  localPos[0] = Coord[0][0];
  localPos[1] = Coord[1][0];
  localPos[2] = DAQID;
  trackerhit->setPosition( &localPos[0] );   
       
  // Set covariance matrix as (c_uu, c_vv, c_uv=0)
  float localCov[3] = {0.,0.,0.}; 
  localCov[0] = Cov[0][0]; 
  localCov[1] = Cov[1][1];  
  localCov[2] = Cov[0][1];    
  trackerhit->setCovMatrix ( &localCov[0]);    
  
  // Set quality  
  trackerhit->setType(Quality);
         
  return trackerhit; 
}

HepVector TBHit::GetLocalSpacePoint()
{
  // Hit coordinates are given wrt. w=0 surface of sensor
  // This is the middle thickness surface containing the 
  // centre of the sensor. 
  HepVector Point(3, 0); 
  Point[0] = Coord[0][0]; 
  Point[1] = Coord[1][0];
  Point[2] = 0;
  return Point; 
} 

PixelCluster TBHit::GetCluster() 
{
  
  if (RawHitPtr == NULL ) return PixelCluster();
  
  LCObjectVec clusterVec = RawHitPtr->getRawHits();
  
  if ( clusterVec.size() == 0 )  return PixelCluster();
  
  TrackerData * clusterDigits = dynamic_cast<TrackerData *> ( clusterVec[0] );
    
  return PixelCluster(clusterDigits,DAQID);  
}

StripCluster TBHit::GetStripCluster() 
{
  
  if (RawHitPtr == NULL ) return StripCluster();
  
  LCObjectVec clusterVec = RawHitPtr->getRawHits();
  
  return StripCluster(clusterVec,DAQID,Quality);  
}


} // Namespace
