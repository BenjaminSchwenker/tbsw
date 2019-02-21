// C++ includes
#include <iostream>
#include <cstdlib>

// TBTools includes
#include "TBHit.h"


using namespace lcio;
using namespace std;

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

TBHit::TBHit() : DAQID(-1)
{
  // Set hit coordinates 
  Coord<<0, 0;
  // Set hit covariance
  Cov<<0, 0,
       0, 0;

  // Set hit id 
  UniqueID = -1; 
  // Set raw data pointer
  RawHitPtr = nullptr;
  // Set quality 
  Quality = 0; 
}
 
TBHit::TBHit(int newdaqid, double u, double v, double cov_u, double cov_v, double cov_uv, int q) : DAQID(newdaqid)
{
  // Set hit coordinates 
  Coord<< u, v;
  // Set hit covariance
  Cov<<   cov_u,    cov_uv,
          cov_uv,   cov_v;

  // Set hit id 
  UniqueID = -1; 
  // Set raw data pointer
  RawHitPtr = nullptr;
  // Set quality 
  Quality = q; 
}


TBHit::TBHit(int newdaqid, const Vector2d &mcoord, const Matrix2d &mcov) : DAQID(newdaqid)
{
  // Set hit coordinates 
  Coord = mcoord;
  // Set hit covariance
  Cov = mcov;
  // Set hit id 
  UniqueID = -1; 
  // Set raw data pointer
  RawHitPtr = nullptr;
  // Set quality 
  Quality = 0; 
}

TBHit::TBHit(lcio::TrackerHit* lciohit) : DAQID(-1) 
{
  // Set DAQID
  DAQID =  (int) lciohit->getPosition()[2];
  // Set hit coordinates 
  Coord<< lciohit->getPosition()[0], lciohit->getPosition()[1];
  // Set covariance
  const auto & cov_m=lciohit->getCovMatrix();

  Cov << cov_m[0], cov_m[2],
         cov_m[2], cov_m[1];
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
  localPos[0] = Coord[0];
  localPos[1] = Coord[1];
  localPos[2] = DAQID;
  trackerhit->setPosition( &localPos[0] );   
       
  // Set covariance matrix as (c_uu, c_vv, c_uv=0)
  float localCov[3] = {0.,0.,0.}; 
  localCov[0] = Cov(0,0);
  localCov[1] = Cov(1,1);
  localCov[2] = Cov(0,1);
  trackerhit->setCovMatrix ( &localCov[0]);    
  
  // Set quality  
  trackerhit->setType(Quality);
         
  return trackerhit; 
}

Vector3d TBHit::GetLocalSpacePoint()
{
  // Hit coordinates are given wrt. w=0 surface of sensor
  // This is the middle thickness surface containing the 
  // centre of the sensor. 
  Vector3d Point;
  Point[0] = Coord[0];
  Point[1] = Coord[1];
  Point[2] = 0;
  return Point; 
} 

PixelCluster TBHit::GetCluster() 
{
  
  if (RawHitPtr == nullptr ) return PixelCluster();
  
  LCObjectVec clusterVec = RawHitPtr->getRawHits();
  
  if ( clusterVec.size() == 0 )  return PixelCluster();
  
  TrackerData * clusterDigits = dynamic_cast<TrackerData *> ( clusterVec[0] );
    
  return PixelCluster(clusterDigits,DAQID);  
}

StripCluster TBHit::GetStripCluster() 
{
  
  if (RawHitPtr == nullptr ) return StripCluster();
  
  LCObjectVec clusterVec = RawHitPtr->getRawHits();
  
  return StripCluster(clusterVec,DAQID,Quality);  
}


} // Namespace
