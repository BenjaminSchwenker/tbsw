#ifndef TBHIT_H
#define TBHIT_H 1

#include <Eigen/Core>

// Include LCIO header files 
#include "lcio.h"
#include <IMPL/TrackerHitImpl.h>

#include "PixelCluster.h"
#include "StripCluster.h"
using Eigen::Matrix2d;
using Eigen::Vector2d;
using Eigen::Vector3d;
namespace depfet { 


/** Class TBHit 
 *  
 *  A reco hit in a pixel detector. The following data 
 *  is persistently stored in the LCIO file
 *  
 *  - Measured hit coordinates (u,v)
 *  - Measurement covariance matrix 
 *  - Hit quality flag 
 *  - DAQ ID of detector owning the hit
 * 
 *  In addition, the following info may available 
 *   
 *  - Pointer to cluster raw data (default nullptr)
 *  - Unique hitid from pattern reco (default -1)
 *  
 *  @Author B. Schwenker, University of Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class TBHit {
   
 private: // members 
 
  // Hit covariance (2x2 matrix)
  Matrix2d Cov;
  // Local hit coordinates (2x1 matrix)
  Vector2d Coord;
  // Pointer to raw TrackerHit (default nullptr)
  lcio::TrackerHit * RawHitPtr;
  // Unique ID (default -1)
  int UniqueID;
  // Sensor ID
  short SensorID;
  // Cluster Quality 
  short Quality;

    
 public: // functions 
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Constructors 
  TBHit();
  TBHit(int newSensorID, double u, double v, double cov_u, double cov_v, double cov_uv, int quality);
  TBHit(int newSensorID, const Vector2d &coord, const Matrix2d &cov);
  TBHit(lcio::TrackerHit* lciohit);

  // Build LCIO TrackerHit
  lcio::TrackerHitImpl * MakeLCIOHit() const;  
  static lcio::TrackerHitImpl * MakeLCIOHit(int newSensorID, double u, double v, double cov_u, double cov_v, double cov_uv, int quality);
  
  // Get/Set sensorID 
  void SetSensorID(int newSensorID) { SensorID = short(newSensorID); }
  int GetSensorID() const { return SensorID; }
  
  // Dimension of a pixel hit
  int GetDim() const { return 2;}
  
  // Get/Set quality of raw cluster
  void SetQuality(int flag) { Quality=short(flag); }
  int GetQuality() const { return Quality; }
  
  // Get/Set measured hit coord
  void SetCoord(const Vector2d& aCoord) { Coord = aCoord; }
  const Vector2d&  GetCoord() const { return Coord; }
  // Get/Set measurment covariance 
  void SetCov(const Matrix2d& aCov ) { Cov = aCov; }
  const Matrix2d&  GetCov() const { return Cov; }

  void ScaleCov(double alpha) { Cov*=alpha; }


  // Get/Set original raw cluster data   
  void SetRawHit(lcio::TrackerHit* lciohit ) { RawHitPtr = lciohit; }
  lcio::TrackerHit* GetRawHit() const  { return RawHitPtr; }

  
  PixelCluster GetCluster() const;

  StripCluster GetStripCluster() const; 
  
  // Get/Set unique ID for hit  
  void SetUniqueID(int ID ) { UniqueID = ID; }
  int GetUniqueID() const { return UniqueID; }
   
  // Get 3D point in sensor coordinates 
  Vector3d GetLocalSpacePoint() const;
  
};

} // Namespace
#include<Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(depfet::TBHit)

#endif // TBHIT_H

