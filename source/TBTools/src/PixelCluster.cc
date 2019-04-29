// PixelCluster implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>



#include "PixelCluster.h"
#include <algorithm>
#include <numeric>
#include <limits>
#include <iterator>
#include <cmath>

using namespace std;
using namespace lcio; 


namespace depfet {
  
    
  /** Constructor */
  PixelCluster::PixelCluster( TrackerData * Digits, unsigned short sensorID) : 
                                  m_sensorID(sensorID), m_clsCharge(0), m_seedCharge(0),
                                  m_clsSize(0), m_uSize(0), m_vSize(0), m_uStart(0), m_vStart(0)
  {
    const EVENT::FloatVec &rawDigits = Digits->getChargeValues();
    int size = rawDigits.size()/3; 
    m_sortedDigits.reserve(size);

    unsigned short vMin = std::numeric_limits<unsigned short>::max();
    unsigned short vMax = 0;
    unsigned short uMin = std::numeric_limits<unsigned short>::max();
    unsigned short uMax = 0; 

    for ( int index=0; index<size;  index++) { 
      unsigned short iU = static_cast<unsigned short> (rawDigits[index * 3]);
      unsigned short iV = static_cast<unsigned short> (rawDigits[index * 3 + 1]);
      unsigned short charge = static_cast<unsigned short> (rawDigits[index * 3 + 2]);   
      
      m_clsCharge += charge;
      
      if ( m_seedCharge < charge ) m_seedCharge = charge;  
      if (iV < vMin) vMin = iV; 
      if (iV > vMax) vMax = iV;
      if (iU < uMin) uMin = iU;
      if (iU > uMax) uMax = iU;  
      
      m_sortedDigits.push_back(RawDigit(iU,iV,charge));
    }
    
    // Compute some variables 
    m_clsSize = size;
    m_uStart = uMin; 
    m_vStart = vMin; 
    m_uSize = uMax-uMin+1;
    m_vSize = vMax-vMin+1; 
    
    std::sort(m_sortedDigits.begin(), m_sortedDigits.end());
  }
  
  void PixelCluster::getCenterOfGravity(const Det& Sensor, double& u, double& v, double& cov_u, double& cov_v, double& cov_uv) const { 
    // Calculate hit coord in local frame in mm
    u = 0; 
    v = 0;  
     
    for (const auto&  digit: m_sortedDigits) {
      u += Sensor.GetPixelCenterCoordU(digit.m_cellIDV, digit.m_cellIDU)*digit.m_charge;
      v += Sensor.GetPixelCenterCoordV(digit.m_cellIDV, digit.m_cellIDU)*digit.m_charge;     
    }
       
    if ( m_clsCharge > 0)  {
      u /= m_clsCharge;
      v /= m_clsCharge; 
    } 
        
    // Compute the 2x2 hit covariance matrix 
    cov_u = pow(getUSize()*Sensor.GetPitchU(getVStart(), getUStart())/sqrt(12),2);
    cov_v = pow(getVSize()*Sensor.GetPitchV(getVStart(), getUStart())/sqrt(12),2); 
    cov_uv = 0.0;   
  }
}

