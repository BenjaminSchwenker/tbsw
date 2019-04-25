// PolyPixelCluster implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>



#include "PolyPixelCluster.h"
#include <algorithm>
#include <numeric>
#include <limits>
#include <iterator>
#include <cmath>
#define FMT_HEADER_ONLY 1
#include <fmt/format.h>
using namespace std;
using namespace lcio; 

namespace depfet {
  
    
  /** Constructor */
  PolyPixelCluster::PolyPixelCluster( TrackerData * Digits, const Det& Sensor) : 
                                  m_sensorID(0), m_clsCharge(0), m_seedCharge(0),
                                  m_clsSize(0), m_uSize(0), m_vSize(0), m_uStart(0), m_vStart(0)
  {
    m_sensorID = Sensor.GetSensorID();
    
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
      
      float posU = Sensor.GetPixelCenterCoordU(iV, iU);
      float posV = Sensor.GetPixelCenterCoordV(iV, iU);
      int pixelType = Sensor.GetPixelType(iV, iU);
      m_sortedDigits.push_back(PolyRawDigit(posU,posV,charge,pixelType));
    }
    
    // Compute some variables 
    m_clsSize = size;
    m_uSize = uMax-uMin+1;
    m_vSize = vMax-vMin+1; 
    m_uStart = uMin; 
    m_vStart = vMin;  
    std::sort(m_sortedDigits.begin(), m_sortedDigits.end());
  }
  
  std::string PolyPixelCluster::getShape(int pixeltype, int vCellPeriod, int uCellPeriod, int etaBin) const 
  {
    float uStartPos = m_sortedDigits[0].m_cellPosU;   
    float vStartPos = m_sortedDigits[0].m_cellPosV;     

    std::string ret;
    ret.reserve(9+12*m_sortedDigits.size());
    ret.append(fmt::format("E{}P{}.{}.{}",etaBin,m_vStart % vCellPeriod,m_uStart % uCellPeriod,pixeltype));
    for (auto digit : m_sortedDigits ) {
        ret.append(fmt::format("D{:.0f}.{:.0f}.{}", 10000*(digit.m_cellPosV - vStartPos), 10000*(digit.m_cellPosU - uStartPos), digit.m_pixelType));
    }
    return ret;
  }
  
  std::string PolyPixelCluster::getType(int pixeltype, int vCellPeriod, int uCellPeriod) const 
  { 
    float uStartPos = m_sortedDigits[0].m_cellPosU;   
    float vStartPos = m_sortedDigits[0].m_cellPosV;     
    
    std::string ret;
    ret.reserve(7+12*m_sortedDigits.size());
    ret.append(fmt::format("P{}.{}.{}",m_vStart % vCellPeriod,m_uStart % uCellPeriod,pixeltype));
    for (auto digit : m_sortedDigits ) {
        ret.append(fmt::format("D{:.0f}.{:.0f}.{}", 10000*(digit.m_cellPosV - vStartPos), 10000*(digit.m_cellPosU - uStartPos), digit.m_pixelType));
    }
    return ret;
  }  
  
  std::string PolyPixelCluster::getEtaBinString(int etaBin) const
  {
    std::string ret;
    ret.reserve(9+12*m_sortedDigits.size());
    ret.append(fmt::format("E{}",etaBin));
    return ret;
  }
  
  PolyPixelCluster::operator std::string() const{
    return getType()+" " + getShape();
  }
  
  double PolyPixelCluster::computeEta(double thetaU, double thetaV) const
  {
    auto headPixelIndex = getHeadPixelIndex(thetaU, thetaV);
    auto tailPixelIndex = getTailPixelIndex(thetaU, thetaV);
    float eta = 0;
    if (headPixelIndex != tailPixelIndex) {
      eta = m_sortedDigits[tailPixelIndex].m_charge;
      eta /= (m_sortedDigits[tailPixelIndex].m_charge + m_sortedDigits[headPixelIndex].m_charge);
    } else {
      eta = 0.5;
    }
    return eta;
  }
  
  int PolyPixelCluster::computeEtaBin(double eta, const vector<double>& etaBinEdges) const
  {
    for (int i = etaBinEdges.size() - 1; i >= 0; --i) {
      if (eta >= etaBinEdges[i])
        return i;
    }
    return 0;
  }
  
  int PolyPixelCluster::getHeadPixelIndex(double thetaU, double thetaV) const
  {
    // This is a simplfied logic compared to PixelCluster. But it should work ok 
    // in most cases ... and can serve as a starting point. 
    if (thetaV >= 0) {
      if (thetaU >= 0) {
        //return getLastPixelWithVOffset(m_vSize - 1);
        return m_clsSize-1;
      } else {
        //return getFirstPixelWithVOffset(m_vSize - 1);
        return m_clsSize-1;
      }
    } else {
      if (thetaU >= 0) {
        //return getLastPixelWithVOffset(0);
        return 0;
      } else {
        //return getFirstPixelWithVOffset(0);
        return 0;
      }
    }
  }
  
  int PolyPixelCluster::getTailPixelIndex(double thetaU, double thetaV) const
  {
    // This is a simplfied logic compared to PixelCluster. But it should work ok 
    // in most cases ... and can serve as a starting point. 
    if (thetaV >= 0) {
      if (thetaU >= 0) {
        //return getFirstPixelWithVOffset(0);
        return 0;
      } else {
        //return getLastPixelWithVOffset(0);
        return 0;
      }
    } else {
      if (thetaU >= 0) {
        //return getFirstPixelWithVOffset(m_vSize - 1);
        return m_clsSize-1;
      } else {
        //return getLastPixelWithVOffset(m_vSize - 1);
        return m_clsSize-1;
      }
    }
  }
    
  void PolyPixelCluster::getCenterOfGravity(const Det& Sensor, double& u, double& v, double& cov_u, double& cov_v, double& cov_uv) const { 
    // Calculate hit coord in local frame in mm
    u = 0; 
    v = 0;  
     
    for (const auto&  digit: m_sortedDigits) {
      u += digit.m_cellPosU*digit.m_charge;
      v += digit.m_cellPosV*digit.m_charge;     
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

