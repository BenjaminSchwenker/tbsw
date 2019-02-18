// PixelCluster implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>



#include "PixelCluster.h"
#include <algorithm>
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
    FloatVec rawDigits = Digits->getChargeValues();
    int size = rawDigits.size()/3; 
    
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
  
  std::string PixelCluster::getShape(int pixeltype, int vCellPeriod, int uCellPeriod, int etaBin) const 
  {
    // Compute cluster label string
    stringstream streamLabel;  
    streamLabel << "E" << etaBin << "P" << m_vStart % vCellPeriod << "." <<  m_uStart % uCellPeriod  
                << "." << pixeltype; 
    
    for (auto digit : m_sortedDigits ) {
      streamLabel << "D" << digit.m_cellIDV - m_vStart  <<  "." << digit.m_cellIDU - m_uStart;  
    } 
    return streamLabel.str();   
  }
   
  std::string PixelCluster::getType(int pixeltype, int vCellPeriod, int uCellPeriod) const 
  { 
    // Compute cluster type string
    stringstream streamType;    
    streamType << "P" << int(m_vStart % vCellPeriod) << "." <<  int(m_uStart % uCellPeriod) << "." << pixeltype;   
    
    for (auto digit : m_sortedDigits ) { 
      streamType << "D" << digit.m_cellIDV - m_vStart  <<  "." << digit.m_cellIDU - m_uStart;  
    } 
    return streamType.str();  
  }  
   PixelCluster::operator std::string() const{
       return getType()+" " + getShape();
   }
  double PixelCluster::computeEta(double thetaU, double thetaV) const
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
  
  int PixelCluster::computeEtaBin(double eta, const vector<double>& etaBinEdges) const
  {
    for (int i = etaBinEdges.size() - 1; i >= 0; --i) {
      if (eta >= etaBinEdges[i])
        return i;
    }
    return 0;
  }
  
  int PixelCluster::getHeadPixelIndex(double thetaU, double thetaV) const
  {
    if (thetaV >= 0) {
      if (thetaU >= 0) {
        return getLastPixelWithVOffset(m_vSize - 1);
      } else {
        return getFirstPixelWithVOffset(m_vSize - 1);
      }
    } else {
      if (thetaU >= 0) {
        return getLastPixelWithVOffset(0);
      } else {
        return getFirstPixelWithVOffset(0);
      }
    }
  }
  
  int PixelCluster::getTailPixelIndex(double thetaU, double thetaV) const
  {
    if (thetaV >= 0) {
      if (thetaU >= 0) {
        return getFirstPixelWithVOffset(0);
      } else {
        return getLastPixelWithVOffset(0);
      }
    } else {
      if (thetaU >= 0) {
        return getFirstPixelWithVOffset(m_vSize - 1);
      } else {
        return getLastPixelWithVOffset(m_vSize - 1);
      }
    }
  }
  
  int PixelCluster::getLastPixelWithVOffset(int vOffset) const
  {
    for (auto index = 0; index < m_sortedDigits.size(); ++index) {
      int v = m_sortedDigits[index].m_cellIDV - m_vStart;
      if (vOffset < v) {
        if (index == 0) {
          return 0; 
        } else {
          return index-1;
        }
      }
    }
    return m_clsSize-1;
  }
  
  int PixelCluster::getFirstPixelWithVOffset(int vOffset) const
  {
    for (auto index = 0; index < m_sortedDigits.size(); ++index) {
      int v = m_sortedDigits[index].m_cellIDV - m_vStart;
      if (vOffset == v) {
        return index;
      }
    }
    return 0;
  }
  
  void PixelCluster::getCenterOfGravity(Det& Sensor, double& u, double& v, double& cov_u, double& cov_v, double& cov_uv) const { 
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
        
    // Compute the diagonal elements of the hit covariance matrix in mm2
    cov_u = pow(getUSize()*Sensor.GetPitchU(getVStart(), getUStart())/sqrt(12),2);
    cov_v = pow(getVSize()*Sensor.GetPitchV(getVStart(), getUStart())/sqrt(12),2); 
    cov_uv = 0.0;   
  }
}

