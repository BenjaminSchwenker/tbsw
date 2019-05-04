// PolyClusterDescriptor implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>



#include "PolyClusterDescriptor.h"
#include <algorithm>
#include <numeric>
#include <limits>
#include <iterator>
#include <cmath>
#define FMT_HEADER_ONLY 1
#include <fmt/format.h>
using namespace std;

namespace depfet {
  
    
  /** Constructor */
  PolyClusterDescriptor::PolyClusterDescriptor( const PixelCluster& Cluster, const Det& Sensor) : 
                                  m_originU(0), m_originV(0), m_uStart(0), m_vStart(0), m_lowerRight(0), m_upperLeft(0)
  {  
    m_uStart = Cluster.getUStart(); 
    m_vStart = Cluster.getVStart(); 

    m_originU = std::numeric_limits<float>::max();
    m_originV = std::numeric_limits<float>::max();
    
    // Create PolyRawDigit vector 
    const vector<RawDigit>& rawDigits = Cluster.getRawDigits();
    size_t size = rawDigits.size(); 
    m_sortedDigits.reserve(size);
    
    for ( auto&& digit : rawDigits ) {
      float posU = Sensor.GetPixelCenterCoordU(digit.m_cellIDV, digit.m_cellIDU);
      float posV = Sensor.GetPixelCenterCoordV(digit.m_cellIDV, digit.m_cellIDU);
      int pixelType = Sensor.GetPixelType(digit.m_cellIDV, digit.m_cellIDU);
      m_sortedDigits.push_back(PolyRawDigit(posU,posV,digit.m_charge,pixelType));
      
      if (posU < m_originU) m_originU = posU; 
      if (posV < m_originV) m_originV = posV; 
    }
    
    // Sort PolyRawDigit according to position of pixel centers
    std::sort(m_sortedDigits.begin(), m_sortedDigits.end());
    
    // Find index of lower right digit 
    m_lowerRight = size-1; 
    for (size_t index = 0; index < size-1; ++index) {
      if (m_sortedDigits[index+1].m_cellPosV > m_sortedDigits[0].m_cellPosV) {
        m_lowerRight = index;
        break;   
      }
    }    
    
    // Find index of upper left digit
    m_upperLeft = 0; 
    for (size_t index = size-1; index >0 ; --index) {
      if (m_sortedDigits[index-1].m_cellPosV < m_sortedDigits[size-1].m_cellPosV) {
        m_upperLeft = index;
        break;   
      }
    }    
  }
  
  std::string PolyClusterDescriptor::getShape(int vCellPeriod, int uCellPeriod, int etaBin) const 
  {
    std::string ret;
    ret.reserve(9+12*m_sortedDigits.size());
    ret.append(fmt::format("E{}P{}.{}.{}",etaBin,m_vStart%vCellPeriod,m_uStart%uCellPeriod,m_sortedDigits[0].m_pixelType));
    for (auto&& digit : m_sortedDigits ) {
      int tmpV = int(std::round( 10000*(digit.m_cellPosV - m_originV) ));
      int tmpU = int(std::round( 10000*(digit.m_cellPosU - m_originU) ));
      ret.append(fmt::format("D{}.{}.{}", tmpV, tmpU, digit.m_pixelType));
    }
    return ret;
  }
  
  std::string PolyClusterDescriptor::getType(int vCellPeriod, int uCellPeriod) const 
  { 
    std::string ret;
    ret.reserve(7+12*m_sortedDigits.size());
    ret.append(fmt::format("P{}.{}.{}",m_vStart%vCellPeriod,m_uStart%uCellPeriod,m_sortedDigits[0].m_pixelType));
    for (auto&& digit : m_sortedDigits ) {
      int tmpV = int(std::round( 10000*(digit.m_cellPosV - m_originV) ));
      int tmpU = int(std::round( 10000*(digit.m_cellPosU - m_originU) ));
      ret.append(fmt::format("D{}.{}.{}", tmpV, tmpU, digit.m_pixelType));
    }
    return ret;
  }  
  
  std::string PolyClusterDescriptor::getEtaBinString(int etaBin) const
  {
    std::string ret;
    ret.reserve(9+12*m_sortedDigits.size());
    ret.append(fmt::format("E{}",etaBin));
    return ret;
  }
   
  double PolyClusterDescriptor::computeEta(double thetaU, double thetaV) const
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
  
  int PolyClusterDescriptor::computeEtaBin(double eta, const vector<double>& etaBinEdges) 
  {
    for (int i = etaBinEdges.size() - 1; i >= 0; --i) {
      if (eta >= etaBinEdges[i])
        return i;
    }
    return 0;
  }
  
  int PolyClusterDescriptor::getHeadPixelIndex(double thetaU, double thetaV) const
  {
    // This logic estimates at which digit the particle leaves the sensor
    // Most meaningfull when thetaU and theta are not both zero.  
    if (thetaV >= 0) {
      if (thetaU >= 0) {
        return m_sortedDigits.size()-1; // This is the index of upper right digit 
      } else {
        return m_upperLeft; 
      }
    } else {
      if (thetaU >= 0) {
        return m_lowerRight; 
      } else {
        return 0; // This is the index of lower left digit
      }
    }
  }
  
  int PolyClusterDescriptor::getTailPixelIndex(double thetaU, double thetaV) const
  {
    // This logic estimates at which digit the particle enters the sensor
    // Most meaningfull when thetaU and theta are not both zero.  
    if (thetaV >= 0) {
      if (thetaU >= 0) {
        return 0; // This is the index of lower left digit
      } else {
        return m_lowerRight;
      }
    } else {
      if (thetaU >= 0) {
        return m_upperLeft;
      } else {
        return m_sortedDigits.size()-1; // This is the index of upper right digit 
      }
    }
  }
}

