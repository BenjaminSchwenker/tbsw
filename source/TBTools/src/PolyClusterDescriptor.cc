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
                                  m_originU(0), m_originV(0)
  {  
    const vector<RawDigit>& rawDigits = Cluster.getRawDigits();
    int size = rawDigits.size(); 
    m_sortedDigits.reserve(size);
    
    m_originU = std::numeric_limits<float>::max();
    m_originV = std::numeric_limits<float>::max();
    
    for ( auto&& digit : rawDigits ) {
      float posU = Sensor.GetPixelCenterCoordU(digit.m_cellIDV, digit.m_cellIDU);
      float posV = Sensor.GetPixelCenterCoordV(digit.m_cellIDV, digit.m_cellIDU);
      int pixelType = Sensor.GetPixelType(digit.m_cellIDV, digit.m_cellIDU);
      m_sortedDigits.push_back(PolyRawDigit(posU,posV,digit.m_charge,pixelType));
      
      if (posU < m_originU) m_originU = posU; 
      if (posV < m_originV) m_originV = posV; 
    }
    
    std::sort(m_sortedDigits.begin(), m_sortedDigits.end());
  }
  
  std::string PolyClusterDescriptor::getShape(int /*vCellPeriod*/, int /*uCellPeriod*/, int etaBin) const 
  {
    std::string ret;
    ret.reserve(9+12*m_sortedDigits.size());
    ret.append(fmt::format("E{}P{}.{}.{}",etaBin,0,0,m_sortedDigits[0].m_pixelType));
    for (auto&& digit : m_sortedDigits ) {
      ret.append(fmt::format("D{}.{}.{}", int(10000*(digit.m_cellPosV - m_originV)), int(10000*(digit.m_cellPosU - m_originU)), digit.m_pixelType));
    }
    return ret;
  }
  
  std::string PolyClusterDescriptor::getType(int /*vCellPeriod*/, int /*uCellPeriod*/) const 
  { 
    std::string ret;
    ret.reserve(7+12*m_sortedDigits.size());
    ret.append(fmt::format("P{}.{}.{}",0,0,m_sortedDigits[0].m_pixelType));
    for (auto&& digit : m_sortedDigits ) {
      ret.append(fmt::format("D{}.{}.{}", int(10000*(digit.m_cellPosV - m_originV)), int(10000*(digit.m_cellPosU - m_originU)), digit.m_pixelType));
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
    // This is a simplfied logic compared to PixelCluster. But it should work ok 
    // in most cases ... and can serve as a starting point. 
    if (thetaV >= 0) {
      if (thetaU >= 0) {
        return m_sortedDigits.size()-1;
      } else {
        return m_sortedDigits.size()-1;
      }
    } else {
      if (thetaU >= 0) {
        return 0;
      } else {
        return 0;
      }
    }
  }
  
  int PolyClusterDescriptor::getTailPixelIndex(double thetaU, double thetaV) const
  {
    // This is a simplfied logic compared to PixelCluster. But it should work ok 
    // in most cases ... and can serve as a starting point. 
    if (thetaV >= 0) {
      if (thetaU >= 0) {
        return 0;
      } else {
        return 0;
      }
    } else {
      if (thetaU >= 0) {
        return m_sortedDigits.size()-1;
      } else {
        return m_sortedDigits.size()-1;
      }
    }
  }
}

