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
using namespace lcio; 

namespace depfet {
  
    
  /** Constructor */
  PolyClusterDescriptor::PolyClusterDescriptor( TrackerData * Digits, const Det& Sensor) : 
                                  m_originU(0), m_originV(0)
  {
    
    const EVENT::FloatVec &rawDigits = Digits->getChargeValues();
    int size = rawDigits.size()/3; 
    m_sortedDigits.reserve(size);

    m_originU = std::numeric_limits<float>::max();
    m_originV = std::numeric_limits<float>::max();
    
    for ( int index=0; index<size;  index++) { 
      unsigned short iU = static_cast<unsigned short> (rawDigits[index * 3]);
      unsigned short iV = static_cast<unsigned short> (rawDigits[index * 3 + 1]);
      unsigned short charge = static_cast<unsigned short> (rawDigits[index * 3 + 2]);   
      
      float posU = Sensor.GetPixelCenterCoordU(iV, iU);
      float posV = Sensor.GetPixelCenterCoordV(iV, iU);
      int pixelType = Sensor.GetPixelType(iV, iU);
      m_sortedDigits.push_back(PolyRawDigit(posU,posV,charge,pixelType));

      if (posU < m_originU) m_originU = posU; 
      if (posV < m_originV) m_originV = posV; 
    }
    
    std::sort(m_sortedDigits.begin(), m_sortedDigits.end());
  }
  
  std::string PolyClusterDescriptor::getShape(int pixeltype, int /*vCellPeriod*/, int /*uCellPeriod*/, int etaBin) const 
  {
    std::string ret;
    ret.reserve(9+12*m_sortedDigits.size());
    ret.append(fmt::format("E{}P{}.{}.{}",etaBin,0,0,pixeltype));
    for (auto digit : m_sortedDigits ) {
        ret.append(fmt::format("D{:.0f}.{:.0f}.{}", 10000*(digit.m_cellPosV - m_originV), 10000*(digit.m_cellPosU - m_originU), digit.m_pixelType));
    }
    return ret;
  }
  
  std::string PolyClusterDescriptor::getType(int pixeltype, int /*vCellPeriod*/, int /*uCellPeriod*/) const 
  { 
    std::string ret;
    ret.reserve(7+12*m_sortedDigits.size());
    ret.append(fmt::format("P{}.{}.{}",0,0,pixeltype));
    for (auto digit : m_sortedDigits ) {
        ret.append(fmt::format("D{:.0f}.{:.0f}.{}", 10000*(digit.m_cellPosV - m_originV), 10000*(digit.m_cellPosU - m_originU), digit.m_pixelType));
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
  
  int PolyClusterDescriptor::computeEtaBin(double eta, const vector<double>& etaBinEdges) const
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

