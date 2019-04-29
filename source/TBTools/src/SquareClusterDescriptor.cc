// SquareClusterDescriptor implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>



#include "SquareClusterDescriptor.h"
#include <iterator>
#include <cmath>
#define FMT_HEADER_ONLY 1
#include <fmt/format.h>
using namespace std;

namespace depfet {
  
    
  /** Constructor */
  SquareClusterDescriptor::SquareClusterDescriptor( const PixelCluster& Cluster, const Det& Sensor) : 
                                  m_originU(0), m_originV(0)
  {  
    // FIXME: There should be a way to avoid copying the raw digits
    m_sortedDigits = Cluster.getRawDigits();
    m_vSize = Cluster.getVSize();  
    m_uStart = Cluster.getUStart();
    m_vStart = Cluster.getVStart();
    m_pixelType = Sensor.GetPixelType(m_vStart, m_uStart);
    m_originU = Sensor.GetPixelCenterCoordU(m_vStart, m_uStart);
    m_originV = Sensor.GetPixelCenterCoordV(m_vStart, m_uStart);   
  }
  
  std::string SquareClusterDescriptor::getShape(int vCellPeriod, int uCellPeriod, int etaBin) const 
  {
    // Compute shape string
    std::string ret;
    ret.reserve(9+4*m_sortedDigits.size());
    ret.append(fmt::format("E{}P{}.{}.{}",etaBin,m_vStart % vCellPeriod,m_uStart % uCellPeriod,m_pixelType));
    for (auto digit : m_sortedDigits ) {
        ret.append(fmt::format("D{}.{}",digit.m_cellIDV - m_vStart ,digit.m_cellIDU - m_uStart));
    }
    return ret;
  }
  
  std::string SquareClusterDescriptor::getType(int vCellPeriod, int uCellPeriod) const 
  { 
    std::string ret;
    ret.reserve(7+4*m_sortedDigits.size());
    ret.append(fmt::format("P{}.{}.{}",m_vStart % vCellPeriod,m_uStart % uCellPeriod,m_pixelType));
    for (auto digit : m_sortedDigits ) {
        ret.append(fmt::format("D{}.{}",digit.m_cellIDV - m_vStart ,digit.m_cellIDU - m_uStart));
    }
    return ret;
  }  

  std::string SquareClusterDescriptor::getEtaBinString(int etaBin) const
  {
    std::string ret;
    ret.reserve(9+4*m_sortedDigits.size());
    ret.append(fmt::format("E{}",etaBin));
    return ret;
  }
  
  double SquareClusterDescriptor::computeEta(double thetaU, double thetaV) const
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
  
  int SquareClusterDescriptor::computeEtaBin(double eta, const vector<double>& etaBinEdges) 
  {
    for (int i = etaBinEdges.size() - 1; i >= 0; --i) {
      if (eta >= etaBinEdges[i])
        return i;
    }
    return 0;
  }
  
  int SquareClusterDescriptor::getHeadPixelIndex(double thetaU, double thetaV) const
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
  
  int SquareClusterDescriptor::getTailPixelIndex(double thetaU, double thetaV) const
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
  
  int SquareClusterDescriptor::getLastPixelWithVOffset(int vOffset) const
  {
    for (size_t index = 0; index < m_sortedDigits.size(); ++index) {
      int v = m_sortedDigits[index].m_cellIDV - m_vStart;
      if (vOffset < v) {
        if (index == 0) {
          return 0; 
        } else {
          return index-1;
        }
      }
    }
    return m_sortedDigits.size()-1;
  }
  
  int SquareClusterDescriptor::getFirstPixelWithVOffset(int vOffset) const
  {
    for (size_t index = 0; index < m_sortedDigits.size(); ++index) {
      int v = m_sortedDigits[index].m_cellIDV - m_vStart;
      if (vOffset == v) {
        return index;
      }
    }
    return 0;
  }  
}

