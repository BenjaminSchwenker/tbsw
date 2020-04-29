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
                                  m_originU(0), m_originV(0), m_uStart(0), m_vStart(0), m_lowerRight(0), m_upperLeft(0), m_upperRight(0),  m_lowerLeft(0)
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
    
    // The digits are sorted from lowerLeft to 
    // upperRight. Easy to find their index. 
    m_lowerLeft = 0;
    m_upperRight = size-1;

    // Find index of lower right digit 
    m_lowerRight = m_lowerLeft; 
    for (size_t index = 0; index < size-1; ++index) {
      if (m_sortedDigits[index+1].m_cellPosV > m_sortedDigits[0].m_cellPosV) {
        m_lowerRight = index;
        break;   
      }
    }    
    
    // Find index of upper left digit
    m_upperLeft = m_upperRight; 
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
   
  double PolyClusterDescriptor::computeEta(double thetaU, double thetaV, int etaIndex) const
  {
    auto heads = getHeadPixels(thetaU, thetaV);
    auto tails = getTailPixels(thetaU, thetaV);
    double charge_head{1}, charge_tail{1}; 
    
    switch (etaIndex) {
      case 0: charge_head=m_sortedDigits[heads[0]].m_charge;
              charge_tail=m_sortedDigits[tails[0]].m_charge;
              break;
      case 1: charge_head=m_sortedDigits[heads[1]].m_charge;
              charge_tail=m_sortedDigits[tails[0]].m_charge;
              break;
      case 2: charge_head=m_sortedDigits[heads[0]].m_charge;
              charge_tail=m_sortedDigits[tails[1]].m_charge;
              break;
      case 3: charge_head=m_sortedDigits[heads[1]].m_charge;
              charge_tail=m_sortedDigits[tails[1]].m_charge;
              break;
      case 4: charge_head=m_sortedDigits[heads[0]].m_charge+m_sortedDigits[heads[1]].m_charge;
              charge_tail=m_sortedDigits[tails[0]].m_charge;
              break;
      case 5: charge_head=m_sortedDigits[heads[0]].m_charge+m_sortedDigits[heads[1]].m_charge;
              charge_tail=m_sortedDigits[tails[1]].m_charge;
              break;
      case 6: charge_head=m_sortedDigits[heads[0]].m_charge;
              charge_tail=m_sortedDigits[tails[0]].m_charge + m_sortedDigits[tails[1]].m_charge;
              break;
      case 7: charge_head=m_sortedDigits[heads[1]].m_charge;
              charge_tail=m_sortedDigits[tails[0]].m_charge + m_sortedDigits[tails[1]].m_charge;
              break;
      case 8: charge_head=m_sortedDigits[heads[0]].m_charge + m_sortedDigits[heads[1]].m_charge;
              charge_tail=m_sortedDigits[tails[0]].m_charge + m_sortedDigits[tails[1]].m_charge;
              break;
    }          
    
    double eta = 0.5;
    if (charge_head != charge_tail) eta=charge_tail / (charge_tail+charge_head);
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
  
  std::vector<int> PolyClusterDescriptor::getHeadPixels(double thetaU, double thetaV) const
  {
    
    // Basically, we need to distingush same sign 
    // and opposite sign cases for angles.  
    std::vector<int> heads{0,0};
   
    if (thetaV*thetaU >= 0) { 
      heads[0] = m_upperRight;
      heads[1] = m_upperRight;
      
      if (m_sortedDigits.size() > 1) {
        heads[1] = m_upperRight-1;     
      }    
    } else {
      heads[0] = m_upperLeft;
      heads[1] = m_upperLeft;
      
      if (m_sortedDigits.size() > 1) {
        if (m_upperLeft < m_upperRight) {
          heads[1] = m_upperLeft+1;
        } else { 
          heads[1] = m_upperLeft-1;
        }  
      }
    }
    return heads; 
  }
  
  std::vector<int> PolyClusterDescriptor::getTailPixels(double thetaU, double thetaV) const
  {
    // Basically, we need to distingush same sign 
    // and opposite sign cases for angles.  
    std::vector<int> tails{0,0};
     
    if (thetaV*thetaU >=0) { 
      tails[0] = m_lowerLeft; 
      tails[1] = m_lowerLeft; 
      
      if (m_sortedDigits.size() > 1) {
        tails[1] = m_lowerLeft+1; 
      }    
    } else {
      tails[0] = m_lowerRight; 
      tails[1] = m_lowerRight; 
      
      if (m_sortedDigits.size() > 1) {
        if (m_lowerRight > m_lowerLeft) {
          tails[1] = m_lowerRight-1; 
        } else { 
          tails[1] = m_lowerRight+1; 
        }   
      }
    }
    return tails; 
  }
}

