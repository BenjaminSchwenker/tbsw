// PixelCluster implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>



#include "PixelCluster.h"
#include <algorithm>
#include <limits>
#include <iterator>

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
  
  
  std::string PixelCluster::getLabel(int scale) const 
  {
    // Compute cluster label string
    stringstream streamLabel;         
    streamLabel << m_clsSize; 
     
    for (auto digit : m_sortedDigits ) {
      streamLabel << "D" << digit.m_cellIDV - m_vStart  <<  "." << digit.m_cellIDU - m_uStart << "." << digit.m_charge/scale;  
    } 
    return streamLabel.str();   
  }

  
  std::string PixelCluster::getLabel(const std::vector<int>& jumps) const 
  {
    // Compute cluster label string
    stringstream streamLabel;         
    streamLabel << m_clsSize; 
     
    for (auto digit : m_sortedDigits ) {
      int mapped_charge = std::distance(jumps.cbegin(), std::upper_bound(jumps.begin(), jumps.end(), digit.m_charge) );
      streamLabel << "D" << digit.m_cellIDV - m_vStart  <<  "." << digit.m_cellIDU - m_uStart << "." << mapped_charge;  
    } 
    return streamLabel.str();   
  }
  

    
  std::string PixelCluster::getType() const 
  { 
    // Compute cluster type string
    stringstream streamType;    
    streamType << m_clsSize;     
    
    for (auto digit : m_sortedDigits ) { 
      streamType << "D" << digit.m_cellIDV - m_vStart  <<  "." << digit.m_cellIDU - m_uStart;  
    } 
    return streamType.str();  
  }

}

