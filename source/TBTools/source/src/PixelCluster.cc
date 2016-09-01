// PixelCluster implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>



#include "PixelCluster.h"

using namespace std;
using namespace lcio; 

namespace depfet {

/** Constructor */
PixelCluster::PixelCluster( TrackerData * Digits, unsigned short sensorID, unsigned short clsType) : 
   m_sensorID(sensorID), m_clsCharge(0), m_seedCharge(0),
   m_clsSize(0), m_uSize(0), m_vSize(0), m_uStart(0), m_vStart(0), m_clsType(clsType)
{
  
  FloatVec sparsePixels = Digits->getChargeValues();
  m_clsSize = sparsePixels.size()/3; 

  int vMin = 9999999;
  int vMax = 0;
  int uMin = 9999999;
  int uMax = 0; 
  
  for ( int index=0; index<m_clsSize;  index++) { 
         
    int uCell = static_cast<int> (sparsePixels[index * 3]);
    int vCell = static_cast<int> (sparsePixels[index * 3 + 1]);
    float charge =  sparsePixels[index * 3 + 2];
    
    m_clsCharge += charge;
    
    if ( m_seedCharge < charge ) m_seedCharge = charge; 
    
    if (vCell < vMin) vMin = vCell; 
    if (vCell > vMax) vMax = vCell;
    if (uCell < uMin) uMin = uCell;
    if (uCell > uMax) uMax = uCell;
          
  }
  
  m_uStart = uMin; 
  m_vStart = vMin; 
  m_uSize = uMax-uMin+1;
  m_vSize = vMax-vMin+1; 
}



}

