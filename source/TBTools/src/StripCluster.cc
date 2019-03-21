// StripCluster implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>



#include "StripCluster.h"

using namespace std;
using namespace lcio; 

namespace depfet {

/** Constructor */
StripCluster::StripCluster( LCObjectVec& DigitVec, unsigned short sensorID, unsigned short clsType) : 
   m_sensorID(sensorID), m_uClsCharge(0), m_uSeedCharge(0), m_vClsCharge(0), m_vSeedCharge(0),
   m_uSize(0), m_vSize(0), m_uStart(0), m_vStart(0), m_clsType(clsType)
{
  
  if ( DigitVec.size() > 0 ) {
    
    // This is the u side cluster
    TrackerData * clusterDigits = dynamic_cast<TrackerData *> ( DigitVec[0] );
    
    FloatVec digits = clusterDigits->getChargeValues();
    m_uSize = digits.size()/3; 
    
    int min = 9999999;
    int max = 0; 
     
    for ( int index=0; index<m_uSize;  index++) { 
         
      int cell = static_cast<int> (digits[index * 3 + 1]);
      float charge =  digits[index * 3 + 2];
      
      m_uClsCharge += charge;
       
      if ( m_uSeedCharge < charge ) m_uSeedCharge = charge; 
       
      if (cell < min) min = cell; 
      if (cell > max) max = cell;
     
          
    }
  
    m_uStart = min;  
    m_uSize = max-min+1;

  } 

  if ( DigitVec.size() > 1 ) {
    
    // This is the v side cluster
    TrackerData * clusterDigits = dynamic_cast<TrackerData *> ( DigitVec[1] );
    
    FloatVec digits = clusterDigits->getChargeValues();
    m_vSize = digits.size()/3; 
    
    int min = 9999999;
    int max = 0; 
     
    for ( int index=0; index<m_vSize;  index++) { 
         
      int cell = static_cast<int> (digits[index * 3 + 1]);
      float charge =  digits[index * 3 + 2];
      
      m_vClsCharge += charge;
       
      if ( m_vSeedCharge < charge ) m_vSeedCharge = charge; 
       
      if (cell < min) min = cell; 
      if (cell > max) max = cell;
               
    }
  
    m_vStart = min; 
    m_vSize = max-min+1;

  } 
}



}

