// Det implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// TBTools includes 
#include "Det.h"
#include "MaterialEffect.h"

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

#include <algorithm>

// Namespaces
using namespace marlin;

namespace depfet {



/** Dummy constructor for pixel modules
 *
 * Creates a dummy pixel module. 
 */
Det::Det() {}


void Det::SetCellsU( std::vector< std::tuple<int,int,double> > uCells)
{ 
  _uCells = uCells;
   
  // First of all, sort the cell vector 
  std::sort(std::begin(_uCells), std::end(_uCells), [](auto const &t1, auto const &t2) {
    return std::get<0>(t1) < std::get<0>(t2); 
  });
  
  // Secondly, avoid an empty cell vector 
  if ( _uCells.size() == 0 ) { 
    streamlog_out(ERROR) << "SensorID=" << _ID << "has no ucells! Add a default ucell"
                         << std::endl;     
    _uCells.push_back( std::tuple<int, int, double>(0, 0, 0) );
  }

  _minCellU = std::get<0>(_uCells.at(0)); 
  _nCellsU = 0;  
  _sensitiveSizeU = 0;   
  
  for (auto group : _uCells ) {
    int minCell = std::get<0>(group);
    int maxCell = std::get<1>(group);
    double pitch = std::get<2>(group);   
     
    if ( maxCell < minCell ) { 
      streamlog_out(ERROR) << "SensorID=" << _ID << "has cell group with maxCell < minCell!"
                           << std::endl;     
    }
     
    if (pitch <= 0) {
     streamlog_out(ERROR) << "SensorID=" << _ID << "has cell group with pitch <= 0!"
                          << std::endl;
    }
    
    // Add cells from group
    _nCellsU += (maxCell - minCell + 1); 
    
    // Add offset for group
    _offsetsU.push_back(_sensitiveSizeU);
    
    // Compute offset for next group
    _sensitiveSizeU += (maxCell - minCell + 1)*pitch; 
  }
  
}

void Det::SetCellsV( std::vector< std::tuple<int,int,double> > vCells)
{ 
  _vCells = vCells;
   
  // First of all, sort the cell vector 
  std::sort(std::begin(_vCells), std::end(_vCells), [](auto const &t1, auto const &t2) {
    return std::get<0>(t1) < std::get<0>(t2); 
  });
  
  // Secondly, avoid an empty cell vector 
  if ( _vCells.size() == 0 ) { 
    streamlog_out(ERROR) << "SensorID=" << _ID << "has no vcells! Add a default ucell"
                         << std::endl;     
    _vCells.push_back( std::tuple<int, int, double>(0, 0, 0) );
  }
  
  _minCellV = std::get<0>(_vCells.at(0)); 
  _nCellsV = 0;  
  _sensitiveSizeV = 0;   
  
  for (auto group : _vCells ) {
    int minCell = std::get<0>(group);
    int maxCell = std::get<1>(group);
    double pitch = std::get<2>(group);   

    if ( maxCell < minCell ) { 
      streamlog_out(ERROR) << "SensorID=" << _ID << "has cell group with maxCell < minCell!"
                           << std::endl;     
    }
     
    if (pitch <= 0) {
     streamlog_out(ERROR) << "SensorID=" << _ID << "has cell group with pitch <= 0!"
                          << std::endl;
    }

    // Add cells from group
    _nCellsV += (maxCell - minCell + 1); 
    
    // Add offset for group
    _offsetsV.push_back(_sensitiveSizeV);
  
    // Compute offset for next group
    _sensitiveSizeV += (maxCell - minCell + 1)*pitch; 
  }

}



int Det::GetNCellsU() const
{
  return _nCellsU; 
}  

int Det::GetNCellsV() const
{
  return _nCellsV; 
}   

double Det::GetSensitiveSizeU() const
{
  return _sensitiveSizeU; 
}  
  
double Det::GetSensitiveSizeV() const
{
  return _sensitiveSizeV;  
} 

int Det::GetPixelTypeU(int vcell, int ucell)  
{
  int i = 0; 
  for (auto group : _uCells ) {
    int minCell = std::get<0>(group);
    int maxCell = std::get<1>(group);
    if (ucell >= minCell && ucell <= maxCell) 
      break;     
    i++;
  }
  return i; 
}

int Det::GetPixelTypeV(int vcell, int ucell)  
{
  int i = 0; 
  for (auto group : _vCells ) {
    int minCell = std::get<0>(group);
    int maxCell = std::get<1>(group);
    if (vcell >= minCell && vcell <= maxCell) 
      break;     
    i++;
  }
  return i; 
} 

int Det::GetPixelType(int vcell, int ucell)   
{ 
  int iu = GetPixelTypeU(vcell, ucell); 
  int iv = GetPixelTypeV(vcell, ucell); 
  int nGroupsU = _uCells.size();
  return (nGroupsU*iv + iu);
}

double Det::GetPitchU(int vcell, int ucell)  
{
  auto group = _uCells.at(GetPixelTypeU(vcell, ucell));   
  return std::get<2>(group); 
} 
  
double Det::GetPitchV(int vcell, int ucell)
{
  auto group = _vCells.at(GetPixelTypeV(vcell, ucell));   
  return std::get<2>(group); 
}  

int Det::encodePixelID(int vcell, int ucell)
{
  return (GetNCellsU()*vcell + ucell);
}


void Det::decodePixelID(int vcell, int ucell, int uniqPixelID)
{
  vcell = uniqPixelID / GetNCellsU();
  ucell = uniqPixelID - vcell*GetNCellsU();
}
 
 	
/**  Check if sensitive volume is crossed
 */
bool Det::SensitiveCrossed(double u, double v, double w) const
{
  if (u < -(GetSensitiveSizeU())/2.  || u > (GetSensitiveSizeU())/2.) {
   return false;
  }
  if (v < -(GetSensitiveSizeV())/2. || v > (GetSensitiveSizeV())/2.) {
    return false;
  }
  if (w < -GetSensitiveThickness()/2. || w > GetSensitiveThickness()/2.) {
    return false;
  }
  return true; 
}


bool Det::isPointOutOfSensor( double u, double v, double w) 
{
  bool isOut = false; 
  
  // Boundary set +- epsilon
  if ( (u < (-GetSensitiveSizeU()/2. -0.005)) || (u > (+GetSensitiveSizeU()/2. + 0.005)) ||
       (v < (-GetSensitiveSizeV()/2. -0.005)) || (v > (+GetSensitiveSizeV()/2. + 0.005)) ||
       (w < (-GetSensitiveThickness()/2. -0.005)) || (w > (+GetSensitiveThickness()/2. + 0.005)) ) isOut = true;

   // Return if out or not
   return isOut;
}
 	
/**  Check if module is crossed (including supports)
 */
bool Det::ModuleCrossed(double u, double v, double w) const
{  
  if (u < -GetLadderSizeU()/2.  || u > GetLadderSizeU()/2.) {
   return false;
  }
  if (v < -GetLadderSizeV()/2. || v > GetLadderSizeV()/2.) {
    return false;
  }
  if (w < -GetLadderThickness()/2. || w > GetLadderThickness()/2.) {
    return false;
  }
  return true; 
}

double Det::GetThickness(double u, double v) const
{
  if ( SensitiveCrossed(u, v) ) {
    return GetSensitiveThickness(); 
  } 
  if ( ModuleCrossed(u, v) ) {
    return GetLadderThickness(); 
  }
  return 0; 
}  

/** Get length of particle track (mm)
 */
double Det::GetTrackLength(double u, double v, double dudw, double dvdw) const
{
  return GetThickness(u,v)*std::sqrt(1 + dudw*dudw + dvdw*dvdw);  
}
    
double Det::GetRadLength(double u, double v) const
{ 
  if ( SensitiveCrossed(u, v) ) {
    return GetSensitiveRadLength(); 
  } 
  if ( ModuleCrossed(u, v) ) {
    return GetLadderRadLength(); 
  }
  return materialeffect::X0_air; 
} 

double Det::GetAtomicNumber(double u, double v) const
{ 
  if ( SensitiveCrossed(u, v) ) {
    return GetSensitiveAtomicNumber(); 
  } 
  if ( ModuleCrossed(u, v) ) {
    return GetLadderAtomicNumber(); 
  }
  return 8; 
} 

double Det::GetAtomicMass(double u, double v) const
{ 
  if ( SensitiveCrossed(u, v) ) {
    return GetSensitiveAtomicMass(); 
  } 
  if ( ModuleCrossed(u, v) ) {
    return GetLadderAtomicMass(); 
  }
  return 16; 
} 

/**  Get v coord of pixel center 
 */
double Det::GetPixelCenterCoordV(int vcell, int ucell)
{    
  int i = GetPixelTypeV(vcell, ucell);   
  double offset = _offsetsV.at(i);
  auto group = _vCells.at(i);
  int minCell = std::get<0>(group);  
  double pitch = std::get<2>(group);  
  
  // V coord measured from lower feft corner
  double v_coord = offset + pitch*(vcell - minCell) + pitch*0.5;
  // Ok, shift coord to sensor center
  v_coord -= 0.5*GetSensitiveSizeV();           
  
  return v_coord;
}
 
/**  Get u coord of pixel center 
 */
double Det::GetPixelCenterCoordU(int vcell, int ucell)
{
  int i = GetPixelTypeU(vcell, ucell);   
  double offset = _offsetsU.at(i);
  auto group = _uCells.at(i);
  int minCell = std::get<0>(group);  
  double pitch = std::get<2>(group);  
  
  // U coord measured from lower feft corner
  double u_coord = offset + pitch*(ucell - minCell) + pitch*0.5;
  // Ok, shift coord to sensor center
  u_coord -= 0.5*GetSensitiveSizeU();           
  
  return u_coord;
}

/**  Calculate pixel column from coord (u,v)
 */      
int Det::GetUCellFromCoord( double u, double v )
{
  if (u < -GetSensitiveSizeU()/2.) {
   return _minCellU;
  } else if (u > GetSensitiveSizeU()/2.) {
   return _nCellsU-1;
  } 
  
  int ucell = _minCellU;
  double offset = -0.5*GetSensitiveSizeU();
  for (auto group : _uCells ) {
    int minCell = std::get<0>(group);
    int maxCell = std::get<1>(group);
    double pitch = std::get<2>(group);   
    
    if ( u >= offset && u <= (offset + (maxCell - minCell + 1)*pitch) ) {
     ucell = floor( (u-offset) / pitch ) + minCell; 
     break;   
    }

    // Compute offset for next group
    offset += (maxCell - minCell + 1)*pitch; 
  }
  
  return ucell;  
} 
   
/**  Calculate pixel column from coord (u,v)
 */   
int Det::GetVCellFromCoord( double u, double v ) 
{
  if (v < -GetSensitiveSizeV()/2.) {
   return _minCellV;
  } else if (v > GetSensitiveSizeV()/2.) {
   return _nCellsV-1;
  } 
  
  int vcell = _minCellV;
  double offset = -0.5*GetSensitiveSizeV();
  for (auto group : _vCells ) {
    int minCell = std::get<0>(group);
    int maxCell = std::get<1>(group);
    double pitch = std::get<2>(group);   
    
    if ( v >= offset && v <= (offset + (maxCell - minCell + 1)*pitch) ) {
     vcell = floor( (v-offset) / pitch ) + minCell; 
     break;   
    }

    // Compute offset for next group
    offset += (maxCell - minCell + 1)*pitch; 
  }
  
  return vcell;  
}


/**  Print methods
 */ 
void Det::Print()
{
  streamlog_out(MESSAGE3) << std::endl;
  streamlog_out(MESSAGE3) << "  Plane Number:    " << GetPlaneNumber()    << std::endl;
  streamlog_out(MESSAGE3) << "  DAQ ID:          " << GetDAQID()          << std::endl;  
  streamlog_out(MESSAGE3) << "  NCellsU:         " << GetNCellsU()        << std::endl;
  streamlog_out(MESSAGE3) << "  NCellsV:         " << GetNCellsV()        << std::endl;  
}
 


} // Namespace;

