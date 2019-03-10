// SquareDet implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// TBTools includes 
#include "SquareDet.h"
#include "MaterialEffect.h"

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

#include <algorithm>

// Namespaces
using namespace marlin;

namespace depfet {



SquareDet::SquareDet(const std::string& typeName, int sensorID, int planeNumber, 
                     double sensThick, double sensRadLenght, double sensAtomicNumber,
                     double sensAtomicMass, double ladderThick, double ladderRadLength, 
                     double ladderAtomicNumber, double ladderAtomicMass, double ladderSizeU, 
                     double ladderSizeV, const std::vector< std::tuple<int,int,double> >& uCells, 
                     const std::vector< std::tuple<int,int,double> >& vCells, 
                     const ReferenceFrame& discrete, const ReferenceFrame& nominal )
  : Det(typeName, sensorID, planeNumber) 
{
  
  // Set u cells 
  SetCellsU( uCells);
      
  // Set v cells 
  SetCellsV( vCells );

  m_sensitiveThickness = sensThick;
  m_sensitiveRadLength = sensRadLenght;
  m_sensitiveAtomicNumber = sensAtomicNumber;
  m_sensitiveAtomicMass = sensAtomicMass;
  m_ladderThickness = ladderThick;
  m_ladderRadLength = ladderRadLength;
  m_ladderAtomicNumber = ladderAtomicNumber;  
  m_ladderAtomicMass = ladderAtomicMass;  
  m_ladderSizeU = ladderSizeU; 
  m_ladderSizeV = ladderSizeV; 
  m_discrete = discrete;
  m_nominal = nominal;
}

SquareDet::SquareDet(const std::string& typeName, int sensorID, int planeNumber)
  : Det(typeName, sensorID, planeNumber) {}



// TODO this code should be put into the class where the SquareDet get constructed from an XML file
void SquareDet::SetCellsU( std::vector< std::tuple<int,int,double> > uCells)
{ 
  m_uCells = uCells;
   
  // First of all, sort the cell vector 
  std::sort(std::begin(m_uCells), std::end(m_uCells), [](auto const &t1, auto const &t2) {
    return std::get<0>(t1) < std::get<0>(t2); 
  });
  
  // Secondly, avoid an empty cell vector 
  if ( m_uCells.size() == 0 ) { 
    streamlog_out(ERROR) << "SensorID=" << GetSensorID() << "has no ucells! Add a default ucell"
                         << std::endl;     
    m_uCells.push_back( std::tuple<int, int, double>(0, 0, 0) );
  }

  m_minCellU = std::get<0>(m_uCells.at(0)); 
  m_nCellsU = 0;  
  m_sensitiveSizeU = 0;   
  
  for (auto group : m_uCells ) {
    int minCell = std::get<0>(group);
    int maxCell = std::get<1>(group);
    double pitch = std::get<2>(group);   
     
    if ( maxCell < minCell ) { 
      streamlog_out(ERROR) << "SensorID=" << GetSensorID() << "has cell group with maxCell < minCell!"
                           << std::endl;     
    }
     
    if (pitch <= 0) {
     streamlog_out(ERROR) << "SensorID=" << GetSensorID() << "has cell group with pitch <= 0!"
                          << std::endl;
    }
    
    // Add cells from group
    m_nCellsU += (maxCell - minCell + 1); 
    
    // Add offset for group
    m_offsetsU.push_back(m_sensitiveSizeU);
    
    // Compute offset for next group
    m_sensitiveSizeU += (maxCell - minCell + 1)*pitch; 
  }
  
}

// TODO this code should be put into the class where the SquareDet get constructed from an XML file
void SquareDet::SetCellsV( std::vector< std::tuple<int,int,double> > vCells)
{ 
  m_vCells = vCells;
   
  // First of all, sort the cell vector 
  std::sort(std::begin(m_vCells), std::end(m_vCells), [](auto const &t1, auto const &t2) {
    return std::get<0>(t1) < std::get<0>(t2); 
  });
  
  // Secondly, avoid an empty cell vector 
  if ( m_vCells.size() == 0 ) { 
    streamlog_out(ERROR) << "SensorID=" << GetSensorID() << "has no vcells! Add a default ucell"
                         << std::endl;     
    m_vCells.push_back( std::tuple<int, int, double>(0, 0, 0) );
  }
  
  m_minCellV = std::get<0>(m_vCells.at(0)); 
  m_nCellsV = 0;  
  m_sensitiveSizeV = 0;   
  
  for (auto group : m_vCells ) {
    int minCell = std::get<0>(group);
    int maxCell = std::get<1>(group);
    double pitch = std::get<2>(group);   

    if ( maxCell < minCell ) { 
      streamlog_out(ERROR) << "SensorID=" << GetSensorID() << "has cell group with maxCell < minCell!"
                           << std::endl;     
    }
     
    if (pitch <= 0) {
     streamlog_out(ERROR) << "SensorID=" << GetSensorID() << "has cell group with pitch <= 0!"
                          << std::endl;
    }

    // Add cells from group
    m_nCellsV += (maxCell - minCell + 1); 
    
    // Add offset for group
    m_offsetsV.push_back(m_sensitiveSizeV);
  
    // Compute offset for next group
    m_sensitiveSizeV += (maxCell - minCell + 1)*pitch; 
  }

}

int SquareDet::GetMaxUCell() const
{
  return m_nCellsU-1;
}  

int SquareDet::GetMaxVCell() const
{
  return m_nCellsV-1;
}  

int SquareDet::GetMinUCell() const
{
  return m_minCellU;
}  

int SquareDet::GetMinVCell() const
{
  return m_minCellV;
}  


double SquareDet::GetSensitiveMaxU() const
{
  return m_sensitiveSizeU/2.; 
}  
  
double SquareDet::GetSensitiveMaxV() const
{
  return m_sensitiveSizeV/2.;  
} 

double SquareDet::GetSensitiveMinU() const
{
  return -m_sensitiveSizeU/2.; 
}  
  
double SquareDet::GetSensitiveMinV() const
{
  return -m_sensitiveSizeV/2.;  
} 

bool SquareDet::areNeighbors(int vcell1, int ucell1, int vcell2, int ucell2) const
{
  int deltaV = abs(vcell1-vcell2);
  int deltaU = abs(ucell1-ucell2);
 
  // Not sure what to do with this variable?
  // Maybe add a distance variable to arguments          
  int accept = 1; 
  
  // A side in common
  if(deltaU+deltaV < 2) return true;
     
  // A corner in common 
  if(accept == 1 && deltaU == 1 && deltaV == 1) return true;
     
  // max distance is 2 pixels (includes cases with missing pixels)
  if(accept == 2 && deltaU+deltaV < 3 ) return true;
     
  // max distance is 2 pixels along a diagonal (includes cases with missing pixels)
  if(accept == 3 && deltaU < 3 && deltaV < 3) return true;
                        
  return false;
}

int SquareDet::GetPixelTypeU(int ucell) const
{
  int i = 0; 
  for (auto group : m_uCells ) {
    int minCell = std::get<0>(group);
    int maxCell = std::get<1>(group);
    if (ucell >= minCell && ucell <= maxCell) 
      break;     
    i++;
  }
  return i; 
}

int SquareDet::GetPixelTypeV(int vcell) const 
{
  int i = 0; 
  for (auto group : m_vCells ) {
    int minCell = std::get<0>(group);
    int maxCell = std::get<1>(group);
    if (vcell >= minCell && vcell <= maxCell) 
      break;     
    i++;
  }
  return i; 
} 

int SquareDet::GetPixelType(int vcell, int ucell) const   
{ 
  int iu = GetPixelTypeU(ucell); 
  int iv = GetPixelTypeV(vcell); 
  int nGroupsU = m_uCells.size();
  return (nGroupsU*iv + iu);
}

double SquareDet::GetPitchU(int /*vcell*/, int ucell) const
{
  auto group = m_uCells.at(GetPixelTypeU(ucell));   
  return std::get<2>(group); 
} 
  
double SquareDet::GetPitchV(int vcell, int /*ucell*/) const
{
  auto group = m_vCells.at(GetPixelTypeV(vcell));   
  return std::get<2>(group); 
}  

int SquareDet::encodePixelID(int vcell, int ucell) const
{
  return (m_nCellsU*vcell + ucell);
}


void SquareDet::decodePixelID(int& vcell, int& ucell, int uniqPixelID) const
{
  vcell = uniqPixelID / m_nCellsU;
  ucell = uniqPixelID - vcell*m_nCellsU;
}
 
 	

bool SquareDet::SensitiveCrossed(double u, double v, double w) const
{
  if (u < -(m_sensitiveSizeU)/2.  || u > (m_sensitiveSizeU)/2.) {
   return false;
  }
  if (v < -(m_sensitiveSizeV)/2. || v > (m_sensitiveSizeV)/2.) {
    return false;
  }
  if (w < -m_sensitiveThickness/2. || w > m_sensitiveThickness/2.) {
    return false;
  }
  return true; 
}


bool SquareDet::isPointOutOfSensor( double u, double v, double w) const
{
  bool isOut = false; 
  
  // Boundary set +- epsilon
  if ( (u < (-m_sensitiveSizeU/2. -0.005)) || (u > (+m_sensitiveSizeU/2. + 0.005)) ||
       (v < (-m_sensitiveSizeV/2. -0.005)) || (v > (+m_sensitiveSizeV/2. + 0.005)) ||
       (w < (-m_sensitiveThickness/2. -0.005)) || (w > (+m_sensitiveThickness/2. + 0.005)) ) isOut = true;

   // Return if out or not
   return isOut;
}
 	

bool SquareDet::ModuleCrossed(double u, double v) const
{  
  if (u < -m_ladderSizeU/2.  || u > m_ladderSizeU/2.) {
   return false;
  }
  if (v < -m_ladderSizeV/2. || v > m_ladderSizeV/2.) {
    return false;
  }
  
  return true; 
}

double SquareDet::GetThickness(double u, double v) const
{
  if ( SensitiveCrossed(u, v) ) {
    return m_sensitiveThickness; 
  } 
  if ( ModuleCrossed(u, v) ) {
    return m_ladderThickness; 
  }
  return 0; 
}  


double SquareDet::GetTrackLength(double u, double v, double dudw, double dvdw) const
{
  return GetThickness(u,v)*std::sqrt(1 + dudw*dudw + dvdw*dvdw);  
}
    
double SquareDet::GetRadLength(double u, double v) const
{ 
  if ( SensitiveCrossed(u, v) ) {
    return m_sensitiveRadLength; 
  } 
  if ( ModuleCrossed(u, v) ) {
    return m_ladderRadLength; 
  }
  return materialeffect::X0_air; 
} 

double SquareDet::GetAtomicNumber(double u, double v) const
{ 
  if ( SensitiveCrossed(u, v) ) {
    return m_sensitiveAtomicNumber; 
  } 
  if ( ModuleCrossed(u, v) ) {
    return m_ladderAtomicNumber; 
  }
  return materialeffect::AtomicNumber_air;
} 

double SquareDet::GetAtomicMass(double u, double v) const
{ 
  if ( SensitiveCrossed(u, v) ) {
    return m_sensitiveAtomicMass;  
  } 
  if ( ModuleCrossed(u, v) ) {
    return m_ladderAtomicMass; 
  }
  return materialeffect::AtomicMass_air;
} 


double SquareDet::GetPixelCenterCoordV(int vcell, int /*ucell*/) const
{    
  int i = GetPixelTypeV(vcell);   
  double offset = m_offsetsV.at(i);
  auto group = m_vCells.at(i);
  int minCell = std::get<0>(group);  
  double pitch = std::get<2>(group);  
  
  // V coord measured from lower feft corner
  double v_coord = offset + pitch*(vcell - minCell) + pitch*0.5;
  // Ok, shift coord to sensor center
  v_coord -= 0.5*m_sensitiveSizeV;           
  
  return v_coord;
}
 

double SquareDet::GetPixelCenterCoordU(int /*vcell*/, int ucell) const
{
  int i = GetPixelTypeU(ucell);   
  double offset = m_offsetsU.at(i);
  auto group = m_uCells.at(i);
  int minCell = std::get<0>(group);  
  double pitch = std::get<2>(group);  
  
  // U coord measured from lower feft corner
  double u_coord = offset + pitch*(ucell - minCell) + pitch*0.5;
  // Ok, shift coord to sensor center
  u_coord -= 0.5*m_sensitiveSizeU;           
  
  return u_coord;
}

     
int SquareDet::GetUCellFromCoord( double u, double /*v*/ ) const
{
  if (u < -m_sensitiveSizeU/2.) {
   return m_minCellU;
  } else if (u > m_sensitiveSizeU/2.) {
   return m_nCellsU-1;
  } 
  
  int ucell = m_minCellU;
  double offset = -0.5*m_sensitiveSizeU;
  for (auto group : m_uCells ) {
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
   
   
int SquareDet::GetVCellFromCoord( double /*u*/, double v ) const 
{
  if (v < -m_sensitiveSizeV/2.) {
   return m_minCellV;
  } else if (v > m_sensitiveSizeV/2.) {
   return m_nCellsV-1;
  } 
  
  int vcell = m_minCellV;
  double offset = -0.5*m_sensitiveSizeV;
  for (auto group : m_vCells ) {
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

} // Namespace;

