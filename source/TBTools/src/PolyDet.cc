// PolyDet implementation file 
// 
// Author: Helge C. Beck, University of Göttingen 
// <mailto:helge-christoph.beck@phys.uni-goettingen.de>

// TBTools includes 
#include "PolyDet.h"
#include "MaterialEffect.h"

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

#include <algorithm>

// Include ROOT
#include <TH2Poly.h>
#include <TGraph.h>

// Namespaces
using namespace marlin;

namespace depfet {

//TODO sensitive size as parameter, there could be a max size from layout but only comparing all points
PolyDet::PolyDet(const std::string& typeName, int sensorID, int planeNumber, 
                 double sensThick, double sensRadLenght, double sensAtomicNumber,
                 double sensAtomicMass, double ladderThick, double ladderRadLength, 
                 double ladderAtomicNumber, double ladderAtomicMass, double ladderSizeU, 
                 double ladderSizeV, const std::vector< std::tuple<int,int,int,double,double> >& cells, 
		 const std::vector< std::tuple<int,double,double,std::vector<std::tuple<double,double>>>> & protocells,
                 const ReferenceFrame& discrete, const ReferenceFrame& nominal )
  : Det(typeName, sensorID, planeNumber) 
{
  
  // Set cells and layout
  SetCells(cells, protocells)

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


PolyDet::PolyDet(const std::string& typeName, int sensorID, int planeNumber) : Det(typeName, sensorID, planeNumber) {}


PolyDet* PolyDet::newDet()
{
  return new PolyDet("PolyDet", -1, -1);
}
// set during TH2Poly creation
int PolyDet::GetMaxUCell()
{
  return m_nCellsU-1;
}  

int PolyDet::GetMaxVCell()
{
  return m_nCellsV-1;
}  

// TODO this code should be put into the class where the PolyDet get constructed from an XML file
void PolyDet::SetCells(const std::vector< std::tuple<int, int, int, double, double> >& cells, const std::vector< std::tuple<int,double,double,std::vector<std::tuple<double,double>>>> & protocells)
{
  m_cells = cells;
  
  for (auto group: protocells){ // not checking if type already in ?
    m_cells_neighb_dist.emplace_back(std::get<0>(group), std::get<1>(group), std::get<2>(group))
  }

  m_layout = new TH2Poly();
  m_layout->SetFloat();
  m_layout->SetStats(0);
  m_layout->SetName((std::string("layout_d"+std::to_string(GetPlaneID()))).c_str());
  // generate the bins in m_layout
  for (auto group: cells){
    TGraph gpixel; // does this need to be a pointer ? does this object need to stay to stay in the layout?
    std::string pixelname = std::to_string(std::get<0>(group)) + "," + std::to_string(std::get<1>(group));
    gpixel.SetName(pixelname.c_str());

    int type = std::get<3>(group);
    double centeru = std::get<4>(group);
    double centerv = std::get<5>(group);
    for (auto protopix: protocells){
      if (type == std::get<0>(protopix)){
        int i = 0;
        for (auto points: std::get<3>(protopix)){
          gpixel.SetPoint(i, std::get<0>(points)+centeru, std::get<1>(points)+centerv);
	  i++;
	}
        break;
      }
    }
    m_layout->AddBin(&gpixel);
  }
  // sensitive area?
}

double PolyDet::GetSensitiveSizeU()
{
  return m_sensitiveSizeU; 
}  
  
double PolyDet::GetSensitiveSizeV()
{
  return m_sensitiveSizeV;  
} 
// TODO different distance comparison?
bool PolyDet::areNeighbors(int vcell1, int ucell1, int vcell2, int ucell2)
{
  // get coord, type, make distance, compare to type neighbour distance
  // initialised with maximum distance, so not automatic neighbours as if assigned 0.
  double vcoord1 = -m_sensitiveSizeV/2.;
  double ucoord1 = -m_sensitiveSizeU/2.;
  double vcoord2 = m_sensitiveSizeV/2.;
  double ucoord2 = m_sensitiveSizeU/2.;
  GetPixelCenterCoodinate(vcoord1, ucoord1, vcell1, ucell1);
  GetPixelCenterCoodinate(vcoord2, ucoord2, vcell2, ucell2);
  
  int type1 = GetPixelType(vcell1, ucell1);
  int type2 = GetPixelType(vcell2, ucell2); 

  distu1 = 0.;
  distv1 = 0.;
  distu2 = 0.;
  distv2 = 0.;
  bool found1 = false;
  bool found2 = false;
  for (auto group: m_cells_neighb_dist){
    if (!found1 && type1 == std::get<0>(group)){
      distu1 = std::get<1>(group);
      distv1 = std::get<2>(group);
      found1 = true;
    }
    if (!found2 && type2 == std::get<0>(group)){
      distu2 = std::get<1>(group);
      distv2 = std::get<2>(group);
      found2 = true;
    }
    if (found1 && found2)
      break;
  }
  if ((abs(coordu1-coordu2) < distu1 || abs(coordu1-coordu2) < distu2) && (abs(coordv1-coordv2) < distv1 || abs(coordv1-coordv2) < distv2))
    return true;

  return false;
}

int PolyDet::GetPixelType(int vcell, int ucell)  
{
  int i = 0; 
  for (auto group : m_cells ) {
    int uCell = std::get<0>(group);
    int vCell = std::get<1>(group);
    if (ucell == uCell && vcell == vCell)
      return std::get<3>(group);
  }
  return -1; // not pixel match found?
} 

double PolyDet::GetPitchU(int vcell, int ucell)  
{
  int type = GetPixelType(vcell, ucell);
  for (auto group: m_pitch){
    if (type = std::get<0>(group))
      return std::get<1>(group);  
  return -1.0; // better return value?
} 

double PolyDet::GetPitchV(int vcell, int ucell)
{
  int type = GetPixelType(vcell, ucell);
  for (auto group: m_pitch){
    if (type = std::get<0>(group))
      return std::get<2>(group);  
  return -1.0;
}  

int PolyDet::encodePixelID(int vcell, int ucell)
{
  return (m_nCellsU*vcell + ucell);
}


void PolyDet::decodePixelID(int& vcell, int& ucell, int uniqPixelID)
{
  vcell = uniqPixelID / m_nCellsU;
  ucell = uniqPixelID - vcell*m_nCellsU;
}


// TODO would be faster with sensitiveSize, no need to search for bin?
bool PolyDet::SensitiveCrossed(double u, double v, double w)
{
  // use TH2Poly isInside(u,v)?
  // catch -5?
  int bin m_layout->FindBin(u, v);
  if (bin < 0 && bin != -5)
    return false;
  /*if (u < -(m_sensitiveSizeU)/2.  || u > (m_sensitiveSizeU)/2.) {
   return false;
  }
  if (v < -(m_sensitiveSizeV)/2. || v > (m_sensitiveSizeV)/2.) {
    return false;
  }
  */
  if (w < -m_sensitiveThickness/2. || w > m_sensitiveThickness/2.) {
    return false;
  }
  return true; 
}

// TODO leave as is, difficult with margin and findBin
bool PolyDet::isPointOutOfSensor( double u, double v, double w) 
{
  bool isOut = false; 
  
  // Boundary set +- epsilon
  if ( (u < (-m_sensitiveSizeU/2. -0.005)) || (u > (+m_sensitiveSizeU/2. + 0.005)) ||
       (v < (-m_sensitiveSizeV/2. -0.005)) || (v > (+m_sensitiveSizeV/2. + 0.005)) ||
       (w < (-m_sensitiveThickness/2. -0.005)) || (w > (+m_sensitiveThickness/2. + 0.005)) ) isOut = true;

   // Return if out or not
   return isOut;
}
	
// TODO leave as is as ladder not described with TH2Poly
bool PolyDet::ModuleCrossed(double u, double v)
{  
  if (u < -m_ladderSizeU/2.  || u > m_ladderSizeU/2.) {
   return false;
  }
  if (v < -m_ladderSizeV/2. || v > m_ladderSizeV/2.) {
    return false;
  }
  
  return true; 
}

double PolyDet::GetThickness(double u, double v)
{
  if ( SensitiveCrossed(u, v) ) {
    return m_sensitiveThickness; 
  } 
  if ( ModuleCrossed(u, v) ) {
    return m_ladderThickness; 
  }
  return 0; 
}  


double PolyDet::GetTrackLength(double u, double v, double dudw, double dvdw)
{
  return GetThickness(u,v)*std::sqrt(1 + dudw*dudw + dvdw*dvdw);  
}
    
double PolyDet::GetRadLength(double u, double v)
{ 
  if ( SensitiveCrossed(u, v) ) {
    return m_sensitiveRadLength; 
  } 
  if ( ModuleCrossed(u, v) ) {
    return m_ladderRadLength; 
  }
  return materialeffect::X0_air; 
} 

double PolyDet::GetAtomicNumber(double u, double v)
{ 
  if ( SensitiveCrossed(u, v) ) {
    return m_sensitiveAtomicNumber; 
  } 
  if ( ModuleCrossed(u, v) ) {
    return m_ladderAtomicNumber; 
  }
  return materialeffect::AtomicNumber_air;
} 

double PolyDet::GetAtomicMass(double u, double v)
{ 
  if ( SensitiveCrossed(u, v) ) {
    return m_sensitiveAtomicMass;  
  } 
  if ( ModuleCrossed(u, v) ) {
    return m_ladderAtomicMass; 
  }
  return materialeffect::AtomicMass_air;
} 


double PolyDet::GetPixelCenterCoordV(int vcell, int ucell)
{    
  // geometric centre of bounding box or charge collection centre as in config file specified?
  for (auto group: m_cells){
    if (ucell == std::get<0>(group) && vcell == std::get<1>(group))
      return std::get<4>(group);
  }
  return -m_sensitiveSizeV/2.-0.1; // -1.0 (mm) could be a actual position, so messing this up with this, better something out of sensitive range?! 100 um out off sensitive?
}
 

double PolyDet::GetPixelCenterCoordU(int vcell, int ucell)
{
  for (auto group: m_cells){
    if (ucell == std::get<0>(group) && vcell == std::get<1>(group))
      return std::get<3>(group);
  }
  return -m_sensistiveSizeU/2.-0.1;
}


void PolyDet::GetPixelCenterCoord(double& vcoord, double& ucoord, int vcell, int ucell)
{
  vcoord = GetPixelCenterCoordV(int vcell, int ucell);
  ucoord = GetPixelCenterCoordU(int vcell, int ucell);
}

int PolyDet::GetUCellFromCoord( double u, double v )
{
  if (u < -m_sensitiveSizeU/2.) {
   return -1;//m_minCellU;
  } else if (u > m_sensitiveSizeU/2.) {
   return -1;//m_nCellsU-1;
  } 
  // ? |
  int bin =  m_layout->FindBin(u, v);
  if (bin < 0) // possible non desribed areas in the sensitive volumne
    return -1;
  std::string binname = m_layout->GetBinName(bin);
  return std::stoi(binname.substr(0, binname.find(",")));
} 
   
   
int PolyDet::GetVCellFromCoord( double u, double v ) 
{
  if (v < -m_sensitiveSizeV/2.) {
   return -1;//m_minCellV; no minCell?
  } else if (v > m_sensitiveSizeV/2.) {
   return -1;//m_nCellsV-1;
  }
  // ? | 
  int bin =  m_layout->FindBin(u, v);
  if (bin < 0) // possible non desribed areas in the sensitive volumne
    return -1;
  std::string binname = m_layout->GetBinName(bin);
  return std::stoi(binname.substr(binname.find(",")+1));
}

} // Namespace;

