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

#include <string>
#include <vector>
#include <tuple>
#include <limits>

// Include ROOT
#include <TH2Poly.h>
#include <TGraph.h>
#include <TList.h>

// Namespaces
using namespace marlin;

namespace depfet {

//TODO sensitive size as parameter, there could be a max size from layout but only comparing all points
PolyDet::PolyDet(const std::string& typeName, int sensorID, int planeNumber, 
                 double sensThick, double sensRadLenght, double sensAtomicNumber,
                 double sensAtomicMass, double ladderThick, double ladderRadLength, 
                 double ladderAtomicNumber, double ladderAtomicMass, double ladderSizeU, 
                 double ladderSizeV, const std::vector< std::tuple<int,int,int,double,double> >& cells, 
		         const std::vector< std::tuple<int,double,double,std::vector<std::tuple<double,double>>> > & protocells,
                 const ReferenceFrame& discrete, const ReferenceFrame& nominal )
  : Det(typeName, sensorID, planeNumber) 
{
  
  // Set cells and layout
  SetCells(cells, protocells);

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


// set during TH2Poly creation
int PolyDet::GetMaxUCell()
{
  return m_maxCellU;
}  

int PolyDet::GetMaxVCell()
{
  return m_maxCellV;
}  

// Creates the TH2Poly layout with the center at 0,0
void PolyDet::SetCells(const std::vector< std::tuple<int, int, int, double, double> >& cells, const std::vector< std::tuple<int,double,double,std::vector<std::tuple<double,double>>>> & protocells)
{
  m_cells = cells;
  
  for (auto group: protocells){ // not checking if type already in ?
    m_cells_neighb_dist.emplace_back(std::get<0>(group), std::get<1>(group), std::get<2>(group));
  }
  TH2Poly *unshiftedLayout = new TH2Poly();
  unshiftedLayout->SetFloat();
  unshiftedLayout->SetName("helperlayout");

  m_layout = new TH2Poly();
  m_layout->SetFloat();
  m_layout->SetStats(0);
  m_layout->SetName((std::string("layout_d"+std::to_string(GetSensorID()))).c_str());

  double sensSizeUmax = -std::numeric_limits<double>::max();
  double sensSizeUmin = std::numeric_limits<double>::max();
  double sensSizeVmax = -std::numeric_limits<double>::max();
  double sensSizeVmin = std::numeric_limits<double>::max();

  m_maxCellU = 0;
  m_maxCellV = 0;
  m_minCellU = std::numeric_limits<int>::max();
  m_minCellV = std::numeric_limits<int>::max();
  // generate the pixel in the layout from the center position and the prototype pixel.
  for (auto group: cells){
    int type = std::get<2>(group);
    double centeru = std::get<3>(group);
    double centerv = std::get<4>(group);
    for (auto protopix: protocells){
      if (type == std::get<0>(protopix)){
        int i = 0;
        TGraph *gpixel = new TGraph(std::get<3>(protopix).size()); 
	std::string pixelname = std::to_string(std::get<0>(group)) + "," + std::to_string(std::get<1>(group));
        gpixel->SetName(pixelname.c_str());
	m_maxCellU = (std::get<0>(group) > m_maxCellU) ? std::get<0>(group) : m_maxCellU;
        m_maxCellV = (std::get<1>(group) > m_maxCellV) ? std::get<1>(group) : m_maxCellV;
        m_minCellU = (std::get<0>(group) < m_minCellU) ? std::get<0>(group) : m_minCellU;
        m_minCellV = (std::get<1>(group) < m_minCellV) ? std::get<1>(group) : m_minCellV;
        for (auto points: std::get<3>(protopix)){
          double x = std::get<0>(points)+centeru;
	  double y = std::get<1>(points)+centerv;
          gpixel->SetPoint(i, x, y);
	  sensSizeUmax = (x > sensSizeUmax) ? x : sensSizeUmax;
	  sensSizeUmin = (x < sensSizeUmin) ? x : sensSizeUmin;
	  sensSizeVmax = (y > sensSizeVmax) ? y : sensSizeVmax;
	  sensSizeVmin = (y < sensSizeVmin) ? y : sensSizeVmin;
	  i++;
	}
        unshiftedLayout->AddBin(gpixel);
        break;
      }
    } 
  }
  // layout not necessary centred around origin, so shift it
  double shiftu = sensSizeUmax - (sensSizeUmax-sensSizeUmin)/2.;
  double shiftv = sensSizeVmax - (sensSizeVmax-sensSizeVmin)/2.;
  
  TList *binList = unshiftedLayout->GetBins();
  TH2PolyBin *polyBin;
  TIter next(binList);
  TObject *obj = 0;
  while ((obj = next())){
    polyBin = (TH2PolyBin*)obj;
    TGraph *gpoly = (TGraph*)polyBin->GetPolygon();
    for (int k = 0; k < gpoly->GetN(); k++){
      double x = 0.0;
      double y = 0.0;
      gpoly->GetPoint(k, x, y);
      gpoly->SetPoint(k, x-shiftu, y-shiftv);
    }
    m_layout->AddBin(gpoly->Clone());    
  }
  // sensitive area
  m_sensitiveSizeU = sensSizeUmax - sensSizeUmin;
  m_sensitiveSizeV = sensSizeVmax - sensSizeVmin;

  // number of u cells and v cells, assuemed u,v start at min then number is max-minu/v +1
  m_nCellsU = m_maxCellU - m_minCellU + 1;
  m_nCellsV = m_maxCellV - m_minCellV + 1;

  delete unshiftedLayout;
  unshiftedLayout = nullptr;

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
  GetPixelCenterCoord(vcoord1, ucoord1, vcell1, ucell1);
  GetPixelCenterCoord(vcoord2, ucoord2, vcell2, ucell2);
  
  int type1 = GetPixelType(vcell1, ucell1);
  int type2 = GetPixelType(vcell2, ucell2); 

  double distu1 = 0.;
  double distv1 = 0.;
  double distu2 = 0.;
  double distv2 = 0.;
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
  if ((abs(ucoord1-ucoord2) < distu1 || abs(ucoord1-ucoord2) < distu2) && (abs(vcoord1-vcoord2) < distv1 || abs(vcoord1-vcoord2) < distv2))
    return true;

  return false;
}

int PolyDet::GetPixelType(int vcell, int ucell)  
{ 
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
    if (type == std::get<0>(group))
      return std::get<1>(group);
  }  
  return -1.0; // better return value?
} 

double PolyDet::GetPitchV(int vcell, int ucell)
{
  int type = GetPixelType(vcell, ucell);
  for (auto group: m_pitch){
    if (type == std::get<0>(group))
      return std::get<2>(group);
  }  
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

bool PolyDet::SensitiveCrossed(double u, double v, double w)
{
  int bin = m_layout->FindBin(u, v);
  if (bin < 0 && bin != -5)
    return false;
  if (w < -m_sensitiveThickness/2. || w > m_sensitiveThickness/2.) {
    return false;
  }
  return true; 
}

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
  return -m_sensitiveSizeU/2.-0.1;
}


void PolyDet::GetPixelCenterCoord(double& vcoord, double& ucoord, int vcell, int ucell)
{
  vcoord = GetPixelCenterCoordV(vcell, ucell);
  ucoord = GetPixelCenterCoordU(vcell, ucell);
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

