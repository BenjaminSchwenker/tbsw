// Det implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// DEPFETTrackTools includes 
#include "Det.h"
#include "MaterialEffect.h"

// Include basic C header files

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

// Namespaces
using namespace CLHEP;
using namespace marlin;

namespace depfet {

// Define constants
#define EPS 0.005           // mm

/** Dummy constructor for pixel modules
 *
 * Creates a dummy pixel module. 
 */
Det::Det() {Name = "DUMMY";}



int Det::encodePixelID(int row, int column)
{
  return (GetNColumns()*row + column);
}


void Det::decodePixelID(int & row, int & column, int uniqPixelID)
{
  row    = uniqPixelID / GetNColumns();
  column = uniqPixelID - row*GetNColumns();
}
 

 	
/**  Check if sensitive volume is crossed
 */
bool Det::SensitiveCrossed(double u, double v, double w)
{

  // Smear the active area boundaries by half pixel pitch 

  if (u < -(GetSensitiveSizeU()+GetPitchU())/2.  || u > (GetSensitiveSizeU()+GetPitchU())/2.) {
   return false;
  }
  if (v < -(GetSensitiveSizeV()+GetPitchV())/2. || v > (GetSensitiveSizeV()+GetPitchV())/2.) {
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
  if ( (u < (-GetSensitiveSizeU()/2. -EPS)) || (u > (+GetSensitiveSizeU()/2. + EPS)) ||
       (v < (-GetSensitiveSizeV()/2. -EPS)) || (v > (+GetSensitiveSizeV()/2. + EPS)) ||
       (w < (-GetSensitiveThickness()/2. -EPS)) || (w > (+GetSensitiveThickness()/2. + EPS)) ) isOut = true;

   // Return if out or not
   return isOut;
}
 	
/**  Check if module box is crossed
 */
bool Det::ModuleCrossed(double u, double v, double w)
{  
  if (u < -GetModuleBoxSizeU()/2.  || u > GetModuleBoxSizeU()/2.) {
   return false;
  }
  if (v < -GetModuleBoxSizeV()/2. || v > GetModuleBoxSizeV()/2.) {
    return false;
  }
  if (w < -GetLadderThickness()/2. || w > GetLadderThickness()/2.) {
    return false;
  }
  return true; 
}

/**  Check if module box is crossed
 */
double Det::GetThickness(double u, double v)
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
double Det::GetTrackLength(double u, double v, double dudw, double dvdw)
{
  return GetThickness(u,v)*std::sqrt(1 + dudw*dudw + dvdw*dvdw);  
}
    
  
/**  Check if module box is crossed
 */
double Det::GetRadLength(double u, double v)
{ 
  if ( SensitiveCrossed(u, v) ) {
    return GetSensitiveRadLength(); 
  } 
  if ( ModuleCrossed(u, v) ) {
    return GetLadderRadLength(); 
  }
  return materialeffect::X0_air; 
} 

/**  Get v coord of pixel center 
 */
double Det::GetPixelCenterCoordV(int row, int column)
{  
  
  
  // Get pitch
  double sensPitch = GetPitchV();
  
  // V coord measured from lower feft corner
  double v_coord = sensPitch*(row + 0.5);
  // Ok, shift coord to sensor center
  v_coord -= 0.5*GetNRows()*sensPitch; 
  
  return v_coord;
}
 
/**  Get u coord of pixel center 
 */
double Det::GetPixelCenterCoordU(int row, int column)
{
   
  double u_coord;
  // Get pitch
  double sensPitch = GetPitchU();
  if ( GetDeviceType() == 0 ){
    // Non Bricked 
    u_coord = sensPitch*(column + 0.5);
  } 
  else {
    // Bricked pixels, every second row shifted
    if (row%2 == 0) {
      u_coord = sensPitch*column;
    }
    else {       
      u_coord = sensPitch*(column + 0.5);
    }
  } 
  
  // Ok, shift coord to sensor center
  u_coord -= 0.5*GetNColumns()*sensPitch; 
  
  return u_coord;
}

/**  Calculate pixel column from coord (u,v)
 */      
int Det::GetColumnFromCoord( double u, double v )
{
  
  int column = -1;
  double sensPitch = GetPitchU();
  int sensNPixels = GetNColumns();   
  
  if (sensPitch == 0) {
    streamlog_out(ERROR) << "GetColumnFromPoint - pitchU is zero!!!"
                         << std::endl;
  }
  
  // For convinience, measure point from lower left corner of sensor
  u += 0.5*GetNColumns()*GetPitchU();
  v += 0.5*GetNRows()*GetPitchV();
  
  if ( GetDeviceType() == 0 ){
    // Non Bricked 
    if (u <= 0.) column = -1;  // overflow
    else {
      column = floor( u / sensPitch );
      if (column >= sensNPixels) column = sensNPixels; // overflow
    }    
  } else {
    // Bricked Pixels, every second row shifted
    if ( GetRowFromCoord(u,v)%2==0) {
      if (u <= 0.) column = 0;
      else {
        column = floor( u / sensPitch );
        if (column >= sensNPixels) column = sensNPixels - 1;
      }    
    }    
    else {   
      // First pixel
      if (u <= sensPitch/2.) column = 0;
      else {
        column = floor( (u + sensPitch/2.) / sensPitch);      
        if (column >= sensNPixels) column = sensNPixels - 1;
      }
    }
  } 
  
  return column;  
} 
   
/**  Calculate pixel column from coord (u,v)
 */   
int Det::GetRowFromCoord( double u, double v ) 
{
  
  int row = -1;
  double sensPitch = GetPitchV();
  int sensNPixels = GetNRows(); 
  
  if (sensPitch == 0) {
    streamlog_out(ERROR) << "Det::GetPixelRow - pitchV is zero!!!"
                         << std::endl;
  }
   
  // For convinience, measure point from lower left corner of sensor
  u += 0.5*GetNColumns()*GetPitchU();
  v += 0.5*GetNRows()*GetPitchV();

  if (v <= 0.) row = -1; // overflow
  else {   
    row = floor( v / sensPitch );
    if (row >= sensNPixels) row = sensNPixels; // overflow
  }
  
  return row;
}


/**  Print methods
 */ 
void Det::Print()
{
  streamlog_out(MESSAGE3) << std::endl;
  streamlog_out(MESSAGE3) << "  Sensor Name:     " << GetName()           << std::endl;
  streamlog_out(MESSAGE3) << "  Plane Number:    " << GetPlaneNumber()    << std::endl;
  streamlog_out(MESSAGE3) << "  DAQ ID:          " << GetDAQID()          << std::endl;  
  streamlog_out(MESSAGE3) << "  PitchU[um]:      " << GetPitchU()*1000    << std::endl;
  streamlog_out(MESSAGE3) << "  PitchV[um]:      " << GetPitchV()*1000    << std::endl;
  streamlog_out(MESSAGE3) << "  NColumns:        " << GetNColumns()       << std::endl;
  streamlog_out(MESSAGE3) << "  NRows:           " << GetNRows()          << std::endl;  
}
 


} // Namespace;

