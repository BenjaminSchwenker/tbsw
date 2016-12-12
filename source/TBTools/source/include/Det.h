#ifndef DET_H
#define DET_H 1

// Include basic C
#include <string>

// Include CLHEP classes
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/Vector.h>

// DEPFETTrackTools include
#include "ReferenceFrame.h"

namespace depfet {
	
/** Class Det
 *  
 *  A pixel detector module in space.
 *  
 *  The Det class represents a pixel module in a pixel tracking telescope. Its 
 *  main purpose is to provide a functional interface to the detector data 
 *  for reconstruction algorithms. The class interfaces all data fields for a 
 *  gear layer in the gear file. Most fields are self explanatory. Some need a
 *  little bit of explanation
 * 
 *  a) DAQ ID: Unique index to match detector to event data
 *  b) Device Type: Unique index to represent sensor layout      
 *   
 *  The position and orientation of the detector is represented by an object
 *  of class ReferenceFrame called 'Nominal'. In short, the 'Nominal' reference 
 *  frame manages the transformation from local UVW to global XYZ coordinates. 
 *  The 'Nominal' reference frame is loaded from the alignment data base, if it
 *  is available. Otherwise, it is loaded from the gear file. 
 *  
 *  A special instance of class ReferenceFrame called 'Discrete' is used to  
 *  represent the discrete 'flips' from the local UVW to global XYZ axes. 
 *  The 'Discrete' reference frame represents the different possibolities 
 *  to install the detector in the beam line,i.e. is the detector front 
 *  side pointing into the beam? Or, is the local u axis pointing upwards? 
 *  
 *  The Det class also provides an advanced interface to detector data for 
 *  reconstruction methodes. This includes: 
 *    
 *  A) Transform between local coord and readout channels 
 *  B) Position resolved detector material profile data  
 *  C) Boundary intersection tests for partical tracking 
 *  
 *  IMPORTANT NOTE TO DEVELOPERS
 *  
 *  Some member functions may be overwritten in daugther classes. This can be 
 *  used by sensor developers to adapt the class behavior to new designs. Their 
 *  must be unique 'DeviceType' for each daugther class. There should also 
 *  be some mechanism to register new Det implementations -> TODO 
 *  
 *  @Author B. Schwenker, University of Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
 
class Det {
   	
 public:
   
  // Constructors 
  Det();
   
  // Get/Set functions to implement interface to gear layer 
  // ------------------------------------------------------
  
  // Get/Set detecor device type
  int GetDeviceType()  { return DeviceType; };
  void SetDeviceType(int type)  { DeviceType = type; };
  
  // Get/Set detector Name 
  std::string GetName() { return Name; };
  
  // Get/Set detector DAQ ID
  int GetDAQID(){ return DAQID; };
  void SetDAQID (int ID){ DAQID = ID; };
  
  // Get/Set plane number   
  int GetPlaneNumber()  { return PlaneNumber; } ;
  void SetPlaneNumber(int ipl)  { PlaneNumber = ipl; } ;
  
  //! Get/Set number of pixel columns  
  double GetNColumns() { return NColumns; };
  void SetNColumns(int cols) { NColumns = cols; };
    
  //! Get/Set number of pixel rows  
  double GetNRows() { return NRows; };
  void SetNRows(int rows )  { NRows = rows; };
     
  //! Get/Set column/u pixel pitch 
  double GetPitchU()  { return PitchU; };
  void SetPitchU(double pitch)  { PitchU = pitch; };
  
  //! Get/Set row/v pixel pitch  
  double GetPitchV()  { return PitchV; };	
  void SetPitchV(double pitch)  { PitchV=pitch; };
  
  // Get/Set thickness of active silicon 
  double GetSensitiveThickness()  { return SensitiveThickness; };
  void SetSensitiveThickness(double thickness)  { SensitiveThickness = thickness; };
  
  // Get/Set thickness of support module box    
  double GetLadderThickness() { return LadderThickness; };
  void SetLadderThickness(double thickness) { LadderThickness=thickness; };
  
  // Get/Set sensitive radiation length
  double GetSensitiveRadLength() { return SensitiveRadLength; };
  void SetSensitiveRadLength(double X0) { SensitiveRadLength=X0; };
  
  // Get/Set ladder radiation length
  double GetLadderRadLength()  { return LadderRadLength; };
  void SetLadderRadLength(double X0)  { LadderRadLength = X0; };
   
  // Get/Set resolution along columns, or U  
  double GetResolutionU()  { return ResolutionU; };
  void SetResolutionU(double res)  { ResolutionU = res; };
  
  // Get/Set resolution along rows, or V  
  double GetResolutionV()  { return ResolutionV; }; 
  void SetResolutionV(double res)  { ResolutionV=res; }; 
    
  // Get/Set size u of module box  
  double GetModuleBoxSizeU()  { return ModuleBoxSizeU; };
  void SetModuleBoxSizeU(double size)  { ModuleBoxSizeU = size; };
   
  // Get/Set size v of module box 
  double GetModuleBoxSizeV()  { return ModuleBoxSizeV; };
  void SetModuleBoxSizeV(double size)  { ModuleBoxSizeV = size; };
  
  // Get/Set size u of active area  
  virtual double GetSensitiveSizeU()  { return NColumns*PitchU; };
  void SetSensitiveSizeU(double size)  { SensitiveSizeU = size; };
  
  // Get/Set size v of active area  
  double GetSensitiveSizeV()  { return NRows*PitchV; };
  void SetSensitiveSizeV(double size)  { SensitiveSizeV = size; };	
   
  // Get/Set nominal sensor frame (i.e. where the detector is supposed to be)
  ReferenceFrame & GetNominal()  { return Nominal; };
  void SetNominalFrame(ReferenceFrame nominal) { Nominal  = nominal; };   
   	
  // Get/Set discrete rotation (mounting orientation of sensor) 
  ReferenceFrame & GetDiscrete()  { return Discrete; }; 
  void SetDiscreteFrame(ReferenceFrame discrete) {Discrete  = discrete;};   
 
  // These functions are relevant for tracking and DUT analysis
  // ----------------------------------------------------------
  // Their implementation will typically depend on the actual 
  // sensor design. These functions may be overwritten depending
  // in the device type. 
 
  // Print methods 
  virtual void Print();
   
  // Check if module box crossed
  virtual bool ModuleCrossed(double u, double v, double w = 0);

  bool isPointOutOfSensor( double u , double v , double w = 0); 
  
  // Check if sensitive volume crossed
  virtual bool SensitiveCrossed(double u, double v, double w = 0);
      
  // Get thickness at coordinates (u,v) 
  virtual double GetThickness(double u, double v);

  // Get length of particle track (mm)
  virtual double GetTrackLength(double u, double v, double dudw, double dvdw);  
    
  // Get radlenght at coordinates (u,v) 
  virtual double GetRadLength(double u, double v);
  
  // Calculate pixel column from coord (u,v)   
  virtual int GetColumnFromCoord( double u, double v );
  
  // Calculate pixel row from coord (u,v)  
  virtual int GetRowFromCoord( double u, double v ); 
  
  // Get u coord of pixel center 
  virtual double GetPixelCenterCoordU(int row, int column);
  
  // Get v coord of pixel center 
  virtual double GetPixelCenterCoordV(int row, int column);

  // Encode pixelID
  virtual int encodePixelID(int row, int col);

  // Decode pixelID
  virtual void decodePixelID(int & row, int & col, int uniqPixelID);
 
 private:
   
  // Name of module: 'DEPFET', 'TAKI', 'MIMOSA', etc.  
  std::string Name;
  // DAQ module "ID"
  int DAQID;
  // Plane number along beam line 
  int PlaneNumber; 
  // Type of sensor: '0' : non bricked; '1' : bricked pixels 
  int DeviceType;
  // Number of pixel columns  
  double NColumns;
  // Number of pixel rows 
  double NRows;
  // Column/U pixel pitch
  double PitchU;
  // Row/V pixel pitch  
  double PitchV;
  // Thickness in sensitive area   
  double SensitiveThickness;
  // Thickness outside active area
  double LadderThickness;
  // Rad. length in active area 
  double SensitiveRadLength;
  // Rad. length outside active area
  double LadderRadLength;
  // Module size along u  
  double ModuleBoxSizeU; 
  // Module size along v 
  double ModuleBoxSizeV; 
  // Sensitive size u
  double SensitiveSizeU; 
  // Sensitive size v 
  double SensitiveSizeV; 
  // Resolution in U    
  double ResolutionU;
  // Resolution in V  
  double ResolutionV;        
  // Nominal reference frame
  ReferenceFrame Nominal;
  // Mounting rotations  
  ReferenceFrame Discrete; 
  
};
 
} // Namespace

#endif // DET_H
