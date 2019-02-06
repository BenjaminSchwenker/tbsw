#ifndef DET_H
#define DET_H 1

// Include basic C
#include <string>
#include <vector>
#include <tuple>

// Include TBTools 
#include "ReferenceFrame.h"

namespace depfet {
	
/** Class Det
 *  
 *  A pixel detector module
 *  
 *  The Det class represents a pixel module in a pixel tracking telescope. Its 
 *  main purpose is to provide a functional interface to the detector data 
 *  for reconstruction algorithms. The class interfaces all data fields for a 
 *  gear layer in the gear file. Most fields are self explanatory. Some need a
 *  little bit of explanation
 * 
 *  ID: Unique index to match geometry information to detector raw data 
 *  UCellID: Index of readout cells along u axis, also called column. Column i neignbors columns i-1 and i+1  
 *  VCellID: Index of readout cells along v axis, also called rows. Row i neignbors rows i-1 and i+1  
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
 *  @Author B. Schwenker, University of Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
 
class Det {
   	
 public:
   
  // Constructors 
  Det();
   
  // Get/Set functions to implement interface to gear layer 
  // ------------------------------------------------------
 
  // Get/Set detector DAQ ID
  int GetDAQID(){ return _ID; };
  void SetDAQID (int ID){ _ID = ID; };
  
  // Get/Set plane number   
  int GetPlaneNumber()  { return _PlaneNumber; } ;
  void SetPlaneNumber(int ipl)  { _PlaneNumber = ipl; } ;

  // Get/Set thickness in sensitive volume
  double GetSensitiveThickness()  { return _SensitiveThickness; };
  void SetSensitiveThickness(double thickness)  { _SensitiveThickness = thickness; };

  // Get/Set radiation length in sensitive volume
  double GetSensitiveRadLength() { return _SensitiveRadLength; };
  void SetSensitiveRadLength(double X0) { _SensitiveRadLength=X0; };

  // Get/Set atomic number in sensitive volume 
  double GetSensitiveAtomicNumber() { return _SensitiveAtomicNumber; };
  void SetSensitiveAtomicNumber(double Z) { _SensitiveAtomicNumber=Z; };

  // Get/Set atomic mass in sensitive volume 
  double GetSensitiveAtomicMass() { return _SensitiveAtomicMass; };
  void SetSensitiveAtomicMass(double A) { _SensitiveAtomicMass=A; };

  // Get/Set size u of ladder (sensitive + supports)  
  double GetLadderSizeU()  { return _LadderSizeU; };
  void SetLadderSizeU(double size)  { _LadderSizeU = size; };
   
  // Get/Set size V of ladder (sensitive + supports)  
  double GetLadderSizeV()  { return _LadderSizeV; };
  void SetLadderSizeV(double size)  { _LadderSizeV = size; };

  // Get/Set thickness of ladder in nonsentive volume    
  double GetLadderThickness() { return _LadderThickness; };
  void SetLadderThickness(double thickness) { _LadderThickness=thickness; };
  
  // Get/Set radiation length in nonsentive volume    
  double GetLadderRadLength()  { return _LadderRadLength; };
  void SetLadderRadLength(double X0)  { _LadderRadLength = X0; };
   
  // Get/Set atomic number in nonsensitive volume 
  double GetLadderAtomicNumber() { return _LadderAtomicNumber; };
  void SetLadderAtomicNumber(double Z) { _LadderAtomicNumber=Z; };

  // Get/Set atomic mass in nonsensitive volume 
  double GetLadderAtomicMass() { return _LadderAtomicMass; };
  void SetLadderAtomicMass(double A) { _LadderAtomicMass=A; };
  
  // Set the u cells  
  void SetCellsU( std::vector< std::tuple<int,int,double> > uCells); 
  
  // Set the v cells  
  void SetCellsV( std::vector< std::tuple<int,int,double> > vCells); 
       
  // Get/Set nominal sensor frame (i.e. where the detector is supposed to be)
  ReferenceFrame & GetNominal()  { return _Nominal; };
  void SetNominalFrame(ReferenceFrame nominal) { _Nominal  = nominal; };   
   	
  // Get/Set discrete rotation (mounting orientation of sensor) 
  ReferenceFrame & GetDiscrete()  { return _Discrete; }; 
  void SetDiscreteFrame(ReferenceFrame discrete) {_Discrete  = discrete;};   
 
  // These functions are relevant for tracking and DUT analysis
  // ----------------------------------------------------------
  // Their implementation will typically depend on the actual 
  // sensor design. 
 
  // Print methods 
  virtual void Print();

  // Get pixel type u axis
  virtual int GetPixelTypeU(int vcell, int ucell);  
   
  // Get pixel type v axis
  virtual int GetPixelTypeV(int vcell, int ucell);  

  // Get pixel type
  virtual int GetPixelType(int vcell, int ucell);  

  // Get number of ucells  
  virtual int GetNCellsU();  

  // Get number of vcells  
  virtual int GetNCellsV();

  // Get number of ucells  
  virtual int GetNColumns() { return _nCellsU; }  

  // Get number of vcells  
  virtual int GetNRows() { return _nCellsV; }

  // Get size u of active area  
  virtual double GetSensitiveSizeU();  
  
  // Get size v of active area  
  virtual double GetSensitiveSizeV(); 
      
  // Get column/u pixel pitch 
  virtual double GetPitchU(int vcell, int ucell); 
  
  // Get row/v pixel pitch  
  virtual double GetPitchV(int vcell, int ucell);  
   
  // Check if module box crossed
  virtual bool ModuleCrossed(double u, double v, double w = 0);

  virtual bool isPointOutOfSensor( double u , double v , double w = 0); 
  
  // Check if sensitive volume crossed
  virtual bool SensitiveCrossed(double u, double v, double w = 0);
      
  // Get thickness at coordinates (u,v) 
  virtual double GetThickness(double u, double v);

  // Get length of particle track (mm)
  virtual double GetTrackLength(double u, double v, double dudw, double dvdw);  
    
  // Get radlenght at coordinates (u,v) 
  virtual double GetRadLength(double u, double v);

  // Get atomic number at coordinates (u,v)
  virtual double GetAtomicNumber(double u, double v); 
 
  // Get atomic mass at coordinates (u,v)
  virtual double GetAtomicMass(double u, double v);   
  
  // Calculate pixel column from coord (u,v)   
  virtual int GetUCellFromCoord( double u, double v );
  
  // Calculate pixel vcell from coord (u,v)  
  virtual int GetVCellFromCoord( double u, double v ); 
  
  // Get u coord of pixel center 
  virtual double GetPixelCenterCoordU(int vcell, int ucell);
  
  // Get v coord of pixel center 
  virtual double GetPixelCenterCoordV(int vcell, int ucell);

  // Encode pixelID
  virtual int encodePixelID(int vcell, int ucell);

  // Decode pixelID
  virtual void decodePixelID(int vcell, int ucell, int uniqPixelID);
 
 private:
   
  // ID
  int _ID;
  // Plane number along beam line 
  int _PlaneNumber; 
  // Cells along sensor u axis  
  std::vector< std::tuple<int,int,double> > _uCells;
  // Cells along sensor v axis   
  std::vector< std::tuple<int,int,double> > _vCells ; 
  // Thickness in sensitive volume
  double _SensitiveThickness;
  // Rad. length in sensitive volume
  double _SensitiveRadLength;
  // Atomic number in sensitive volume
  double _SensitiveAtomicNumber;  
  // Atomic mass in sensitive volume
  double _SensitiveAtomicMass; 
  // Size along u of sensitive + supports  
  double _LadderSizeU; 
  // Size along v of sensitive + supports 
  double _LadderSizeV; 
  // Thickness in nonsensitive volume
  double _LadderThickness;
  // Rad. length in nonsensitive volume
  double _LadderRadLength;
  // Atomic number in nonsensitive volume
  double _LadderAtomicNumber;  
  // Atomic mass in nonsensitive volume
  double _LadderAtomicMass; 
       
  // Nominal reference frame
  ReferenceFrame _Nominal;
  // Mounting rotations  
  ReferenceFrame _Discrete; 

  // Some helper variables
  int _nCellsU; 
  int _nCellsV; 
  int _minCellU;
  int _minCellV; 
  double _sensitiveSizeU;  
  double _sensitiveSizeV; 
  std::vector<double> _offsetsU;
  std::vector<double> _offsetsV; 
};
 
} // Namespace

#endif // DET_H
