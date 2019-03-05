#ifndef POLYDET_H
#define POLYDET_H 1

// Include basic C
#include <vector>
#include <tuple>

// Include TBTools 
#include "Det.h"

// Include ROOT
#include <TH2Poly.h>
namespace depfet {
	
/** Class PolyDet
 *  TODO
 *  A pixel detector module
 *  
 *  The SquareDet class represents a pixel module in a pixel tracking telescope. Its 
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
 *  The SquareDet class also provides an advanced interface to detector data for 
 *  reconstruction methodes. This includes: 
 *    
 *  A) Transform between local coord and readout channels 
 *  B) Position resolved detector material profile data  
 *  C) Boundary intersection tests for partical tracking 
 *  
 *  @Author H. C. Beck, University of Göttingen
 *  <mailto:helge-christoph.beck@phys.uni-goettingen.de>
 */
 
class PolyDet : public Det {
   	
 public:
//  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
 
   
  /** Constructor
   */
  PolyDet(const std::string& typeName, int sensorID, int planeNumber) ; 
  
  /** Constructor
   */
  PolyDet::PolyDet(const std::string& typeName, int sensorID, int planeNumber, 
                 double sensThick, double sensRadLenght, double sensAtomicNumber,
                 double sensAtomicMass, double ladderThick, double ladderRadLength, 
                 double ladderAtomicNumber, double ladderAtomicMass, double ladderSizeU, 
                 double ladderSizeV, const std::vector< std::tuple<int,int,int,double,double> >& cells, 
		 const std::vector< std::tuple<int,double,double,std::vector<std::tuple<double,double>>>> & protocells,
                 const ReferenceFrame& discrete, const ReferenceFrame& nominal );

  /** Get pixel type for pixel at position vcell and ucell. 
   */
  int GetPixelType(int vcell, int ucell) override;  

  /** Get the maximum uCell on the pixel matrix.   
   *  The uCell numbers of pixels are in the intervall [0,max].
   */  
  int GetMaxUCell() override;  

  /** Get the maximum vCell on the pixel matrix.   
   *  The vCell numbers of pixels are in the intervall [0,max].
   */  
  int GetMaxVCell() override;  
  
  /** Get size u of sensitive area. 
   *  All pixels must be inside the box 
   *  [-SensSizeU/2,SensSizeU/2]x[-SensSizeV/2,SensSizeV/2].
   */  
  double GetSensitiveSizeU() override;  
  
  /** Get size v of sensitive area.  
   *  All pixels must be inside the box 
   *  [-SensSizeU/2,SensSizeU/2]x[-SensSizeV/2,SensSizeV/2].
   */  
  double GetSensitiveSizeV() override; 
      
  /** Get u pitch for pixel at position vcell,ucell.   
   */ 
  double GetPitchU(int vcell, int ucell) override; 
  
  /** Get v pitch for pixel at position vcell,ucell.   
   */ 
  double GetPitchV(int vcell, int ucell) override;  
   
  /** Returns true if point (u,v,w) is not inside the sensitive volume.  
   */ 
  bool isPointOutOfSensor( double u , double v , double w = 0) override; 
  
  /** Returns true if point (u,v,w) is inside the sensitive volume. 
   */ 
  bool SensitiveCrossed(double u, double v, double w = 0) override;
      
  /** Get thickness of detector at position (u,v). 
   */ 
  double GetThickness(double u, double v) override;

  /** Get length of particle track intersecting the detector.  
   */
  double GetTrackLength(double u, double v, double dudw, double dvdw) override;  
    
  /** Get radlenght at position (u,v).
   */ 
  double GetRadLength(double u, double v) override;

  /** Get atomic number at position (u,v).
   */
  double GetAtomicNumber(double u, double v) override; 
 
  /** Get atomic mass at position (u,v).
   */
  double GetAtomicMass(double u, double v) override;   
  
  /** Returns uCell of pixel at position (u,v). Returns -1 if there is no pixel.    
   */    
  int GetUCellFromCoord( double u, double v ) override;
  
  /** Returns vCell of pixel at position (u,v). Returns -1 if there is no pixel.  
   */    
  int GetVCellFromCoord( double u, double v ) override; 
  
  /** Returns u position of pixel center. 
   */
  double GetPixelCenterCoordU(int vcell, int ucell) override;
  
  /** Returns v position of pixel center.   
   */
  double GetPixelCenterCoordV(int vcell, int ucell) override;
  
  /** Return unique ID for pixel at position vcell, ucell.  
   */
  int encodePixelID(int vcell, int ucell) override;
  
  /** Compute ucell and vcell from pixelID.
   */ 
  void decodePixelID(int& vcell, int& ucell, int uniqPixelID) override;

  /** Returns true if pixels at position (vcell1,ucell1) and (vcell2,ucell2) are neighbors.
   */ 
  bool areNeighbors(int vcell1, int ucell1, int vcell2, int ucell2) override;

  /** Get nominal sensor frame (i.e. where the detector is supposed to be)
   */
  void SetNominalFrame(const ReferenceFrame& nominal) { m_nominal = nominal; }
  
  /** Get nominal sensor frame (i.e. where the detector is supposed to be)
   */
  const ReferenceFrame & GetNominal() const override { return m_nominal; }  
	
  /** Get discrete rotation (mounting orientation of sensor) 
   */
  const ReferenceFrame & GetDiscrete() const override { return m_discrete; }
  
 private:

  /** Set the TH2Poly bins aka pixel, tuple< u,v,type,centreu,centrev >, tuple<type, centreu,centrev, points vector (coordu, coordv)
   */
  void SetCells(const std::vector< std::tuple<int, int, int, double, double> >& cells, const std::vector<std::tuple<int, double, double, std::vector<std::tuple<double, double>>>>& protopix);
 
  /** Gets the coordinates of the pixel vcell, ucell, helper function combining the get functions for the individual coordinates.
   */
  void PolyDet::GetPixelCenterCoord(double& vcoord, double& ucoord, int vcell, int ucell)
  /** Check if module is crossed (including supports)
   */
  bool ModuleCrossed(double u, double v);

  /** Get u type for pixel at position vcell and ucell. 
   */
  int GetPixelTypeU(int vcell, int ucell);   
   
  /** Get v type for pixel at position vcell and ucell. 
   */
  int GetPixelTypeV(int vcell, int ucell);   

  // TH2Poly object that defines the class, describing the pixel layout
  TH2Poly *m_layout;
  // Cells, tuple< u,v,type,centreu, centrev >
  std::vector< std::tuple< int, int, int, double, double > > m_cells;
  // Distance from cell to an other to be counted as neighbour odered by pixel type. tuple< type, distu, distv >
  std::vector< std::tuple< int, double, double> > m_cells_neighb_dist;
  // Pitch/bounding box of pixeltype
  std::vector< std::tuple< int, double, double > > m_pitch;
  // Thickness in sensitive volume
  double m_sensitiveThickness;
  // Rad. length in sensitive volume
  double m_sensitiveRadLength;
  // Atomic number in sensitive volume
  double m_sensitiveAtomicNumber;  
  // Atomic mass in sensitive volume
  double m_sensitiveAtomicMass; 
  // Size along u of sensitive  volume
  double m_sensitiveSizeU;  
  // Size along v of sensitive volume  
  double m_sensitiveSizeV; 
  // Size along u of sensitive + supports  
  double m_ladderSizeU; 
  // Size along v of sensitive + supports 
  double m_ladderSizeV; 
  // Thickness in nonsensitive volume
  double m_ladderThickness;
  // Rad. length in nonsensitive volume
  double m_ladderRadLength;
  // Atomic number in nonsensitive volume
  double m_ladderAtomicNumber;  
  // Atomic mass in nonsensitive volume
  double m_ladderAtomicMass; 
       
  // Nominal reference frame
  ReferenceFrame m_nominal;
  // Mounting rotations  
  ReferenceFrame m_discrete; 

  // Some helper variables
  int m_nCellsU; 
  int m_nCellsV; 
  int m_minCellU;
  int m_minCellV; 
  std::vector<double> m_offsetsU;
  std::vector<double> m_offsetsV; 
};
 
} // Namespace
//#include<Eigen/StdVector>
//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(depfet::PolyDet)
#endif // POLYDET_H
