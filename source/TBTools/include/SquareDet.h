#ifndef SQUAREDET_H
#define SQUAREDET_H 1

// Include basic C
#include <vector>
#include <tuple>

// Include TBTools 
#include "Det.h"

namespace depfet {
	
/** Class SquareDet
 *    
 *  The SquareDet class is a subclass of the Det class and represents a pixel detector 
 *  with a checkerboard type pixel matrix. The pixel pitch may change along the local 
 *  u-axis or v-axis.  
 *  
 *  The pixel matrix has pixel at positions vCell, uCell in the range [0,maxUCell] and 
 *  [0,maxVCell]. The center of the pixel is the geometrical center of the square area 
 *  with height GetPitchV(vCell,uCell) and width GetPitch(vCell,uCell) in local sensor
 *  coordinates. The total sensitve area of the entire pixel matrix has height 
 *  GetSensitiveSizeV() and width GetSensitiveSizeU(). 
 *  
 *  @Author B. Schwenker, University of Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
 
class SquareDet : public Det {
   	
 public:
//  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
 
  
      
  /** Constructor
   */
  SquareDet(const std::string& typeName, int sensorID, int planeNumber) ; 
  
  /** Constructor
   */
  SquareDet(const std::string& typeName, int sensorID, int planeNumber, 
                     double sensThick, double sensRadLenght, double sensAtomicNumber,
                     double sensAtomicMass, double ladderThick, double ladderRadLength, 
                     double ladderAtomicNumber, double ladderAtomicMass, double ladderSizeU, 
                     double LayerSizeV, const std::vector< std::tuple<int,int,double> >& uCells, 
                     const std::vector< std::tuple<int,int,double> >& vCells, 
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

  /** Set nominal sensor frame. This is needed for applying alignment corrections. 
   */
  void SetNominalFrame(const ReferenceFrame& nominal) { m_nominal  = nominal; }   
  
  /** Get nominal sensor frame (i.e. where the detector is supposed to be)
   */
  const ReferenceFrame & GetNominal() const override { return m_nominal; }  
	
  /** Get discrete rotation (mounting orientation of sensor) 
   */
  const ReferenceFrame & GetDiscrete() const override { return m_discrete; }
  
 private:

  /** Set the u cells 
   */ 
  void SetCellsU( std::vector< std::tuple<int,int,double> > uCells); 
  
  /** Set the v cells 
   */ 
  void SetCellsV( std::vector< std::tuple<int,int,double> > vCells); 
 
  /** Check if module is crossed (including supports)
   */
  bool ModuleCrossed(double u, double v);

  /** Get u type for pixel at position vcell and ucell. 
   */
  int GetPixelTypeU(int vcell, int ucell);   
   
  /** Get v type for pixel at position vcell and ucell. 
   */
  int GetPixelTypeV(int vcell, int ucell);   

  
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
  // Cells along sensor u axis  
  std::vector< std::tuple<int,int,double> > m_uCells;
  // Cells along sensor v axis   
  std::vector< std::tuple<int,int,double> > m_vCells ; 
  std::vector<double> m_offsetsU;
  std::vector<double> m_offsetsV; 
};
 
} // Namespace
//#include<Eigen/StdVector>
//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(depfet::SquareDet)
#endif // SQUAREDET_H
