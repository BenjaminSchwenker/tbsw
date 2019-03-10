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
 * 
 *  The PolyDet class is a subclass of the Det class and represents a pixel detector 
 *  with a pixel matrix of polygon shaped pixels.
 *  
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
  PolyDet(const std::string& typeName, int sensorID, int planeNumber, 
                 double sensThick, double sensRadLenght, double sensAtomicNumber,
                 double sensAtomicMass, double ladderThick, double ladderRadLength, 
                 double ladderAtomicNumber, double ladderAtomicMass, double ladderSizeU, 
                 double ladderSizeV, const std::vector< std::tuple<int,int,int,double,double> >& cells, 
		         const std::vector< std::tuple<int,double,double,std::vector<std::tuple<double,double>>>> & protocells,
                 const ReferenceFrame& discrete, const ReferenceFrame& nominal );

  /** Get pixel type for pixel at position vcell and ucell. 
   */
  int GetPixelType(int vcell, int ucell) const override;  

  /** Get the maximum uCell on the pixel matrix.   
   */  
  int GetMaxUCell() const override;  

  /** Get the maximum vCell on the pixel matrix.   
   */  
  int GetMaxVCell() const override;  

  /** Get the minimum uCell on the pixel matrix.   
   */  
  int GetMinUCell() const override;  

  /** Get the minimum vCell on the pixel matrix.   
   */  
  int GetMinVCell() const override;  
  
  /** Get maximum u position of sensitive area. 
   */  
  double GetSensitiveMaxU() const override;  
  
  /** Get maximum v position of sensitive area. 
   */  
  double GetSensitiveMaxV() const override;  

  /** Get minimum u position of sensitive area. 
   */  
  double GetSensitiveMinU() const override;  
  
  /** Get minimum v position of sensitive area. 
   */  
  double GetSensitiveMinV() const override;  
      
  /** Get u pitch for pixel at position vcell,ucell.   
   */ 
  double GetPitchU(int vcell, int ucell) const override; 
  
  /** Get v pitch for pixel at position vcell,ucell.   
   */ 
  double GetPitchV(int vcell, int ucell) const override;  
   
  /** Returns true if point (u,v,w) is not inside the sensitive volume.  
   */ 
  bool isPointOutOfSensor( double u , double v , double w = 0) const override; 
  
  /** Returns true if point (u,v,w) is inside the sensitive volume. 
   */ 
  bool SensitiveCrossed(double u, double v, double w = 0) const override;
      
  /** Get thickness of detector at position (u,v). 
   */ 
  double GetThickness(double u, double v) const override;

  /** Get length of particle track intersecting the detector.  
   */
  double GetTrackLength(double u, double v, double dudw, double dvdw) const override;  
    
  /** Get radlenght at position (u,v).
   */ 
  double GetRadLength(double u, double v) const override;

  /** Get atomic number at position (u,v).
   */
  double GetAtomicNumber(double u, double v) const override; 
 
  /** Get atomic mass at position (u,v).
   */
  double GetAtomicMass(double u, double v) const override;   
  
  /** Returns uCell of pixel at position (u,v). Returns -1 if there is no pixel.    
   */    
  int GetUCellFromCoord( double u, double v ) const override;
  
  /** Returns vCell of pixel at position (u,v). Returns -1 if there is no pixel.  
   */    
  int GetVCellFromCoord( double u, double v ) const override; 
  
  /** Returns u position of pixel center. 
   */
  double GetPixelCenterCoordU(int vcell, int ucell) const override;
  
  /** Returns v position of pixel center.   
   */
  double GetPixelCenterCoordV(int vcell, int ucell) const override;
  
  /** Return unique ID for pixel at position vcell, ucell.  
   */
  int encodePixelID(int vcell, int ucell) const override;
  
  /** Compute ucell and vcell from pixelID.
   */ 
  void decodePixelID(int& vcell, int& ucell, int uniqPixelID) const override;

  /** Returns true if pixels at position (vcell1,ucell1) and (vcell2,ucell2) are neighbors.
   */ 
  bool areNeighbors(int vcell1, int ucell1, int vcell2, int ucell2) const override;

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

  /** Set the TH2Poly bins aka pixel, tuple< u,v,type,centreu,centrev >, tuple<type, centreu,centrev, points vector (coordu, coordv)
   */
  void SetCells(const std::vector< std::tuple<int, int, int, double, double> >& cells, const std::vector<std::tuple<int, double, double, std::vector<std::tuple<double, double>>>>& protopix);
 
  /** Gets the coordinates of the pixel vcell, ucell, helper function combining the get functions for the individual coordinates.
   */
  void GetPixelCenterCoord(double& vcoord, double& ucoord, int vcell, int ucell) const;
  
  /** Check if module is crossed (including supports)
   */
  bool ModuleCrossed(double u, double v) const;

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
  int m_maxCellU;
  int m_maxCellV; 
  std::vector<double> m_offsetsU;
  std::vector<double> m_offsetsV; 
};
 
} // Namespace
//#include<Eigen/StdVector>
//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(depfet::PolyDet)
#endif // POLYDET_H
