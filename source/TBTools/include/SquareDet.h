#ifndef SQUAREDET_H
#define SQUAREDET_H 1

// Include basic C
#include <vector>
#include <tuple>
#include <map>

// Include TBTools 
#include "Det.h"

// ROOT includes 
#include <TH2D.h>

namespace depfet {
	
/** Class SquareDet
 *    
 *  The SquareDet class is a subclass of the Det class and represents a pixel detector 
 *  with a checkerboard type pixel matrix. The pixel pitch may change along the local 
 *  u-axis or v-axis. See also Det.h and TBDetector.h files for more information. 
 *  
 *  The geometrical 2D layout of pixels is constructed from a product of two 1D 'strip' layouts 
 *  for the u and v axis. As the 'strip' pitch in the 1D layouts changes, also the area of the 
 *  2D pixel cells will change as well. Neighboring 'strips' with the same pitch are called a 
 *  CellGroup.     
 *  
 *  The definition of the 1D layouts is read from the gear file. For example, a Belle II pixel 
 *  sensor has 768 (vCells) and 250 columns (uCells). The row pitch is large for rows [0,511] 
 *  and fine for rows [512,767]. The column pitch is constant. 
 *  
 *  <pre>
 *  <vCellGroup
 *    minCell="0"
 *    maxCell="511"
 *    pitch="0.060		 		
 *  /vCellGroup>
 *  <vCellGroup
 *    minCell="512"
 *    maxCell="767"
 *    pitch="0.055"
 *  /vCellGroup>
 *  <uCellGroup
 *    minCell="0"
 *    maxCell="249"
 *    pitch="0.050"
 *  /uCellGroup>
 *  </pre>
 *
 *  Note that there is no restriction on the number of vCellGroups or uCellGroups and the above 
 *  example just demonstrates the simplest non trivial case. 
 *
 *  The material budget, or radiation length X0, for a SquareDet can optionally be specified by 
 *  providing a (TH2D) histogram with X0 numbers for bins in the local sensor u/v plane. The root 
 *  file containing the histogram must be placed in the steerfiles folder next to the gear files. 
 *  Outside the range of this 2d histogram, default X0 numbers defined in the xml nodes 'sensitive'
 *  and 'ladder' will be used. 
 *  
 *  In order to use this feature, add the following xml node to your sensor layer in the gear file. 
 *  
 *  <pre>
 *  <x0map 
 *    x0file="EoS2_f.root"
 *	  x0object="x0_image"
 *	/>
 *  </pre>
 *  
 *  @Author B. Schwenker, University of Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
 
class SquareDet : public Det {
   	
 public:
   
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
                     const ReferenceFrame& discrete, const ReferenceFrame& nominal,
					 const std::string& x0FileName, const std::string& x0ObjName );
  
    
  /** Get pixel type for pixel at position vcell and ucell. 
   */
  int GetPixelType(int vcell, int ucell) const override;  

  /** Get map of protopixels. The map keys are the pixeltypes and values are the vectors of polygon edges.  
   */
  const std::map<int, std::vector<std::tuple<double,double>>> & GetProtopixels() const override {return m_protopixels;}  

  /** Get the maximum uCell on the pixel matrix.   
   *  The uCell numbers of pixels are in the intervall [min,max].
   */  
  int GetMaxUCell() const override;  

  /** Get the minimum uCell on the pixel matrix.   
   *  The uCell numbers of pixels are in the intervall [min,max].
   */  
  int GetMinUCell() const override;  

  /** Get the maximum vCell on the pixel matrix.   
   *  The vCell numbers of pixels are in the intervall [min,max].
   */  
  int GetMaxVCell() const override;  

  /** Get the minimum vCell on the pixel matrix.   
   *  The vCell numbers of pixels are in the intervall [min,max].
   */  
  int GetMinVCell() const override;  
  
  /** Get maximum u position of sensitive area. 
   */  
  double GetSensitiveMaxU() const override;  

  /** Get minimum u position of sensitive area. 
   */  
  double GetSensitiveMinU() const override;  
  
  /** Get maximum v position of sensitive area. 
   */  
  double GetSensitiveMaxV() const override;  

  /** Get minimum v position of sensitive area. 
   */  
  double GetSensitiveMinV() const override;  
      
  /** Get u pitch for pixel at position vcell,ucell.
   *  For a SquareDet, the u pitch depends only on the ucell.   
   */ 
  double GetPitchU(int /*vcell*/, int ucell) const override; 
  
  /** Get v pitch for pixel at position vcell,ucell.   
   *  For a SquareDet, the v pitch depends only on the vcell.   
   */ 
  double GetPitchV(int vcell, int /*ucell*/) const override;  
   
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
   *  For a SquareDet, the ucell depends only on the u position.      
   */    
  int GetUCellFromCoord( double u, double /*v*/ ) const override;
  
  /** Returns vCell of pixel at position (u,v). Returns -1 if there is no pixel.  
   *  For a SquareDet, the vcell depends only on the v position.   
   */    
  int GetVCellFromCoord( double /*u*/, double v ) const override; 
  
  /** Returns u position of pixel center. 
   *  For a SquareDet, the u position depends only on the ucell.   
   */
  double GetPixelCenterCoordU(int /*vcell*/, int ucell) const override;
  
  /** Returns v position of pixel center.  
   *  For a SquareDet, the v position depends only on the vcell.    
   */
  double GetPixelCenterCoordV(int vcell, int /*ucell*/) const override;
  
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

  /** Set the u cells 
   */ 
  void SetCellsU( std::vector< std::tuple<int,int,double> > uCells); 
  
  /** Set the v cells 
   */ 
  void SetCellsV( std::vector< std::tuple<int,int,double> > vCells); 
 
  /** Check if module is crossed (including supports)
   */
  bool ModuleCrossed(double u, double v) const ;

  /** Get u type for pixel at position ucell. 
   */
  int GetPixelTypeU(int ucell) const;   
   
  /** Get v type for pixel at position vcell. 
   */
  int GetPixelTypeV(int vcell) const;   

  /** Get default radlenght at position (u,v) from sensitive / ladder specifications.
   */ 
  double GetRadLengthDefault(double u, double v) const;  

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
  // Map containing the protopixels of the layout
  std::map<int, std::vector<std::tuple<double,double>>> m_protopixels;

  // For using user specified X0 maps 
  TH2D * m_x0map {nullptr};
};
 
} // Namespace
#endif // SQUAREDET_H
