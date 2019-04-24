#ifndef DET_H
#define DET_H 1

#include <string>
#include <vector>
#include <tuple>
#include <map>

// Include TBTools 
#include "ReferenceFrame.h"

namespace depfet {

/** Class Det
 *  
 *  Base class for pixel detector planes in tbsw.
 *  
 *  The Det class provides a user interface for detector information to 
 *  reconstruction algorithms and Marlin processors. In order to provide 
 *  support for new detector types, developers must override all  
 *  pure virtual members.   
 *    
 *  The Det class contains some objects with a special meaning in tbsw: 
 *     
 *  1) The SensorID index is unique for a pixel detector in a telescope. This number must be the same as the sensorID 
 *     stored in LCIO::TrackerRawData and LCIO::TrackerData collections.  
 * 
 *  2) The sensitive area is divided in non overlapping geometrical pixels addressed by an unique vCell/uCell index 
 *     pair. The indices must be integers in the bounding box [minU,maxU]x[minV,maxV]. 
 *     
 *  3) The vCell (uCell) ID increments in the along the v (u) axis of the local uvw sensor coordinate 
 *     system.  
 *  
 *  4) All pixels have an integer valued pixeltype. Pixels with identical type are copies of a protopixel
 *     placed a certain position on the sensor area, called the pixel center. 
 *  
 *  5) The outline (circumference) of a protopixel on the sensor are is approximated as a 2d polygon. In the 
 *     simplest case, this is just a square. 
 *   
 *  6) The nominal reference frame encodes the trafo from global xyz to local uvw coordinates. The nomal reference
 *     frame is first loaded from the gear file and, if available, updated from the alignment data base.
 *  
 *  Other global conventions in tbsw regarding the telescope geometry are detailed in TBDetector.h. 
 * 
 *  
 *  @Author B. Schwenker, University of Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
 
class Det {
   	
 public:
  
  /** Default constructor - subclasses need to call this in their
   *  default constructor.
   */
  Det(const std::string& typeName, int sensorID, int planeNumber) ; 
  
  /** Destructor */
  virtual ~Det() {}
  
  /** Get pixel type for pixel at position vcell and ucell. 
   *  Has to be implemented by subclasses.
   */
  virtual int GetPixelType(int vcell, int ucell) const = 0;  
  
  /** Get map of protopixels. The map keys are the pixeltypes and values are the vectors of polygon edges.  
   */
  virtual const std::map<int, std::vector<std::tuple<double,double>>> & GetProtopixels() const = 0;  
  
  /** Get the maximum uCell on the pixel matrix.  
   *  Has to be implemented by subclasses. 
   *  The uCell numbers of pixels are in the intervall [min,max].
   */  
  virtual int GetMaxUCell() const = 0;  

  /** Get the minimum uCell on the pixel matrix.  
   *  Has to be implemented by subclasses. 
   *  The uCell numbers of pixels are in the intervall [min,max].
   */  
  virtual int GetMinUCell() const = 0;  

  /** Get the maximum vCell on the pixel matrix.  
   *  Has to be implemented by subclasses. 
   *  The vCell numbers of pixels are in the intervall [min,max].
   */  
  virtual int GetMaxVCell() const = 0;  
  
  /** Get the minimum vCell on the pixel matrix.  
   *  Has to be implemented by subclasses. 
   *  The vCell numbers of pixels are in the intervall [min,max].
   */  
  virtual int GetMinVCell() const = 0;  
  
  /** Get the maximun u position of sensitive area. 
   *  Has to be implemented by subclasses. 
   *  All pixels must be inside the box 
   *  [-SensMinU,SensMaxU]x[-SensMinV,SensMaxV].
   */  
  virtual double GetSensitiveMaxU() const = 0;  
  
  /** Get minimum u position of sensitive area. 
   *  Has to be implemented by subclasses. 
   *  All pixels must be inside the box 
   *  [-SensMinU,SensMaxU]x[-SensMinV,SensMaxV].
   */  
  virtual double GetSensitiveMinU() const = 0;  
  
  /** Get maximum v position of sensitive area. 
   *  Has to be implemented by subclasses. 
   *  All pixels must be inside the box 
   *  [-SensMinU,SensMaxU]x[-SensMinV,SensMaxV].
   */  
  virtual double GetSensitiveMaxV() const = 0; 
  
  /** Get minimum v position of sensitive area. 
   *  Has to be implemented by subclasses. 
   *  All pixels must be inside the box 
   *  [-SensMinU,SensMaxU]x[-SensMinV,SensMaxV].
   */  
  virtual double GetSensitiveMinV() const = 0; 
      
  /** Get u pitch for pixel at position vcell,ucell. 
   *  Has to be implemented by subclasses.  
   */ 
  virtual double GetPitchU(int vcell, int ucell) const = 0; 
  
  /** Get v pitch for pixel at position vcell,ucell. 
   *  Has to be implemented by subclasses.  
   */ 
  virtual double GetPitchV(int vcell, int ucell) const = 0;  
   
  /** Returns true if point (u,v,w) is not inside the sensitive volume. 
   *  Has to be implemented by subclasses.  
   */ 
  virtual bool isPointOutOfSensor( double u , double v , double w = 0) const = 0; 
  
  /** Returns true if point (u,v,w) is inside the sensitive volume. 
   *  Has to be implemented by subclasses.  
   */ 
  virtual bool SensitiveCrossed(double u, double v, double w = 0) const = 0;
      
  /** Get thickness of detector at position (u,v).
   *  Has to be implemented by subclasses.  
   */ 
  virtual double GetThickness(double u, double v) const = 0;

  /** Get length of particle track intersecting the detector.
   *  Has to be implemented by subclasses.  
   */
  virtual double GetTrackLength(double u, double v, double dudw, double dvdw) const = 0;  
    
  /** Get radlenght at position (u,v).
   *  Has to be implemented by subclasses.  
   */ 
  virtual double GetRadLength(double u, double v) const = 0;

  /** Get atomic number at position (u,v).
   *  Has to be implemented by subclasses.  
   */
  virtual double GetAtomicNumber(double u, double v) const = 0; 
 
  /** Get atomic mass at position (u,v).
   *  Has to be implemented by subclasses.  
   */
  virtual double GetAtomicMass(double u, double v) const = 0;   
  
  /** Returns uCell of pixel at position (u,v). Returns -1 if there is no pixel.  
   *  Has to be implemented by subclasses.  
   */    
  virtual int GetUCellFromCoord( double u, double v ) const = 0;
  
  /** Returns vCell of pixel at position (u,v). Returns -1 if there is no pixel.  
   *  Has to be implemented by subclasses.  
   */    
  virtual int GetVCellFromCoord( double u, double v ) const = 0; 
  
  /** Returns u position of pixel center. 
   *  Has to be implemented by subclasses.  
   */
  virtual double GetPixelCenterCoordU(int vcell, int ucell) const = 0;
  
  /** Returns v position of pixel center. 
   *  Has to be implemented by subclasses.  
   */
  virtual double GetPixelCenterCoordV(int vcell, int ucell) const = 0;
  
  /** Return unique ID for pixel at position vcell, ucell.  
   *  Has to be implemented by subclasses.  
   */
  virtual int encodePixelID(int vcell, int ucell) const = 0;
  
  /** Compute ucell and vcell from pixelID.
   *  Has to be implemented by subclasses.  
   */ 
  virtual void decodePixelID(int& vcell, int& ucell, int uniqPixelID) const = 0;

  /** Returns true if pixels at position (vcell1,ucell1) and (vcell2,ucell2) are neighbors.
   *  Has to be implemented by subclasses.  
   */ 
  virtual bool areNeighbors(int vcell1, int ucell1, int vcell2, int ucell2) const = 0;
  
  /** Get nominal sensor frame (i.e. where the detector is supposed to be)
   */
  virtual const ReferenceFrame & GetNominal() const = 0;    

  /** Set nominal sensor frame. This is needed for applying alignment corrections. 
   */
  virtual void SetNominalFrame(const ReferenceFrame& nominal) = 0;  
	
  /** Get discrete rotation (mounting orientation of sensor) 
   */
  virtual const ReferenceFrame & GetDiscrete() const = 0;
  
  /** Print method.   
   */ 
  virtual void Print() const;
  
  /** Get SensorID 
  */
  int GetSensorID() const { return m_sensorID; }
  
  /** Get plane number
   */   
  int GetPlaneNumber() const { return m_planeNumber; }

  /** Set plane number
   */   
  void SetPlaneNumber(int planeNumber)  { m_planeNumber=planeNumber; }
  
  /** Get type name (as set in constructor)
   */
  const std::string & GetType() const { return m_typeName ; } 
       
 private:
  
  // Type of sensor 
  std::string m_typeName; 
  // SensorID
  int m_sensorID;
  // Plane number along beam line 
  int m_planeNumber; 
  
};
 
} // Namespace

#endif // DET_H
