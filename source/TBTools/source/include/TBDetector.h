#ifndef TBDETECTOR_H
#define TBDETECTOR_H 1

// Local includes
#include "Det.h"

// Include basic C
#include <vector>
#include <map>
#include <math.h>

// Include CLHEP classes
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/Vector.h>

namespace depfet {

  //! Class TBDetector   
  /*! 
   *  The TBDetecor class holds all information about the geometry of a pixel tracking 
   *  telescope. This includes detector layout data, detector alignment data and 
   *  material budget data. The class defines a common interface to be used in all 
   *  reconstruction modules for analysis of test beam data.  
   *  
   *  GLOBAL COORDINATE FRAME 
   *  
   *  The global frame of reference XYZ is a right-handed cartesian coordinate system. 
   *  It is choosen such that the  global Z axis points along the beam axis. 
   *  The Y axis is vertical going from the bottom to the top. The X axis is horizontal, 
   *  choosen to complete a right handed coordinate system. 
   *   
   *  LOCAL COORDINATE FRAME
   *  
   *  A local frame of reference UVW is attached to each detector in the 
   *  telescope setup. The origin of the local frame is the center of the active 
   *  volume, at half-width, half-lenght and half-depth of the detector. The local
   *  u and v axes are parallel to the sensitive detector plane and the w- axis is perp.
   *  to the plane. Here, the direction of the  u- axis is defined by increasing 
   *  pixel columns and the v-axis points in the direction of increasing rows. 
   *  It means that the pixel at (column=0,row=0) is in the lower left corner of 
   *  the local coordinate frame.   
   *  
   *  The transformation law from local q=(u,v,w) to global r=(x,y,z) coordinates 
   *  is in matrix notation
   *  
   *  q = R0 * (r-r0) , det(R0) = 1
   *  
   *  The translation vector r0=(y0,y0,z0) points from the global origin to the local 
   *  sensor origin. The 3x3 rotation matrix is decomposed into a product
   *  
   *  R0 = R*D ,   det(D) = 1 and det(R) = 1
   * 
   *  of a discrete rotation matrix D and a continous rotation matrix R. The discrete 
   *  matrix D takes into account all big rotations defined by the sensor mounting. It 
   *  has four discete rotation parameters and is given by  
   *  
   *  D(1,1) =  rotation1  ,  D(1,2) =  rotation2  , D(1,3) =  0
   *  D(2,1) =  rotation3  ,  D(2,2) =  rotation4  , D(2,3) =  0
   *  D(3,1) =  0          ,  D(3,2) =  0          , D(3,3) =  sign
   *  
   *  The parameters are either +/-1 or 0. The component D(3,3) is determined by 
   *  det(D)=1 condition. The matrix D reflects the different possibilities to 
   *  install a detector in the beam.
   *  
   *  Finally, the continous rotation matrix R is parametrized by three Euler angles
   *  around local x-, y-, and z- axis.
   *  
   *  R = R_g(gamma) * R_b(beta) *  R_a(alpha).
   *  
   *  For the telescope sensors, apha, beta and gamma should always be small (~1mrad). The 
   *  X and Y components of r0 should also be small for a good telescope (~1mm). For 
   *  the DUT, alpha and/or beta may be large in so called tilt angle scans using rotation
   *  stages.  
   *   
   *  PIXEL DATA
   *  
   *  The pixel data coming from the DAQ is addressed by an integer column and row 
   *  number counting from zero. The measured raw signals are integer valued. 
   *  Calibrated signal values may be doubles or integers. 
   *  
   *  The position of the center of pixel p=(col,row,0) in local coordinates depends on 
   *  the physical arrangement of pixel. In case the pixels form a rectangular grid with 
   *  pixel size pitch_u and pitch_v it is  
   *  
   *  q = PITCH * (p - p0 ) ; p = (col, row, 0) 
   *   
   *  where PITCH=diag(pitch_u,pitch_v,0) and shift vector p0 = (nxpixels/2,nypixels/2,0)  
   *  pointing from lower left edge (0,0) to the center of the active sensor area.  
   *  
   *  COORDINATE TRANSFORMATIONS 
   *  
   *  The concept of a reference frame (RotationMatrix + TranslationVector) is 
   *  used to do transformations between coordinate systems. The so called NOMINAL 
   *  reference frame transforms from 'local' to 'global' frame. 
   *  
   *  SENSOR PLANE NUMBERING 
   *  
   *  Each sensor has a unique DAQ ID which is always attached to the readout 
   *  stream from global DAQ system. In addition, there is the geometrical order
   *  of the detectors along the beam axis called plane number. The ordering of 
   *  detector is determined from the sensor's z0 positions in the gear file. The 
   *  detector with smallest z0 gets plane number zero. 
   *     
   *  PHYSICAL UNITS
   *  
   *  The unit system is mm for length and rad for angles. The unit for Particle
   *  energy and momentum is GeV. The unit of mass is g. The radiation length is 
   *  given in units of mm. Angles are given in units of radians.
   *  
   *  @Author B. Schwenker, University of Göttingen (October 2010)
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
  
class TBDetector {
  
 public:
  
  //!Constructor 
  TBDetector( ); 
  
  //!Destructor
  ~TBDetector( );
  
  //! Write alignment data base file - overwrites old DB file
  void WriteAlignmentDB( );
    
  //! Build detector from gear file 
  void ReadGearConfiguration( );
  
  //! Read alignment data base file 
  void ReadAlignmentDB( std::string FileName );

  //! Read alignment data base file name 
  void SetAlignmentDBName( std::string FileName );
  
  //! Get number of pixel sensors 
  int GetNSensors() { return _numberOfSensors; };
  
  //! Translate DAQ (Sensor) ID to plane number  
  int GetPlaneNumber(int sensorID);
  
  //! Get detector at plane number ipl
  Det& GetDet(int ipl);

  //! Get components of magnetic field in Tesla
  double GetBx() { return _Bx;};
  double GetBy() { return _By;};
  double GetBz() { return _Bz;};
  
  //! Method printing general Gear parameters
  void Print() ;
     
 private:
       
  // Cartesian components of magnetic field
  double _Bx; 
  double _By; 
  double _Bz; 

  // Name of LCIO data base file 
  std::string _alignmentDBFileName;
  
  // Total number of sensors  
  int _numberOfSensors;       
  
  // Maps sensorID to planeNumber 
  std::map< int, int > _indexMap;   
  
  // Vector of pixel detector in beam 
  std::vector<Det> _DetVec;      
  
}; // End class TBDetector
 
} // Namespace

#endif // TBDETECTOR_H
