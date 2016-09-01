#ifndef SIPXLGEOM_H
#define SIPXLGEOM_H 1

// Include basic C
#include <vector>

// Include CLHEP classes
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Vector/ThreeVector.h>

// Include Gear header files
#include <gear/GEAR.h>
#include <gear/GearParameters.h>

namespace sipxl {

// Define constants
#define EPS 5           // um
#define LAYERCOD   1000 // Const to encode and decode layers
#define LADDERCOD    10 // Const to encode and decode ladders
#define SENSORCOD     1 // Const to encode and decode sensors

// Define pixel type
enum PixelType { RPhi = 0, Z = 1 };

// Define layer type
enum LayerType { pixel = 0, stripB = 1, stripF = 2};

//! Gear geometry class - holds all geometry information about silicon pixel
//! sensors. The data are taken directly from Gear xml file and values are
//! saved in the system of units defined in PhysicalConstants.h. The local
//! reference system for VXD is defined as follows: X-axis is perpendicular to the
//! beam axis and to the ladder (sensor) plane; Y-axis lies in a ladder
//! (sensor) plane and is also perpendicular to the beam axis; Z-axis is
//! parallel to the beam axis (for zero theta); for TBSiDet it's defined such as
//! X-axis represents the beam axis and Y-axis (RPhi) and Z-axis are parallel,
//! the reason is to have the same geometry picture as for VXD; (0,0,0) point
//! is positioned for both geometries such as X, Y, Z coordinates are always
//! positive.
//!
//! @author Z. Drasal, Charles University Prague
//!

class SiPxlGeom {

 public:

//!Constructor
   SiPxlGeom()  : _gearType(""), _gearParamsGen(0), _gearParams(0), _numberOfLayers(0) {;}

//!Destructor
	~SiPxlGeom();

//!Method initializing class - reads Gear parameters from XML file
	void initGearParams();


// MAGNETIC FIELD

//!Get magnetic field - x
	double getBx() { return _magField.getX(); }

//!Get magnetic field - y
	double getBy() { return _magField.getY(); }

//!Get magnetic field - z
	double getBz() { return _magField.getZ(); }

//!Get magnetic field - Three vector
	CLHEP::Hep3Vector get3MagField() { return _magField; }


// GEOMETRY PROPERTIES

//!Get gear type
   std::string getGearType() { return _gearType; }

// ENCODING

//!Encode sensorID
   int encodeSensorID(short int layerID, short int ladderID, short int sensorID) const;

//!Decode sensorID
   void decodeSensorID(short int & layerID, short int & ladderID, short int & sensorID, int uniqSensorID) const;

//!Encode pixelID
   int encodePixelID(short int layerID, short int sensorID, int pixelIDRPhi, int pixelIDZ) const;

//!Decode pixelID
   void decodePixelID(int & pixelIDRPhi, int & pixelIDZ, short int layerID, short int sensorID, int uniqPixelID) const;

// LAYER PROPERTIES

//!Get number of layers
   short int getNLayers() const {return _numberOfLayers;}

//!Get number of PXD layers
   short int getNPXDLayers() const;

//!Get layer real ID
   int getLayerRealID(short int layerID) const;

//!Transform real layer ID to C-type numbering 0 - n ...
   short int getLayerIDCTypeNo(int realLayerID) const;

//!Get layer type
   short int getLayerType(short int layerID) const;

// LADDER PROPERTIES

//!Get number of ladders
	short int getNLadders(short int layerID) const;

// SENSOR PROPERTIES

//!Get number of sensors for given ladder
	short int getNSensors(short int layerID) const;

//!Get number of pixels in Z axis (in each sensor)
	int getSensorNPixelsInZ(short int layerID, short int sensorID) const;

//!Get number of pixels in R-Phi (in each sensor)
	int getSensorNPixelsInRPhi(short int layerID, short int sensorID) const;

//!Get sensor pitch in Z axis
	double getSensorPitchInZ(short int layerID, short int sensorID) const;

//!Get sensor pitch in R-Phi
	double getSensorPitchInRPhi(short int layerID, short int sensorID) const;

//!Get sensor thickness
	double getSensorThick(short int layerID, short int sensorID) const;

//!Get sensor width
	double getSensorWidth(short int layerID, short int sensorID) const;

//!Get sensor length
	double getSensorLength(short int layerID, short int sensorID) const;


// TRANSFORMATION METHODS - GLOBAL REF. SYSTEM

//!Transform given point from global ref. system to local ref. system (sensor)
	CLHEP::Hep3Vector transformPointToLocal(short int layerID, short int ladderID, short int sensorID, const CLHEP::Hep3Vector & point);

//!Transform given vector from global ref. system to local ref. system (sensor)
	CLHEP::Hep3Vector transformVecToLocal(short int layerID, short int ladderID, short int sensorID, const CLHEP::Hep3Vector & vec);

//!Transform given matrix from global ref. system to local ref. system (sensor)
   CLHEP::HepMatrix transformMatxToLocal(short int layerID, short int ladderID, short int sensorID, const CLHEP::HepMatrix & matrix);


// TRANSFORMATION METHODS - LOCAL REF. SYSTEM

//!Transform given point from local ref. system (sensor) to global ref. system
	CLHEP::Hep3Vector transformPointToGlobal(short int layerID, short int ladderID, short int sensorID, const CLHEP::Hep3Vector & point);

//!Transform given vector from local ref. system (sensor) to global ref. system
	CLHEP::Hep3Vector transformVecToGlobal(short int layerID, short int ladderID, short int sensorID, const CLHEP::Hep3Vector & vec);

//!Transform given diagonal matrix from local ref. system (sensor) to global ref. system
   CLHEP::HepMatrix transformMatxToGlobal(short int layerID, short int ladderID, short int sensorID, const CLHEP::HepMatrix & matrix);


// OTHER METHODS - LOCAL REF. SYSTEM

//!Get info whether the given point is out of Si sensor
	bool isPointOutOfSensor(short int layerID, short int sensorID, const CLHEP::Hep3Vector & point) const;

//!Correct given point if out of sensor
	void correctPointOutOfSensor(short int layerID, short int sensorID, CLHEP::Hep3Vector & point) const;

//!Get pixel position in Z direction
	double getPixelPosInZ(short int layerID, short int sensorID, int pixelIDZ) const;

//!Get pixel position in R-Phi direction
   double getPixelPosInRPhi(short int layerID, short int sensorID, bool bricked, int pixelIDRPhi, int pixelIDZ) const;

//!Get pixel ID in Z direction
   int getPixelIDInZ(short int layerID, short int sensorID, double posZ) const;

//!Get pixel ID in R-Phi direction
   int getPixelIDInRPhi(short int layerID, short int sensorID, bool bricked, double posRPhi, double posZ) const;

// PRINT METHODS

//!Method printing general Gear parameters
	void printGearParams() const;

//!Method printing sensor Gear parameters
	void printSensorParams(short int layerID) const;

 protected:

// LAYER PROPERTIES - PROTECTED

//!Get layer radius
   double getLayerRadius(short int layerID) const;

//!Get layer phi zero angle
   double getLayerPhi0(short int layerID) const;

// LADDER PROPERTIES - PROTECTED

//!Get ladder thickness
   double getLadderThick(short int layerID) const;

//!Get ladder width
   double getLadderWidth(short int layerID) const;

//!Get ladder length
   double getLadderLength(short int layerID) const;

//!Get ladder offset in Y
   double getLadderOffsetY(short int layerID) const;

//!Get ladder offset in Z
   double getLadderOffsetZ(short int layerID) const;

//!Get ladder rotation - phi angle
   double getLadderPhi(short int layerID, short int ladderID) const;

//!Get ladder rotation - theta angle
   double getLadderTheta(short int layerID) const;

// SENSOR PROPERTIES - PROTECTED

//!Get sensor pos in X
   double getSensorPosX(short int sensorID) const;

//!Get sensor pos in Y
   double getSensorPosY(short int sensorID) const;

//!Get sensor pos in Y
   double getSensorPosZ(short int sensorID) const;

//!Get sensor rotation around X axis
   double getSensorRotX(short int sensorID) const;

//!Get gap size inbetween sensors
   double getSensorGapInBetween(short int layerID) const;

//!Get width of sensor rim in Z (passive part of silicon)
   double getSensorRimWidthInZ(short int layerID) const;

//!Get width of sensor rim in R-Phi (passive part of silicon)
   double getSensorRimWidthInRPhi(short int layerID) const;

 private:

 	std::string _gearType;                       //!< GearType
   const gear::GearParameters * _gearParamsGen; //!< General gear parameters
   const gear::GearParameters * _gearParams;    //!< Concrete gear parameters

   // Magnetic field
   CLHEP::Hep3Vector _magField;

   // Geometry parameters - layers
   short int _numberOfLayers;

   std::vector<int> _layerRealID;
   std::vector<int> _layerType;

   std::vector<double> _layerRadius;
   std::vector<double> _layerPhi0;
   std::vector<double> _layerTheta;

   // Geometry parameters - ladders
   std::vector<int> _numberOfLadders;

   std::vector<double> _ladderThick;
   std::vector<double> _ladderWidth;
   std::vector<double> _ladderLength;
   std::vector<double> _ladderOffsetY;
   std::vector<double> _ladderOffsetZ;

   // Geometry - sensors
   std::vector<int> _numberOfSensors;

   std::vector<int> _sensorNPixelsInZ;
   std::vector<int> _sensorNPixelsInRPhi;

   std::vector<double> _sensorPitchInZ;
   std::vector<double> _sensorPitchInRPhi;
   std::vector<double> _sensorThick;
   std::vector<double> _sensorWidth;
   std::vector<double> _sensorLength;
   std::vector<double> _sensorPosX;
   std::vector<double> _sensorPosY;
   std::vector<double> _sensorPosZ;
   std::vector<double> _sensorRotX;
   std::vector<double> _sensorGapInBtw;
   std::vector<double> _sensorRimWidthInZ;
   std::vector<double> _sensorRimWidthInRPhi;

}; // Class

} // Namespace

#endif // SIPXLGEOM_H
