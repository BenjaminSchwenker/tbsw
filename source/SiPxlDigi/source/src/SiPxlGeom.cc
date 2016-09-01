#include "SiPxlGeom.h"
#include "Colours.h"
#include "PhysicalConstants.h"

// Include basic C header files
#include <cstdlib>
#include <iomanip>
#include <algorithm>

// Include CLHEP header files
#include <CLHEP/Vector/EulerAngles.h>
#include <CLHEP/Vector/Rotation.h>

// Include Gear header files
#include <gear/BField.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>
#include <gear/TBSiParameters.h>
#include <gearimpl/Vector3D.h>

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

// Namespaces
using namespace CLHEP;
using namespace marlin;

namespace sipxl {

//
// Destructor
//
SiPxlGeom::~SiPxlGeom()
{
//	std::cout << "Deleting SiPxlGeom" << std::endl;
}

//
// Method initializing this class - reads Gear parameters from XML file
//
void SiPxlGeom::initGearParams()
{
	// BField
   try {

      const gear::BField & bField = Global::GEAR->getBField();

      _magField.setX( (bField.at( gear::Vector3D(0.,0.,0.) )).x() * T);
      _magField.setY( (bField.at( gear::Vector3D(0.,0.,0.) )).y() * T);
      _magField.setZ( (bField.at( gear::Vector3D(0.,0.,0.) )).z() * T);
   }
   catch (gear::UnknownParameterException& e) {}

   // VXDParameters
   try {

      const gear::GearParameters & paramsGen = Global::GEAR->getGearParameters("VXDParameters");
      const gear::GearParameters & params    = Global::GEAR->getVXDParameters();

      _gearParamsGen = &paramsGen;
      _gearParams    = &params;
      _gearType      = "VXD";

   }
   catch (gear::UnknownParameterException& e) {}

   // TBSiParameters
   try {

      const gear::GearParameters & params    = Global::GEAR->getTBSiParameters();

      _gearParamsGen = 0;
      _gearParams    = &params;
      _gearType      = "TBSiDet";
   }
   catch (gear::UnknownParameterException& e) {}

   // Error - no Gear parameters saved
   if (_gearType == "") {
      streamlog_out(ERROR) << "Couldn't read Gear parameters from given xml "
                           << "file, check if file of correct type!!! "
       	                  << std::endl;
      exit(0);
   }

   // No error - read variables
   else {

      //
      // Gear type: VXD
      if (_gearType == "VXD") {

         const gear::VXDParameters  * gearVXD  = dynamic_cast<const gear::VXDParameters*>(_gearParams);
         const gear::VXDLayerLayout & layerVXD = gearVXD->getVXDLayerLayout();

         // Define number of layers
         _numberOfLayers  = (short int)layerVXD.getNLayers();

         // Define layer parameters
         _layerRadius.resize(_numberOfLayers);
         _layerPhi0.resize(_numberOfLayers);
         _layerTheta  = _gearParamsGen->getDoubleVals("ActiveLayerTheta");

         _layerRealID = _gearParamsGen->getIntVals("ActiveLayerID");
         _layerType   = _gearParamsGen->getIntVals("ActiveLayerType");

         // Define ladder parameters
         _numberOfLadders.resize(_numberOfLayers);

         _ladderThick.resize(_numberOfLayers);
         _ladderWidth.resize(_numberOfLayers);
         _ladderLength.resize(_numberOfLayers);
         _ladderOffsetY.resize(_numberOfLayers);
         _ladderOffsetZ = _gearParamsGen->getDoubleVals("ActiveLayerOffsetZ");

         // Define sensor parameters
         _numberOfSensors      = _gearParamsGen->getIntVals("ActiveLadderNSensors");

         _sensorNPixelsInZ.resize(_numberOfLayers);
         _sensorNPixelsInRPhi.resize(_numberOfLayers);

         _sensorThick.resize(_numberOfLayers);

         _sensorWidth          = _gearParamsGen->getDoubleVals("ActiveSensorWidth");
         _sensorLength         = _gearParamsGen->getDoubleVals("ActiveSensorLength");
         _sensorGapInBtw       = _gearParamsGen->getDoubleVals("SensorGapInBetween");
         _sensorPitchInZ       = _gearParamsGen->getDoubleVals("ActiveSensorPadSizeZ");
         _sensorPitchInRPhi    = _gearParamsGen->getDoubleVals("ActiveSensorPadSizeRPhi");

         _sensorRimWidthInZ.resize(_numberOfLayers);
         _sensorRimWidthInRPhi.resize(_numberOfLayers);

         // Set parameters in correct units
         for (int iLayer = 0; iLayer < _numberOfLayers; iLayer++) {

            _numberOfLadders[iLayer]      = layerVXD.getNLadders(iLayer);

            _layerRadius[iLayer]          = layerVXD.getSensitiveDistance(iLayer) * mm;
            _layerPhi0[iLayer]            = layerVXD.getPhi0(iLayer)/180.         * pi;
            _layerTheta[iLayer]           = _layerTheta[iLayer]/180.              * pi;

            _ladderThick[iLayer]          = layerVXD.getSensitiveThickness(iLayer)* mm;
            _ladderWidth[iLayer]          = layerVXD.getSensitiveWidth(iLayer)    * mm;
            _ladderLength[iLayer]         = layerVXD.getSensitiveLength(iLayer)   * mm;
            _ladderOffsetY[iLayer]        = layerVXD.getSensitiveOffset(iLayer)   * mm;
            _ladderOffsetZ[iLayer]        = _ladderOffsetZ[iLayer]                * mm;

            _sensorThick[iLayer]          = _ladderThick[iLayer];
            _sensorWidth[iLayer]          = _sensorWidth[iLayer]                  * mm;
            _sensorLength[iLayer]         = _sensorLength[iLayer]                 * mm;

            _sensorGapInBtw[iLayer]       = _sensorGapInBtw[iLayer]               * mm;

            _sensorPitchInZ[iLayer]       = _sensorPitchInZ[iLayer]               * mm;
            _sensorPitchInRPhi[iLayer]    = _sensorPitchInRPhi[iLayer]            * mm;

            _sensorRimWidthInZ[iLayer]    = ((_ladderLength[iLayer] - (_numberOfSensors[iLayer]-1)*_sensorGapInBtw[iLayer])/_numberOfSensors[iLayer] -
                                            _sensorLength[iLayer])/2.;
            _sensorRimWidthInRPhi[iLayer] = (_ladderWidth[iLayer] - _sensorWidth[iLayer])/2.;

            if (_sensorPitchInZ[iLayer]    != 0.) _sensorNPixelsInZ[iLayer]    = floor(_sensorLength[iLayer]/_sensorPitchInZ[iLayer] + 0.5);
            else                                  _sensorNPixelsInZ[iLayer]    = 0.;
            if (_sensorPitchInRPhi[iLayer] != 0.) _sensorNPixelsInRPhi[iLayer] = floor(_sensorWidth[iLayer]/_sensorPitchInRPhi[iLayer] + 0.5);
            else                                  _sensorNPixelsInRPhi[iLayer] = 0.;
         }
      }

      //
      // Gear type: TBSiDet
      if (_gearType == "TBSiDet") {

         const gear::TBSiParameters  * gearTBSiDet  = dynamic_cast<const gear::TBSiParameters*>(_gearParams);
          
         // Define number of layers - for compatibility with cylindrical structure
         _numberOfLayers  = 1;

         // Define layer parameters - for compatibility with cylindrical structure
         _layerRealID.resize(1);
         _layerRealID[0] = 0;
         _layerType.resize(1);
         _layerType[0]   = 0;

         // Define ladder parameters - for compatibility with cylindrial structure
         _numberOfLadders.resize(1);
         _numberOfLadders[0] = 1;

         // Define sensor parameters
         _numberOfSensors.resize(1);
         _numberOfSensors[0] = gearTBSiDet->getNModules();

         _sensorWidth.resize(_numberOfSensors[0]);
         _sensorLength.resize(_numberOfSensors[0]);
         _sensorThick.resize(_numberOfSensors[0]);

         _sensorPosX.resize(_numberOfSensors[0]);
         _sensorPosY.resize(_numberOfSensors[0]);
         _sensorPosZ.resize(_numberOfSensors[0]);

         _sensorRotX.resize(_numberOfSensors[0]);

         _sensorPitchInRPhi.resize(_numberOfSensors[0]);
         _sensorPitchInZ.resize(_numberOfSensors[0]);

         _sensorNPixelsInRPhi.resize(_numberOfSensors[0]);
         _sensorNPixelsInZ.resize(_numberOfSensors[0]);

         for (int iSensor=0; iSensor<_numberOfSensors[0]; iSensor++) {

            _sensorWidth[iSensor]         = gearTBSiDet->getSensorSizeX(iSensor) * mm;
            _sensorLength[iSensor]        = gearTBSiDet->getSensorSizeY(iSensor) * mm;
            _sensorThick[iSensor]         = gearTBSiDet->getSensorSizeZ(iSensor) * mm;

            _sensorPosX[iSensor]          = gearTBSiDet->getSensorPosX(iSensor) * mm;
            _sensorPosY[iSensor]          = gearTBSiDet->getSensorPosY(iSensor) * mm;
            _sensorPosZ[iSensor]          = gearTBSiDet->getSensorPosZ(iSensor) * mm;

            _sensorRotX[iSensor]          = gearTBSiDet->getSensorRotX(iSensor)/180. * pi;

            _sensorPitchInRPhi[iSensor]   = gearTBSiDet->getSensorPadSizeX(iSensor) * mm;
            _sensorPitchInZ[iSensor]      = gearTBSiDet->getSensorPadSizeY(iSensor) * mm;

            _sensorNPixelsInRPhi[iSensor] = floor(_sensorWidth[iSensor]/_sensorPitchInRPhi[iSensor] + 0.5);
            _sensorNPixelsInZ[iSensor]    = floor(_sensorLength[iSensor]/_sensorPitchInZ[iSensor]  + 0.5);

         }

      }
   }

   // Print Gear parameters
   //printGearParams();
}


// GEOMETRY PROPERTIES

//
// Encode uniquely sensor ID
//
int SiPxlGeom::encodeSensorID(short int layerID, short int ladderID, short int sensorID) const
{
   return layerID*LAYERCOD + ladderID*LADDERCOD + sensorID*SENSORCOD;
}

//
// Decode unique sensor ID
//
void SiPxlGeom::decodeSensorID(short int & layerID, short int & ladderID, short int & sensorID, int uniqSensorID) const
{
   layerID  =  uniqSensorID / LAYERCOD;
   ladderID = (uniqSensorID - layerID*LAYERCOD) / LADDERCOD;
   sensorID = (uniqSensorID - layerID*LAYERCOD - ladderID*LADDERCOD) / SENSORCOD;
}

//
// Encode uniquely pixel ID
//
int SiPxlGeom::encodePixelID(short int layerID, short int sensorID, int pixelIDRPhi, int pixelIDZ) const
{
   return (SiPxlGeom::getSensorNPixelsInRPhi(layerID, sensorID)*pixelIDZ + pixelIDRPhi);
}

//
// Decode unique pixel ID
//
void SiPxlGeom::decodePixelID(int & pixelIDRPhi, int & pixelIDZ, short int layerID, short int sensorID, int uniqPixelID) const
{
   pixelIDZ    = uniqPixelID / SiPxlGeom::getSensorNPixelsInRPhi(layerID, sensorID);
   pixelIDRPhi = uniqPixelID - pixelIDZ*SiPxlGeom::getSensorNPixelsInRPhi(layerID, sensorID);
}

//
// Get number of PXD layers
//
short int SiPxlGeom::getNPXDLayers() const
{
   short int nLayers = 0;

   for (int i=0; i<getNLayers(); i++) {

      if (getLayerType(i) == pixel) nLayers++;
   }

   return nLayers;
}

//
// Get layer real ID
//
int SiPxlGeom::getLayerRealID(short int layerID) const
{
   if (_layerRealID.size()>(unsigned short int)layerID) return _layerRealID[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getLayerRealID - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Transform real layer ID to C-type numbering 0 - n ...
//
short int SiPxlGeom::getLayerIDCTypeNo(int realLayerID) const
{
   for (unsigned int i=0; i<_layerRealID.size(); i++) {

      if (realLayerID == _layerRealID[i]) return i;
   }

   streamlog_out(ERROR) << "SiPxlGeom::getLayerIDCTypeNo - layer: " << realLayerID << " not found in layer list!!!" << std::endl;

   exit(0);
}

//
// Get layer type
//
short int SiPxlGeom::getLayerType(short int layerID) const
{
   if (_layerType.size()>(unsigned short int)layerID) return _layerType[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getLayerType - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get layer radius
//
double SiPxlGeom::getLayerRadius(short int layerID) const
{
   if (_layerRadius.size()>(unsigned short int)layerID) return _layerRadius[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getLayerRadius - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get layer phi zero angle
//
double SiPxlGeom::getLayerPhi0(short int layerID) const
{
   if (_layerPhi0.size()>(unsigned short int)layerID) return _layerPhi0[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getLayerPhi0 - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get number of ladders
//
short int SiPxlGeom::getNLadders(short int layerID) const
{
   if (_numberOfLadders.size()>(unsigned short int)layerID) return _numberOfLadders[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getNLadders - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder thickness
//
double SiPxlGeom::getLadderThick(short int layerID) const
{
   if (_ladderThick.size()>(unsigned short int)layerID) return _ladderThick[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getLadderThick - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder width
//
double SiPxlGeom::getLadderWidth(short int layerID) const
{
   if (_ladderWidth.size()>(unsigned short int)layerID) return _ladderWidth[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getLadderWidth - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder length
//
double SiPxlGeom::getLadderLength(short int layerID) const
{
   if (_ladderLength.size()>(unsigned short int)layerID) return _ladderLength[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getLadderLength - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder offset in Y
//
double SiPxlGeom::getLadderOffsetY(short int layerID) const
{
   if (_ladderOffsetY.size()>(unsigned short int)layerID) return _ladderOffsetY[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getLadderOffsetY - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder offset in Z
//
double SiPxlGeom::getLadderOffsetZ(short int layerID) const
{
   if (_ladderOffsetZ.size()>(unsigned short int)layerID) return _ladderOffsetZ[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getLadderOffsetZ - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder rotation - phi angle (in system of units defined in PhysicalConstants.h)
//
double SiPxlGeom::getLadderPhi(short int layerID, short int ladderID) const
{
   if (ladderID<getNLadders(layerID)) return (getLayerPhi0(layerID) + 2*pi/getNLadders(layerID)*ladderID);
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getLadderPhi - ladderID: " << ladderID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder rotation - theta angle
//
double SiPxlGeom::getLadderTheta(short int layerID) const
{
   if (_layerTheta.size()>(unsigned short int)layerID) return _layerTheta[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getLadderTheta - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get number of sensors for given ladder
//
short int SiPxlGeom::getNSensors(short int layerID) const
{
   if (_numberOfSensors.size()>(unsigned short int)layerID) return _numberOfSensors[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getNSensors - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get number of pixels in Z axis (in each sensor)
//
int SiPxlGeom::getSensorNPixelsInZ(short int layerID, short int sensorID) const
{
   // Gear type: VXD
   if (_gearType == "VXD") {

      if (_sensorNPixelsInZ.size()>(unsigned short int)layerID) return _sensorNPixelsInZ[layerID];
      else {

         streamlog_out(ERROR) << "SiPxlGeom::getSensorNPixelsInZ - layerID: " << layerID << " out of range!!!"
                              << std::endl;
         exit(0);
      }
   }

   // Gear type: TBSiDet
   else if (_gearType == "TBSiDet") {

      if (_sensorNPixelsInZ.size()>(unsigned short int)sensorID) return _sensorNPixelsInZ[sensorID];
      else {

         streamlog_out(ERROR) << "SiPxlGeom::getSensorNPixelsInZ - sensorID: " << sensorID << " out of range!!!"
                              << std::endl;
         exit(0);
      }
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "SiPxlGeom::getSensorNPixelsInZ - unknown gear type!"
                           << std::endl;

      exit(0);
   }
}

//
// Get number of pixels in R-Phi (in each sensor)
//
int SiPxlGeom::getSensorNPixelsInRPhi(short int layerID, short int sensorID) const
{
   // Gear type: VXD
   if (_gearType == "VXD") {

      if (_sensorNPixelsInRPhi.size()>(unsigned short int)layerID) return _sensorNPixelsInRPhi[layerID];
      else {

         streamlog_out(ERROR) << "SiPxlGeom::getSensorNPixelsInRPhi - layerID: " << layerID << " out of range!!!"
                              << std::endl;
         exit(0);
      }
   }

   // Gear type: TBSiDet
   else if (_gearType == "TBSiDet") {

      if (_sensorNPixelsInRPhi.size()>(unsigned short int)sensorID) return _sensorNPixelsInRPhi[sensorID];
      else {

         streamlog_out(ERROR) << "SiPxlGeom::getSensorNPixelsInRPhi - sensorID: " << sensorID << " out of range!!!"
                              << std::endl;
         exit(0);
      }
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "SiPxlGeom::getSensorNPixelsInRPhi - unknown gear type!"
                           << std::endl;

      exit(0);
   }
}

//
// Get sensor pitch in Z axis
//
double SiPxlGeom::getSensorPitchInZ(short int layerID, short int sensorID) const
{
   // Gear type: VXD
   if (_gearType == "VXD") {

      if (_sensorPitchInZ.size()>(unsigned short int)layerID) return _sensorPitchInZ[layerID];
      else {

         streamlog_out(ERROR) << "SiPxlGeom::getSensorPitchInZ - layerID: " << layerID << " out of range!!!"
                              << std::endl;
         exit(0);
      }
   }

   // Gear type: TBSiDet
   else if (_gearType == "TBSiDet") {

      if (_sensorPitchInZ.size()>(unsigned short int)sensorID) return _sensorPitchInZ[sensorID];
      else {

         streamlog_out(ERROR) << "SiPxlGeom::getSensorPitchInZ - sensorID: " << sensorID << " out of range!!!"
                              << std::endl;
         exit(0);
      }
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "SiPxlGeom::getSensorPitchInZ - unknown gear type!"
                           << std::endl;

      exit(0);
   }
}

//
// Get sensor pitch in R-Phi
//
double SiPxlGeom::getSensorPitchInRPhi(short int layerID, short int sensorID) const
{
   // Gear type: VXD
   if (_gearType == "VXD") {

      if (_sensorPitchInRPhi.size()>(unsigned short int)layerID) return _sensorPitchInRPhi[layerID];
      else {

         streamlog_out(ERROR) << "SiPxlGeom::getSensorPitchInRPhi - layerID: " << layerID << " out of range!!!"
                              << std::endl;
         exit(0);
      }
   }

   // Gear type: TBSiDet
   else if (_gearType == "TBSiDet") {

      if (_sensorPitchInRPhi.size()>(unsigned short int)sensorID) return _sensorPitchInRPhi[sensorID];
      else {

         streamlog_out(ERROR) << "SiPxlGeom::getSensorPitchInRPhi - sensorID: " << sensorID << " out of range!!!"
                              << std::endl;
         exit(0);
      }
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "SiPxlGeom::getSensorPitchInRPhi - unknown gear type!"
                           << std::endl;

      exit(0);
   }
}

//
// Get sensor thickness
//
double SiPxlGeom::getSensorThick(short int layerID, short int sensorID) const
{
   // Gear type: VXD
   if (_gearType == "VXD") {

      if (_sensorThick.size()>(unsigned short int)layerID) return _sensorThick[layerID];
      else {

         streamlog_out(ERROR) << "SiPxlGeom::getSensorThick - layerID: " << layerID << " out of range!!!"
                              << std::endl;
         exit(0);
      }
   }

   // Gear type: TBSiDet
   else if (_gearType == "TBSiDet") {

      if (_sensorThick.size()>(unsigned short int)sensorID) return _sensorThick[sensorID];
      else {

         streamlog_out(ERROR) << "SiPxlGeom::getSensorThick - sensorID: " << sensorID << " out of range!!!"
                              << std::endl;
         exit(0);
      }
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "SiPxlGeom::getSensorThick - unknown gear type!"
                           << std::endl;

      exit(0);
   }
}

//
// Get sensor width
//
double SiPxlGeom::getSensorWidth(short int layerID, short int sensorID) const
{
   // Gear type: VXD
   if (_gearType == "VXD") {

      if (_sensorWidth.size()>(unsigned short int)layerID) return _sensorWidth[layerID];
      else {

         streamlog_out(ERROR) << "SiPxlGeom::getSensorWidth - layerID: " << layerID << " out of range!!!"
                              << std::endl;
         exit(0);
      }
   }

   // Gear type: TBSiDet
   else if (_gearType == "TBSiDet") {

      if (_sensorWidth.size()>(unsigned short int)sensorID) return _sensorWidth[sensorID];
      else {

         streamlog_out(ERROR) << "SiPxlGeom::getSensorWidth - sensorID: " << sensorID << " out of range!!!"
                              << std::endl;
         exit(0);
      }
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "SiPxlGeom::getSensorWidth - unknown gear type!"
                           << std::endl;

      exit(0);
   }
}

//
// Get sensor length
//
double SiPxlGeom::getSensorLength(short int layerID, short int sensorID) const
{
   // Gear type: VXD
   if (_gearType == "VXD") {

      if (_sensorLength.size()>(unsigned short int)layerID) return _sensorLength[layerID];
      else {

         streamlog_out(ERROR) << "SiPxlGeom::getSensorLength - layerID: " << layerID << " out of range!!!"
                              << std::endl;
         exit(0);
      }
   }

   // Gear type: TBSiDet
   else if (_gearType == "TBSiDet") {

      if (_sensorLength.size()>(unsigned short int)sensorID) return _sensorLength[sensorID];
      else {

         streamlog_out(ERROR) << "SiPxlGeom::getSensorLength - sensorID: " << sensorID << " out of range!!!"
                              << std::endl;
         exit(0);
      }
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "SiPxlGeom::getSensorLength - unknown gear type!"
                           << std::endl;

      exit(0);
   }
}

//
// Get sensor position in X
//
double SiPxlGeom::getSensorPosX(short int sensorID) const
{
   if (_sensorPosX.size()>(unsigned short int)sensorID) return _sensorPosX[sensorID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getSensorPosX - sensorID: " << sensorID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get sensor position in Y
//
double SiPxlGeom::getSensorPosY(short int sensorID) const
{
   if (_sensorPosY.size()>(unsigned short int)sensorID) return _sensorPosY[sensorID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getSensorPosY - sensorID: " << sensorID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get sensor position in Z
//
double SiPxlGeom::getSensorPosZ(short int sensorID) const
{
   if (_sensorPosZ.size()>(unsigned short int)sensorID) return _sensorPosZ[sensorID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getSensorPosZ - sensorID: " << sensorID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get sensor rotation in X
//
double SiPxlGeom::getSensorRotX(short int sensorID) const
{
   if (_sensorRotX.size()>(unsigned short int)sensorID) return _sensorRotX[sensorID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getSensorRotX - sensorID: " << sensorID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get gap size inbetween sensors
//
double SiPxlGeom::getSensorGapInBetween(short int layerID) const
{
   if (_sensorGapInBtw.size()>(unsigned short int)layerID) return _sensorGapInBtw[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getSensorGapInBetween - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get width of sensor rim in Z (passive part of silicon)
//
double SiPxlGeom::getSensorRimWidthInZ(short int layerID) const
{
   if (_sensorRimWidthInZ.size()>(unsigned short int)layerID) return _sensorRimWidthInZ[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getSensorRimWidthInZ - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get width of sensor rim in R-Phi (passive part of silicon)
//
double SiPxlGeom::getSensorRimWidthInRPhi(short int layerID) const
{
   if (_sensorRimWidthInRPhi.size()>(unsigned short int)layerID) return _sensorRimWidthInRPhi[layerID];
   else {

      streamlog_out(ERROR) << "SiPxlGeom::getSensorRimWidthInRPhi - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}


// TRANSFORMATION METHODS - GLOBAL REF. SYSTEM

//
// Method transforming given point from global ref. system to local ref. system
// (parameters: layerID, ladderID, sensorID and space point in global ref. system)
//
Hep3Vector SiPxlGeom::transformPointToLocal(short int layerID, short int ladderID, short int sensorID, const Hep3Vector & globalPoint)
{
   // Initialize local point
	Hep3Vector localPoint(globalPoint);

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Calculate rotation angles
       double theta = getLadderTheta(layerID);
       double phi   = getLadderPhi(layerID, ladderID);

       // Find (0,0,0) position of local coordinate system (defined in Geant4)
       Hep3Vector localOrigin(getLayerRadius(layerID), getLadderOffsetY(layerID), getLadderOffsetZ(layerID));
       localOrigin.rotateZ(+phi);

       // Perform translation - to the center of a ladder
       localPoint -= localOrigin;

       // Perform rotation - to the center of a ladder
       localPoint.rotateZ(-phi);
       localPoint.rotateY(+theta);

       // Perform translation to the centre of local coordinate system (defined in SiPxlGeom)
       localPoint += Hep3Vector(+getSensorThick(layerID, sensorID)/2., +getSensorWidth(layerID, sensorID)/2.,
                                +0.5*getLadderLength(layerID) - (2*sensorID + 1)*getSensorRimWidthInZ(layerID) - sensorID*getSensorGapInBetween(layerID)
                                -sensorID*getSensorLength(layerID, sensorID));

      // Check if local point within sensor boundaries +- epsilon
      if (isPointOutOfSensor(layerID, sensorID, localPoint)) {

         streamlog_out(ERROR) << std::setprecision(3) << "SiPxlGeom::transformPointToLocal - point: "
                              << localPoint           << " is out of sensor!!!"
                              << std::setprecision(0) << std::endl;
         exit(0);
      }
   }

   // Gear type: TBSiParameters
   else if (_gearType == "TBSiDet") {

      // Calculate rotation angles
      double phi = getSensorRotX(sensorID);

      // Find (0,0,0) position of local coordinate system (defined in Geant4)
      Hep3Vector localOrigin(getSensorPosX(sensorID), getSensorPosY(sensorID), getSensorPosZ(sensorID));

      // Perform translation - to the center of a sensor
      localPoint -= localOrigin;

      // Perform rotation - to the center of a sensor
      localPoint.rotateX(+phi);

      // Perform rotation Z->X, X->Y, Y->Z (Geant4 local to SiPxlGeom local)
      Hep3Vector tmpVec(localPoint);

      localPoint.setX(tmpVec.getZ());
      localPoint.setY(tmpVec.getX());
      localPoint.setZ(tmpVec.getY());

      // Perform translation to the centre of local coordinate system (defined in SiPxlGeom)
      localPoint += Hep3Vector(+getSensorThick(layerID, sensorID)/2., +getSensorWidth(layerID, sensorID)/2., +getSensorLength(layerID, sensorID)/2.);

      // Check if local point within sensor boundaries +- epsilon
      if (isPointOutOfSensor(layerID, sensorID, localPoint)) {

         streamlog_out(ERROR) << std::setprecision(3) << "SiPxlGeom::transformPointToLocal - point: "
                              << localPoint           << " is out of sensor!!!"
                              << std::setprecision(0) << std::endl;
         exit(0);
      }
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "SiPxlGeom::transformPointToLocal - unknown gear type!"
                           << std::endl;

      exit(0);
   }

   // Return space point in local ref. system
   return localPoint;
}

//
// Method transforming given vector from global ref. system to local ref. system
// (parameters: layerID, ladderID, sensorID and vector in global ref. system)
//
Hep3Vector SiPxlGeom::transformVecToLocal(short int layerID, short int ladderID, short int sensorID, const Hep3Vector & globalVec)
{
   // Initialize local vector
	Hep3Vector localVec(globalVec);

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Calculate rotation angles
      double theta = getLadderTheta(layerID);
      double phi   = getLadderPhi(layerID, ladderID);

      // Perform rotation - to the center of a ladder
      localVec.rotateZ(-phi);
      localVec.rotateY(+theta);

   }

   // Gear type: TBSiParameters
   else if (_gearType == "TBSiDet") {

      // Calculate rotation angles
      double phi = getSensorRotX(sensorID);

      // Perform rotation - to the center of a sensor
      localVec.rotateX(+phi);

      // Perform rotation Z->X, X->Y, Y->Z (Geant4 local to SiPxlGeom local)
      Hep3Vector tmpVec(localVec);

      localVec.setX(tmpVec.getZ());
      localVec.setY(tmpVec.getX());
      localVec.setZ(tmpVec.getY());
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "SiPxlGeom::transformVecToLocal - unknown gear type!"
                           << std::endl;

      exit(0);
   }

   // Return vector in local ref. system
   return localVec;
}

//
// Method transforming given matrix 3x3 from global ref. system to local ref. system
// (parameters: layerID, ladderID, sensorID and matrix in global ref. system)
//
HepMatrix SiPxlGeom::transformMatxToLocal(short int layerID, short int ladderID, short int sensorID, const HepMatrix & globalMatrix)
{
   // Initialize local matrix 3x3 to zero values
	HepMatrix localMatrix(3,3,0);

   // Initialize rotation matrices: R, R^T (transposition)
	HepMatrix rotMatrix(3,3,0);
	HepMatrix rotMatrixT(3,3,0);

	// Gear type: VXD
	if (_gearType == "VXD") {

	   // Calculate rotation angles
	   double theta = getLadderTheta(layerID);
	   double phi   = getLadderPhi(layerID, ladderID);

	   // Calculate rotation matrices - help
	   HepRotation rotMatrixZ(Hep3Vector(0,0,1),-phi);
	   HepRotation rotMatrixY(Hep3Vector(0,1,0),+theta);
	   HepRotation rotMatrixHelp(rotMatrixY*rotMatrixZ);

	   HepRotation rotMatrixHelpT(rotMatrixHelp.inverse());

	   for (int i=0; i<3; i++) {
	      for (int j=0; j<3; j++) {
	         rotMatrix[i][j]  = rotMatrixHelp[i][j];
	         rotMatrixT[i][j] = rotMatrixHelpT[i][j];
	      }
	   }
	}

	// Gear type: unknown - error
	else {
	   streamlog_out(ERROR) << "SiPxlGeom::transformMatxToLocal - unknown gear type!"
	                        << std::endl;

	   exit(0);
	}

   // Transform given matrix - rotation wrt global ref. system
   localMatrix = rotMatrix*globalMatrix*rotMatrixT;

   // Return matrix in local ref. system
   return localMatrix;
}


// TRANSFORMATION METHODS - LOCAL REF. SYSTEM

//
// Method transforming given point from local ref. system to global ref. system
// (parameters: layerID, ladderID, sensorID and space point in local ref. system)
//
Hep3Vector SiPxlGeom::transformPointToGlobal(short int layerID, short int ladderID, short int sensorID, const Hep3Vector & localPoint)
{
   // Initialize global point
   Hep3Vector globalPoint(localPoint);

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Calculate rotation angles
      double theta = getLadderTheta(layerID);
      double phi   = getLadderPhi(layerID, ladderID);

      // Find (0,0,0) position of local coordinate system (defined in Geant4)
      Hep3Vector localOrigin(getLayerRadius(layerID), getLadderOffsetY(layerID), getLadderOffsetZ(layerID));
      localOrigin.rotateZ(+phi);

      // Perform translation - to the center of a ladder
      globalPoint -= Hep3Vector(+getSensorThick(layerID, sensorID)/2., +getSensorWidth(layerID, sensorID)/2.,
                                +0.5*getLadderLength(layerID) - (2*sensorID + 1)*getSensorRimWidthInZ(layerID) - sensorID*getSensorGapInBetween(layerID)
                                -sensorID*getSensorLength(layerID, sensorID));

      // Perform rotation - with respect to the local centre of a ladder
      globalPoint.rotateY(-theta);
      globalPoint.rotateZ(+phi);

      // Perform translation - to the global system
      globalPoint += localOrigin;
   }

   // Gear type: TBSiParameters
   else if (_gearType == "TBSiDet") {

      // Calculate rotation angles
      double phi = getSensorRotX(sensorID);

      // Find (0,0,0) position of local coordinate system (defined in Geant4)
      Hep3Vector localOrigin(getSensorPosX(sensorID), getSensorPosY(sensorID), getSensorPosZ(sensorID));

      // Perform translation - to the centre of a sensor
      globalPoint -= Hep3Vector(+getSensorThick(layerID, sensorID)/2., +getSensorWidth(layerID, sensorID)/2., +getSensorLength(layerID, sensorID)/2.);

      // Perform rotation X->Z, Y->X, Z->Y (SiPxlGeom local to Geant4 local)
      Hep3Vector tmpVec(globalPoint);

      globalPoint.setX(tmpVec.getY());
      globalPoint.setY(tmpVec.getZ());
      globalPoint.setZ(tmpVec.getX());

      // Perform rotation - with respect to the local centre of a sensor
      globalPoint.rotateX(-phi);

      // Perform translation - to the global system
      globalPoint += localOrigin;
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "SiPxlGeom::transformPointToGlobal - unknown gear type!"
                           << std::endl;

      exit(0);
   }

   // Return space point in global ref. system
	return globalPoint;
}

//
// Method transforming given vector from local ref. system to global ref. system
// (parameters: layerID, ladderID, sensorID and vector in local ref. system)
//
Hep3Vector SiPxlGeom::transformVecToGlobal(short int layerID, short int ladderID, short int sensorID, const Hep3Vector & localVec)
{
   // Initialize global vector
	Hep3Vector globalVec(localVec);

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Calculate rotation angles
      double theta = getLadderTheta(layerID);
      double phi   = getLadderPhi(layerID, ladderID);

      // Perform rotation - to the global system
      globalVec.rotateY(-theta);
      globalVec.rotateZ(+phi);
   }

   // Gear type: TBSiParameters
   else if (_gearType == "TBSiDet") {

      // Calculate rotation angles
      double phi = getSensorRotX(sensorID);

      // Perform rotation X->Z, Y->X, Z->Y (SiPxlGeom local to Geant4 local)
      Hep3Vector tmpVec(globalVec);

      globalVec.setX(tmpVec.getY());
      globalVec.setY(tmpVec.getZ());
      globalVec.setZ(tmpVec.getX());

      // Perform rotation - to the global system
      globalVec.rotateX(-phi);
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "SiPxlGeom::transformVecToGlobal - unknown gear type!"
                           << std::endl;

      exit(0);
   }

   // Return vector in global ref. system
   return globalVec;
}

//
// Method transforming given matrix 3x3 from local ref. system to global ref. system
// (parameters: layerID, ladderID, sensorID and matrix in local ref. system)
//
HepMatrix SiPxlGeom::transformMatxToGlobal(short int layerID, short int ladderID, short int sensorID, const HepMatrix & localMatrix)
{
   // Initialize local matrix 3x3 to zero values
   HepMatrix globalMatrix(3,3,0);

   // Initialize rotation matrices: R, R^T (transposition)
   HepMatrix rotMatrix(3,3,0);
   HepMatrix rotMatrixT(3,3,0);

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Calculate rotation angles
      double theta = getLadderTheta(layerID);
      double phi   = getLadderPhi(layerID, ladderID);

      // Calculate rotation matrices - help
      HepRotation rotMatrixY(Hep3Vector(0,1,0),-theta);
      HepRotation rotMatrixZ(Hep3Vector(0,0,1),+phi);
      HepRotation rotMatrixHelp(rotMatrixZ*rotMatrixY);

      HepRotation rotMatrixHelpT(rotMatrixHelp.inverse());

      for (int i=0; i<3; i++) {
         for (int j=0; j<3; j++) {
            rotMatrix[i][j]  = rotMatrixHelp[i][j];
            rotMatrixT[i][j] = rotMatrixHelpT[i][j];
         }
      }
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "SiPxlGeom::transformMatxToGlobal - unknown gear type!"
                           << std::endl;

      exit(0);
   }

   // Transform given matrix - rotation wrt local ref. system
   globalMatrix = rotMatrix*localMatrix*rotMatrixT;

   // Return matrix in global ref. system
   return globalMatrix;
}


// OTHER METHODS - LOCAL REF. SYSTEM

//
// Get info whether the given point is out of Si sensor (parameters: layerID,
// space point in local ref. system)
//
bool SiPxlGeom::isPointOutOfSensor(short int layerID, short int sensorID, const Hep3Vector & point) const
{
	bool isOut = false;

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Boundary set +- epsilon
      if ( (point.getX() < -(EPS*um)) || (point.getX() > +(getSensorThick(layerID, sensorID) + EPS*um)) ||
           (point.getY() < -(EPS*um)) || (point.getY() > +(getSensorWidth(layerID, sensorID) + EPS*um)) ||
           (point.getZ() < -(EPS*um)) || (point.getZ() > +(getSensorLength(layerID,sensorID) + EPS*um)) ) isOut = true;

   }

   // Gear type: TBSiDet
   else if (_gearType == "TBSiDet") {

      // Boundary set +- epsilon
      if ( (point.getX() < -(EPS*um)) || (point.getX() > +(getSensorThick(layerID, sensorID) + EPS*um)) ||
           (point.getY() < -(EPS*um)) || (point.getY() > +(getSensorWidth(layerID, sensorID) + EPS*um)) ||
           (point.getZ() < -(EPS*um)) || (point.getZ() > +(getSensorLength(layerID,sensorID) + EPS*um)) ) isOut = true;

   }
   // Gear type: unknown - error
   else {
   	streamlog_out(ERROR) << "SiPxlGeom::isPointOutOfSensor - unknown gear type!"
   	                     << std::endl;

      exit(0);
   }

   // Return if out or not
   return isOut;
}

//
// Correct the given point if out of Si sensor (parameters: layerID,
// space point in local ref. system)
//
void SiPxlGeom::correctPointOutOfSensor(short int layerID, short int sensorID, Hep3Vector & point) const
{
   // Gear type: VXD
   if (_gearType == "VXD") {

      // Boundary set +- epsilon
      if (point.getX() < - EPS*um) point.setX(+ EPS*um);
      if (point.getX() > +(getSensorThick(layerID, sensorID) + EPS*um)) point.setX(getSensorThick(layerID, sensorID) - EPS*um);
      if (point.getY() < - EPS*um) point.setY(+ EPS*um);
      if (point.getY() > +(getSensorWidth(layerID, sensorID) + EPS*um)) point.setY(getSensorWidth(layerID, sensorID) - EPS*um);
      if (point.getZ() < - EPS*um) point.setZ(+ EPS*um);
      if (point.getZ() > +(getSensorLength(layerID, sensorID)+ EPS*um)) point.setZ(getSensorLength(layerID, sensorID)- EPS*um);
   }

   // Gear type: TBSiDet
   else if (_gearType == "TBSiDet") {

      // Boundary set +- epsilon
      if (point.getX() < - EPS*um) point.setX(+ EPS*um);
      if (point.getX() > +(getSensorThick(layerID, sensorID) + EPS*um)) point.setX(getSensorThick(layerID, sensorID) - EPS*um);
      if (point.getY() < - EPS*um) point.setY(+ EPS*um);
      if (point.getY() > +(getSensorWidth(layerID, sensorID) + EPS*um)) point.setY(getSensorWidth(layerID, sensorID) - EPS*um);
      if (point.getZ() < - EPS*um) point.setZ(+ EPS*um);
      if (point.getZ() > +(getSensorLength(layerID, sensorID)+ EPS*um)) point.setZ(getSensorLength(layerID, sensorID)- EPS*um);
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "Unknown gear type!"
                           << std::endl;

      exit(0);
   }

}

//
// Get Z-position of given pixel in local ref system (in system of units defined in PhysicalConstants.h);
//
double SiPxlGeom::getPixelPosInZ(short int layerID, short int sensorID, int pixelIDZ) const
{
	double posZ;

	// Get pitch
	double sensPitch = getSensorPitchInZ(layerID, sensorID);

	// Get length
	double sensLength = getSensorLength(layerID, sensorID);

   // Calculate position
	posZ  = sensPitch*(pixelIDZ + 0.5);

   // Error
   if ((posZ<0.) || (posZ>+sensLength)) {

      streamlog_out(ERROR) << "SiPxlGeom::getPixelPosInZ - position out of sensor!!!"
                           << std::endl;
      exit(0);

   }

	// Return Z position of given pixel in local ref. system
	return posZ;
}

//
// Get R-Phi position of given pixel in local ref system (in system of units defined in PhysicalConstants.h);
//
double SiPxlGeom::getPixelPosInRPhi(short int layerID, short int sensorID, bool bricked, int pixelIDRPhi, int pixelIDZ) const
{
   double posRPhi;

   // Get pitch
   double sensPitch = getSensorPitchInRPhi(layerID, sensorID);

   // Get width
   double sensWidth = getSensorWidth(layerID, sensorID);

   // Calculate position
   if (!bricked) posRPhi = sensPitch*(pixelIDRPhi + 0.5);

   else {

      // Even rows
      if (pixelIDZ%2 == 0) {

         posRPhi = sensPitch*(pixelIDRPhi + 0.5);
      }
      // Odd rows
      else {

         // First pixel
         if (pixelIDRPhi == 0) posRPhi = sensPitch*(0.25);

         // Other pixels
         else  {

            posRPhi = sensPitch*pixelIDRPhi;

            // Last pixel
            if ( (posRPhi>(sensWidth - 0.5*sensPitch)) && (posRPhi<(sensWidth + 0.5*sensPitch)) ) posRPhi = sensWidth - sensPitch*0.25;

         }
      }
   } // Calculate position


   // Error
   if ((posRPhi<0.) || (posRPhi>+sensWidth)) {

      streamlog_out(ERROR) << "SiPxlGeom::getPixelPosInRPhi - position out of sensor!!!"
                           << std::endl;
      exit(0);

   }

   // Return R-Phi position of given pixel in local ref. system
   return posRPhi;
}

//
// Get pixel ID (in Z), position is given in local ref. system;
//
int SiPxlGeom::getPixelIDInZ(short int layerID, short int sensorID, double posZ) const
{
	int    pixelID;

	// Get pitch
	double sensPitch = getSensorPitchInZ(layerID, sensorID);

	if (sensPitch == 0) {
      streamlog_out(ERROR) << "SiPxlGeom::getPixelIDInZ - division by zero (sensPitch is zero)!!!"
                           << std::endl;
      exit(0);

	}

	// Get number of Pixels
	int sensNPixels = getSensorNPixelsInZ(layerID, sensorID);

	// Calculate pixelID
	if (posZ <= 0.) pixelID = 0;
	else {

	   pixelID = floor( posZ / sensPitch );

	   if (pixelID >= sensNPixels) pixelID = sensNPixels - 1;
	}

	// Error
	if (pixelID >= sensNPixels) {

	   streamlog_out(ERROR) << "SiPxlGeom::getPixelIDInZ - pixelID in Z greater than number of pixels!!!"
	                        << std::endl;
	   exit(0);

	}

   // Return pixelID
   return pixelID;
}

//
// Get pixel ID (in R-Phi), position is given in local ref. system
//
int SiPxlGeom::getPixelIDInRPhi(short int layerID, short int sensorID, bool bricked, double posRPhi, double posZ) const
{
   int    pixelID;

   // Get pitch
   double sensPitch = getSensorPitchInRPhi(layerID, sensorID);

   if (sensPitch == 0) {
      streamlog_out(ERROR) << "SiPxlGeom::getPixelIDInRPhi - division by zero (sensPitch is zero)!!!"
                           << std::endl;
      exit(0);

   }

   // Get number of Pixels
   int sensNPixels = getSensorNPixelsInRPhi(layerID, sensorID);

   // Calculate pixelID
   if (posRPhi <= 0.) pixelID = 0;
   else {

      // Not bricked pixels
      if (!bricked) pixelID = floor(posRPhi / sensPitch);
      // Bricked pixels
      else {

         // Even rows
         if (getPixelIDInZ(layerID, sensorID, posZ)%2 == 0) pixelID = floor(posRPhi / sensPitch);

         // Odd rows
         else {

            // First pixels
            if (posRPhi <= sensPitch/2.) pixelID = 0;

            // Other pixels
            else pixelID = floor( (posRPhi + sensPitch/2.) / sensPitch);
         }

      }

      if (pixelID >= sensNPixels) pixelID = sensNPixels - 1;
   }

   // Error
   if (pixelID >= sensNPixels) {

      streamlog_out(ERROR) << "SiPxlGeom::getPixelIDInRPhi - pixelID in RPhi greater than number of pixels!!!"
                           << std::endl;
      exit(0);

   }

   // Return pixelID
   return pixelID;
}


// PRINT METHODS

//
// Method printing general Gear parameters
//
void SiPxlGeom::printGearParams() const
{
   streamlog_out(MESSAGE3) << std::endl
                           << " "
	                   		<< DUNDERL
                           << DBLUE
                           << "Gear parameters:"
                           << ENDCOLOR
                           << " "
                           << std::endl  << std::endl;

   // Gear type: BField
   streamlog_out(MESSAGE3) << "  B field [T]:       "
                           << _magField/T
                           << std::endl << std::endl;

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Print general info
      for (int i=0; i<_numberOfLayers; ++i){

         if (_layerType[i]==pixel) {

            streamlog_out(MESSAGE2) << std::endl;
            streamlog_out(MESSAGE2) << "  Layer: "           << _layerRealID[i]      << std::endl;
            streamlog_out(MESSAGE2) << "  LayerType:       " << "pixel"              << std::endl;
            streamlog_out(MESSAGE2) << "  NumberOfLadders: " << _numberOfLadders[i]  << std::endl;
            streamlog_out(MESSAGE2) << "  Radius[mm]:      " << _layerRadius[i]/mm   << std::endl;
            streamlog_out(MESSAGE2) << "  Width[mm]:       " << _ladderWidth[i]/mm   << std::endl;
            streamlog_out(MESSAGE2) << "  Length[mm]:      " << _ladderLength[i]/mm  << std::endl;
            streamlog_out(MESSAGE2) << "  Phi0:            " << _layerPhi0[i]        << std::endl;
            streamlog_out(MESSAGE2) << "  Theta:           " << _layerTheta[i]       << std::endl;
            streamlog_out(MESSAGE2) << "  OffsetY[mm]:     " << _ladderOffsetY[i]/mm << std::endl;
            streamlog_out(MESSAGE2) << "  OffsetZ[mm]:     " << _ladderOffsetZ[i]/mm << std::endl;

         }
      }
   }

   // Gear type: TBSiDet
   else if (_gearType == "TBSiDet") {

      // Print general info
      for (int i=0; i<_numberOfSensors[0]; ++i){

         streamlog_out(MESSAGE2) << std::endl;
         streamlog_out(MESSAGE2) << "  Sensor: "          << i                          << std::endl;
         streamlog_out(MESSAGE2) << "  SensorType:      " << "pixel"                    << std::endl;
         streamlog_out(MESSAGE2) << "  PitchInRPhi[um]: " << _sensorPitchInRPhi[i]/um   << std::endl;
         streamlog_out(MESSAGE2) << "  PitchInZ[um]:    " << _sensorPitchInZ[i]/um      << std::endl;
         streamlog_out(MESSAGE2) << "  NPixelsInRPhi:   " << _sensorNPixelsInRPhi[i]    << std::endl;
         streamlog_out(MESSAGE2) << "  NPixelsInZ:      " << _sensorNPixelsInZ[i]       << std::endl;
         streamlog_out(MESSAGE2) << "  Width[mm]:       " << _sensorWidth[i]/mm         << std::endl;
         streamlog_out(MESSAGE2) << "  Length[mm]:      " << _sensorLength[i]/mm        << std::endl;
         streamlog_out(MESSAGE2) << "  Thick[um]:       " << _sensorThick[i]/um         << std::endl;
         streamlog_out(MESSAGE2) << "  PosZ:            " << _sensorPosZ[i]/mm          << std::endl;
         streamlog_out(MESSAGE2) << "  RotX[deg]:       " << _sensorRotX[i]/pi*180      << std::endl;
      }
   }
   // Gear type: unknown - error
   else {
   	streamlog_out(ERROR) << "SiPxlGeom::printGearParams - unknown gear type!"
   	                     << std::endl;

      exit(0);
   }
}

//
// Method printing sensor Gear parameters (parameters: layerID)
//
void SiPxlGeom::printSensorParams(short int layerID) const
{

	// Gear type: VXD
   if (_gearType == "VXD") {

      // Print sensor parameters
      streamlog_out(MESSAGE2) << "    Parameters: " << _sensorThick[layerID]/um       << "um thick, "
                              << "with "            << _sensorPitchInZ[layerID]/um    << "um pitch "
                              << "and "             << _sensorNPixelsInZ[layerID]     << " pixels in Z"
                              << ", resp. "         << _sensorPitchInRPhi[layerID]/um << "um pitch "
                              << "and "             << _sensorNPixelsInRPhi[layerID]  << " pixels in R-Phi."
                              << std::endl;
   }
   // Gear type: unknown - error
   else {
   	streamlog_out(ERROR) << "SiPxlGeom::printSensorParams - unknown gear type!"
   	                     << std::endl;

      exit(0);
   }
}

} // Namespace;

