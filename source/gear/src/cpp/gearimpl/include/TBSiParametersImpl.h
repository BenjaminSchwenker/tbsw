#ifndef GEAR_TBSIPARAMETERSIMPL_H
#define GEAR_TBSIPARAMETERSIMPL_H 1

#include "gear/TBSiParameters.h"
#include "gearimpl/GearParametersImpl.h"

#include <vector>
#include <string>

namespace gear {

//! Concrete implementation of test beam (TBSiDet) layout
//!
//! @author Z.Drasal, Charles University Prague
//! @version 00
//!

class TBSiParametersImpl : public GearParametersImpl, public TBSiParameters {

 public:

//! Structure that defines module properties
    struct Module {

       short id;              // Module ID
       std::string name;      // Name

       double sensorSizeX;    // Sensor size in X (mm)
       double sensorSizeY;    // Sensor size in Y (mm)
       double sensorSizeZ;    // Sensor size in Z (mm)

       double sensorPosX;     // Sensor position X in (mm)
       double sensorPosY;     // Sensor position Y in (mm)
       double sensorPosZ;     // Sensor position Z in mm

       double sensorRotX;     // Sensor rotation angle - around X axis

       double sensorRadLength;// Sensor radiation length in mm

       double sensorPadSizeX; // Sensor pad size in X in mm
       double sensorPadSizeY; // Sensor pad size in Y in mm

       double windowSizeX;    // Window size in X (mm)
       double windowSizeY;    // Window size in Y (mm)
       double windowSizeZ;    // Window size in Z (mm)

       double windowRelPosZ;  // Window relative position (wrt to sensor) in Z (mm)

       double windowRadLength;// Window radiation length (mm)
   };

   typedef std::vector<Module> ModuleVec;

//! Constructor
    TBSiParametersImpl(std::string setupName) {_setupName = setupName;}

//! Destructor
    virtual ~TBSiParametersImpl() {}

//! Get setup name
    virtual std::string getSetupName() const {return _setupName;}

//! Get total number of modules
    virtual short getNModules() const {return _moduleVec.size();}

//! Get module name
    virtual std::string getModuleName(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).name; else return "";}

//! Get sensor size in X in mm
    virtual double getSensorSizeX(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).sensorSizeX; else return 0;}

//! Get sensor size in Y in mm
    virtual double getSensorSizeY(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).sensorSizeY; else return 0;}

//! Get sensor size in Z in mm
    virtual double getSensorSizeZ(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).sensorSizeZ; else return 0;}

//! Get sensor position X in mm
    virtual double getSensorPosX(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).sensorPosX; else return 0;}

//! Get sensor position Y in mm
    virtual double getSensorPosY(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).sensorPosY; else return 0;}

//! Get sensor position Z in mm
    virtual double getSensorPosZ(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).sensorPosZ; else return 0;}

//! Get sensor rotation along X axis (in degrees)
    virtual double getSensorRotX(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).sensorRotX; else return 0;}

//! Get sensor radiation length (in mm)
    virtual double getSensorRadLength(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).sensorRadLength; else return 0;}

//! Get sensor pad size in X (in mm)
    virtual double getSensorPadSizeX(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).sensorPadSizeX; else return 0;}

//! Get sensor pad size in Y (in mm)
    virtual double getSensorPadSizeY(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).sensorPadSizeY; else return 0;}

//! Get window size in X in mm
    virtual double getWindowSizeX(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).windowSizeX; else return 0;}

//! Get window size in Y in mm
    virtual double getWindowSizeY(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).windowSizeY; else return 0;}

//! Get window size in Z in mm
    virtual double getWindowSizeZ(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).windowSizeZ; else return 0;}

//! Get window radiation length (in mm)
    virtual double getWindowRadLength(short moduleID) const {
       if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).windowRadLength; else return 0;}

//! Get window relative position in Z (wrt sensor) (in mm)
    virtual double getWindowRelPosZ(short moduleID) const {
        if ((unsigned short)moduleID<=_moduleVec.size()) return _moduleVec.at(moduleID).windowRelPosZ; else return 0;}

//! Add new module to the given position
    void addModule(const short & id             , const std::string & name,
                   const double & sensorSizeX   , const double & sensorSizeY    , const double & sensorSizeZ,
                   const double & sensorPosX    , const double & sensorPosY     , const double & sensorPosZ ,
                   const double & sensorRotX    , const double & sensorRadLength,
                   const double & sensorPadSizeX, const double & sensorPadSizeY ,
                   const double & windowSizeX   , const double & windowSizeY    , const double & windowSizeZ,
                   const double & windowRelPosZ , const double & windowRadLength);

    void addModule(const Module & rModule);

 protected:

   std::string _setupName; //!< setup name
   ModuleVec   _moduleVec; //!< vector of individual modules

}; // Class

} // Namespace gear

#endif // GEAR_TBSIPARAMATERSIMPL_H
