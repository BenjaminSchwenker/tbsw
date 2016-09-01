#include "gearimpl/TBSiParametersImpl.h"

namespace gear {

//
// Method adding new module
//
void TBSiParametersImpl::addModule(const short & id             , const std::string & name,
                                   const double & sensorSizeX   , const double & sensorSizeY    , const double & sensorSizeZ,
                                   const double & sensorPosX    , const double & sensorPosY     , const double & sensorPosZ ,
                                   const double & sensorRotX    , const double & sensorRadLength,
                                   const double & sensorPadSizeX, const double & sensorPadSizeY ,
                                   const double & windowSizeX   , const double & windowSizeY    , const double & windowSizeZ,
                                   const double & windowRelPosZ , const double & windowRadLength)
{
   Module module;

// Set sensor properties
   module.id              = id;
   module.name            = name;

   module.sensorSizeX     = sensorSizeX;
   module.sensorSizeY     = sensorSizeY;
   module.sensorSizeZ     = sensorSizeZ;
   module.sensorPosX      = sensorPosX;
   module.sensorPosY      = sensorPosY;
   module.sensorPosZ      = sensorPosZ;
   module.sensorRotX      = sensorRotX;
   module.sensorRadLength = sensorRadLength;
   module.sensorPadSizeX  = sensorPadSizeX;
   module.sensorPadSizeY  = sensorPadSizeY;

   module.windowSizeX     = windowSizeX;
   module.windowSizeY     = windowSizeY;
   module.windowSizeZ     = windowSizeZ;
   module.windowRelPosZ   = windowRelPosZ;
   module.windowRadLength = windowRadLength;

// Save sensor in vector
   _moduleVec.push_back(module);
}

//
// Method adding new module
//
void TBSiParametersImpl::addModule(const Module & rModule)
{
   Module module = rModule;

   // Save sensor in vector
   _moduleVec.push_back(module);
}

} // Namespace gear

