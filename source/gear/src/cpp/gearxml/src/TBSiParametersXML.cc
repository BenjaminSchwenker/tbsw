#include "gearxml/TBSiParametersXML.h"

#include "gearxml/XMLHandlerMgr.h"
#include "gearxml/GearParametersXML.h"
#include "gearxml/tinyxml.h"
#include "gear/TBSiParameters.h"
#include "gearimpl/TBSiParametersImpl.h"

#include <string>
#include <sstream>
#include <iostream>
#include <vector>

namespace gear {

//
// Method creating an XML node from the given parameters
//
TiXmlElement TBSiParametersXML::toXML(const GearParameters & parameters) const
{

//
// Check if parameters are of correct type
   const TBSiParameters * param = dynamic_cast<const TBSiParameters *> (&parameters);

   if (param == 0) {

      throw Exception( "TBSiParametersXML::toXML: Given parameters are not of correct type!");
   }
//
// XML node - general info about testbeam setup
   TiXmlElement detector("detector");
   TiXmlElement setup("setup");

   // Attributes
   setup.SetAttribute("name"    , param->getSetupName());
   setup.SetAttribute("nModules", param->getNModules() );

//
// XML nodes - individual modules
   for (int i=0; i<(param->getNModules()); i++) {

      // Creat XML node and set attributes
      TiXmlElement module("module");
      module.SetAttribute("name", param->getModuleName(i));
      module.SetAttribute("id"  , i                       );

      //
      // Describe sensor features
      //
      // Sensor
      TiXmlElement sensor("sensor");
      sensor.SetDoubleAttribute("radLength", param->getSensorRadLength(i));
      sensor.SetDoubleAttribute("padSizeX" , param->getSensorPadSizeX(i) );
      sensor.SetDoubleAttribute("padSizeY" , param->getSensorPadSizeY(i) );

      // Sensor Dimensions
      TiXmlElement sensorDim("dimensions");
      sensorDim.SetDoubleAttribute("sizeX", param->getSensorSizeX(i));
      sensorDim.SetDoubleAttribute("sizeY", param->getSensorSizeY(i));
      sensorDim.SetDoubleAttribute("sizeZ", param->getSensorSizeZ(i));
      sensor.InsertEndChild(sensorDim);

      // Sensor position and rotation
      TiXmlElement sensorPos("position");
      sensorPos.SetDoubleAttribute("posX", param->getSensorPosX(i));
      sensorPos.SetDoubleAttribute("posY", param->getSensorPosY(i));
      sensorPos.SetDoubleAttribute("posZ", param->getSensorPosZ(i));
      sensorPos.SetDoubleAttribute("rotX", param->getSensorRotX(i) );
      sensor.InsertEndChild(sensorPos);

      module.InsertEndChild(sensor);

      // Window
      TiXmlElement window("window");
      window.SetDoubleAttribute("radLength", param->getWindowRadLength(i));

      // Window Dimensions
      TiXmlElement winDim("dimensions");
      winDim.SetDoubleAttribute("sizeX", param->getWindowSizeX(i));
      winDim.SetDoubleAttribute("sizeY", param->getWindowSizeY(i));
      winDim.SetDoubleAttribute("sizeZ", param->getWindowSizeZ(i));
      window.InsertEndChild(winDim);

      // Window position and rotation
      TiXmlElement winPos("position");
      winPos.SetDoubleAttribute("relPosZ", param->getWindowRelPosZ(i));
      window.InsertEndChild(winPos);

      module.InsertEndChild(window);

      setup.InsertEndChild(module);

   } // For

// Create an XML file
   GearParametersXML::getXMLForParameters( &setup , &parameters ) ;

   detector.InsertEndChild(setup);

// Return
   return detector;

}

//
// Method extracting GearParameters subclass from the given XML element (node)
//
GearParameters * TBSiParametersXML::fromXML(const TiXmlElement * xmlElement, GearMgr * gearMgr) const {

   // XML value
   std::stringstream value;

   // Obtain basic information about setup (name, ...)
   const TiXmlElement * xmlSetup = xmlElement->FirstChildElement("setup");
   std::string setupName     =  getXMLAttribute(xmlSetup, "name");

   value.str("")     ; value << getXMLAttribute(xmlSetup, "nModules");
   short nModules; value >> nModules;

   // Create parameters - TBSiParameters
   TBSiParametersImpl * paramImpl = new TBSiParametersImpl(setupName);

   // Obtain basic information about modules
   const TiXmlNode * xmlModule = 0;

   while((xmlModule = xmlSetup->IterateChildren("module" , xmlModule)) != 0 ) {

      // Name & id
      value.clear(); value.str(""); value << getXMLAttribute(xmlModule, "name") << std::ends;
      std::string moduleName      ; value >> moduleName;
      value.clear(); value.str(""); value << getXMLAttribute(xmlModule, "id")   << std::ends;
      short moduleID              ; value >> moduleID;

      // Sensor
      const TiXmlElement * xmlSensor = xmlModule->FirstChildElement("sensor");

      value.clear(); value.str(""); value << getXMLAttribute(xmlSensor, "radLength") << std::ends;
      double sensorRadLength      ; value >> sensorRadLength;
      value.clear(); value.str(""); value << getXMLAttribute(xmlSensor, "padSizeX")  << std::ends;
      double sensorPadSizeX       ; value >> sensorPadSizeX;
      value.clear(); value.str(""); value << getXMLAttribute(xmlSensor, "padSizeY")  << std::ends;
      double sensorPadSizeY       ; value >> sensorPadSizeY;

      const TiXmlElement * xmlSensorDim = xmlSensor->FirstChildElement("dimensions");
      value.clear(); value.str(""); value << getXMLAttribute(xmlSensorDim, "sizeX")  << std::ends;
      double sensorSizeX          ; value >> sensorSizeX;
      value.clear(); value.str(""); value << getXMLAttribute(xmlSensorDim, "sizeY")  << std::ends;
      double sensorSizeY          ; value >> sensorSizeY;
      value.clear(); value.str(""); value << getXMLAttribute(xmlSensorDim, "sizeZ")  << std::ends;
      double sensorSizeZ          ; value >> sensorSizeZ;

      const TiXmlElement * xmlSensorPos = xmlSensor->FirstChildElement("position");
      value.clear(); value.str(""); value << getXMLAttribute(xmlSensorPos, "posX")   << std::ends;
      double sensorPosX     ; value >> sensorPosX;
      value.clear(); value.str(""); value << getXMLAttribute(xmlSensorPos, "posY")   << std::ends;
      double sensorPosY           ; value >> sensorPosY;
      value.clear(); value.str(""); value << getXMLAttribute(xmlSensorPos, "posZ")   << std::ends;
      double sensorPosZ           ; value >> sensorPosZ;
      value.clear(); value.str(""); value << getXMLAttribute(xmlSensorPos, "rotX")   << std::ends;
      double sensorRotX           ; value >> sensorRotX;

      // Window
      const TiXmlElement * xmlWindow = xmlModule->FirstChildElement("window");

      value.clear(); value.str(""); value << getXMLAttribute(xmlWindow, "radLength") << std::ends;
      double windowRadLength; value >> windowRadLength;

      const TiXmlElement * xmlWindowDim = xmlWindow->FirstChildElement("dimensions");
      value.clear(); value.str(""); value << getXMLAttribute(xmlWindowDim, "sizeX") << std::ends;
      double windowSizeX          ; value >> windowSizeX;
      value.clear(); value.str(""); value << getXMLAttribute(xmlWindowDim, "sizeY") << std::ends;
      double windowSizeY          ; value >> windowSizeY;
      value.clear(); value.str(""); value << getXMLAttribute(xmlWindowDim, "sizeZ") << std::ends;
      double windowSizeZ          ; value >> windowSizeZ;

      const TiXmlElement * xmlWindowPos = xmlWindow->FirstChildElement("position");
      value.clear(); value.str(""); value << getXMLAttribute(xmlWindowPos, "relPosZ") << std::ends;
      double windowRelPosZ        ; value >> windowRelPosZ;

      // Use all these parameters and set a new module
      paramImpl->addModule(moduleID  , moduleName     , sensorSizeX   , sensorSizeY   , sensorSizeZ, sensorPosX , sensorPosY , sensorPosZ,
                           sensorRotX, sensorRadLength, sensorPadSizeX, sensorPadSizeY, windowSizeX, windowSizeY, windowSizeZ, windowRelPosZ,
                           windowRadLength);

   } // Module loop

   // Add parameters to GearMgr
   if( gearMgr != 0 ) gearMgr->setTBSiParameters(paramImpl);

   return paramImpl;
}

} // Namespace


