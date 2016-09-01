#ifndef GEAR_TBSIPARAMETERSXML_H
#define GEAR_TBSIPARAMETERSXML_H 1

#include "gearxml/XMLHandler.h"
#include "gear/GearParameters.h"
#include "gear/GearMgr.h"

namespace gear {

class TiXmlNode;

//! XML handler for test beam geometry TBSiDet
//!
//! @author Z.Drasal, Charles University Prague
//! @version 00
//!

class TBSiParametersXML : public XMLHandler {

 public:

   //! Create an XML node from the given parameters
   virtual TiXmlElement toXML(const GearParameters & parameters) const;

   //! Extract GearParameters subclass from the given XML element (node)
   virtual GearParameters * fromXML(const TiXmlElement * xmlElement, GearMgr * gearMgr=0) const;

}; // Class

} // Namespace

#endif // GEAR_TBSIPARAMETERSXML_H
