#include "gearxml/FixedPadSizeDiskLayoutXML.h"

#include "gearxml/XMLHandler.h"
// #include "gearxml/GearParametersXML.h"

#include "gearxml/tinyxml.h"
#include "gearimpl/FixedPadSizeDiskLayout.h"
#include "gear/GearMgr.h"

// #include <vector>
// #include <algorithm>
// #include <sstream>

namespace gear {
  

  TiXmlElement FixedPadSizeDiskLayoutXML::toXML( const PadRowLayout2D* layout ) const {

    // Check if type is FixedPadSizeDiskLayout
    
    const FixedPadSizeDiskLayout* fixedPadLayout=dynamic_cast<const FixedPadSizeDiskLayout*>(layout);
    
    if (fixedPadLayout==NULL) {
      
      throw ParseException("FixedPadSizeDiskLayoutXML::toXML wrong type !");
    }
    
    // append data to PadRowLayout2D
    TiXmlElement padRowLayout2DXML("PadRowLayout2D");    
    
    padRowLayout2DXML.SetAttribute("type","FixedPadSizeDiskLayout"); 
    
    padRowLayout2DXML.SetDoubleAttribute("rMin",(fixedPadLayout->getPlaneExtent())[0]);
    
    padRowLayout2DXML.SetDoubleAttribute("rMax",(fixedPadLayout->getPlaneExtent())[1]);
    
    padRowLayout2DXML.SetDoubleAttribute("padHeight",fixedPadLayout->getPadHeight(1));
    
    padRowLayout2DXML.SetDoubleAttribute("padWidth", fixedPadLayout->getFixedPadWidth());
    
    padRowLayout2DXML.SetAttribute("maxRow",fixedPadLayout->getNRows());
    
    padRowLayout2DXML.SetDoubleAttribute("padGap",fixedPadLayout->getPadGap());
    
    return padRowLayout2DXML ;
  }
    
    
  PadRowLayout2D* FixedPadSizeDiskLayoutXML::fromXML( const TiXmlElement* layout ) const {


    double rMin      =  atof(  getXMLAttribute( layout , "rMin" ) .c_str() ) ;
    double rMax      =  atof(  getXMLAttribute( layout , "rMax" ) .c_str() ) ;
    double padHeight =  atof(  getXMLAttribute( layout , "padHeight" ) .c_str() ) ;
    double padWidth  =  atof(  getXMLAttribute( layout , "padWidth" ) .c_str() ) ;
    int    maxRow    =  atoi(  getXMLAttribute( layout , "maxRow" ) .c_str() ) ;
    double padGap    =  atof(  getXMLAttribute( layout , "padGap" ) .c_str() ) ;
    
    
    return new FixedPadSizeDiskLayout( rMin, rMax, padHeight, padWidth, maxRow, padGap ) ;
  }

}

