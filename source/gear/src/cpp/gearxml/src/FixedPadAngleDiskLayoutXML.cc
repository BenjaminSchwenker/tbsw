#include "gearxml/FixedPadAngleDiskLayoutXML.h"

#include "gearxml/XMLHandler.h"
// #include "gearxml/GearParametersXML.h"

#include "gearxml/tinyxml.h"
#include "gearimpl/FixedPadAngleDiskLayout.h"
#include "gear/GearMgr.h"

// #include <vector>
// #include <algorithm>
// #include <sstream>

namespace gear {
  

  TiXmlElement FixedPadAngleDiskLayoutXML::toXML( const PadRowLayout2D* layout ) const {

    // Check if type is FixedPadAngleDiskLayout
    
    const FixedPadAngleDiskLayout* fixedPadLayout=dynamic_cast<const FixedPadAngleDiskLayout*>(layout);
    
    if (fixedPadLayout==NULL) {
      
      throw ParseException("FixedPadAngleDiskLayoutXML::toXML wrong type !");
    }
    
    // append data to PadRowLayout2D
    TiXmlElement padRowLayout2DXML("PadRowLayout2D");    
    
    padRowLayout2DXML.SetAttribute("type","FixedPadAngleDiskLayout"); 
    
    padRowLayout2DXML.SetDoubleAttribute("rMin",(fixedPadLayout->getPlaneExtent())[0]);
    
    padRowLayout2DXML.SetDoubleAttribute("rMax",(fixedPadLayout->getPlaneExtent())[1]);
    
    padRowLayout2DXML.SetAttribute("nRow",fixedPadLayout->getNRows());
    
    padRowLayout2DXML.SetDoubleAttribute("phiMin", (fixedPadLayout->getPlaneExtent())[2]);
    
    padRowLayout2DXML.SetDoubleAttribute("phiMax", (fixedPadLayout->getPlaneExtent())[3]);
    
    padRowLayout2DXML.SetAttribute("nPadInRow",fixedPadLayout->getNPadsInRow());
    
    return padRowLayout2DXML ;
  }
    
    
  PadRowLayout2D* FixedPadAngleDiskLayoutXML::fromXML( const TiXmlElement* layout ) const {


    double rMin       =  atof(  getXMLAttribute( layout , "rMin" ) .c_str() ) ;
    double rMax       =  atof(  getXMLAttribute( layout , "rMax" ) .c_str() ) ;
    int nRow       =  atoi(  getXMLAttribute( layout , "nRow" ) .c_str() ) ;
    double phiMin     =  atof(  getXMLAttribute( layout , "phiMin" ) .c_str() ) ;
    int    phiMax     =  atoi(  getXMLAttribute( layout , "phiMax" ) .c_str() ) ;
    int nPadsInRow =  atoi(  getXMLAttribute( layout , "nPadsInRow" ) .c_str() ) ;
    
    
    return new FixedPadAngleDiskLayout( rMin, rMax, nRow, phiMin, phiMax, nPadsInRow ) ;
  }

}

