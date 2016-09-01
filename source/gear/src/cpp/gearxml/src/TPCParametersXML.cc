#include "gearxml/TPCParametersXML.h"

#include "gearxml/XMLHandlerMgr.h"
#include "gearxml/GearParametersXML.h"
#include "gearxml/PadRowLayout2DXML.h"

#include "gearxml/tinyxml.h"
#include "gearimpl/TPCParametersImpl.h"
#include "gearimpl/FixedPadSizeDiskLayout.h"
#include "gear/GearMgr.h"

#include <vector>
#include <algorithm>
#include <sstream>

 
#include <typeinfo>


namespace gear {


  TiXmlElement TPCParametersXML::toXML( const GearParameters & parameters ) const {

    // Check whether parameters is valid TPCParameters
    
    const TPCParameters* param=dynamic_cast<const TPCParameters*>(&parameters);

    if(param==NULL) {
            
      throw ParseException ("XML TPCParametersXML::toXML given parameter not of correct type. " 
			    "Needs to be a gear::TPCParameter");
    }
    
    // Set up TPC Detector as Element
    TiXmlElement det("detector");
    
    
    TiXmlElement driftVelocity("driftVelocity");
    driftVelocity.SetDoubleAttribute("value", param->getDriftVelocity());
    
    
    TiXmlElement maxDriftLength("maxDriftLength");
    maxDriftLength.SetDoubleAttribute("value", param->getMaxDriftLength());
    
    
    TiXmlElement readoutFrequency("readoutFrequency");
    readoutFrequency.SetDoubleAttribute("value", param->getReadoutFrequency());


    // Set up PadRowLayout2D
    const PadRowLayout2D* padLayout = &( param->getPadLayout() ) ;
  
    PadRowLayout2DXML* layoutXML = PadRowLayout2DXML::getHandler( padLayout ) ;

    if( layoutXML == 0 ) {

      throw ParseException( "TPCParametersXML::toXML: no handler for " + std::string(typeid(padLayout).name()) + " found !" ) ;
    }
    
    TiXmlElement padRowLayout2DXML = layoutXML->toXML( padLayout ) ;


    // Assemble all items to detector
    det.InsertEndChild(driftVelocity);
    
    det.InsertEndChild(maxDriftLength);

    det.InsertEndChild(readoutFrequency);

    det.InsertEndChild(padRowLayout2DXML);


    // Write all other parameters to detecotor as attributes
    GearParametersXML::getXMLForParameters( &det ,  param ) ;

    
    // return XMLElement
    return det;
     
 }
  
  
  GearParameters* TPCParametersXML::fromXML( const TiXmlElement* xmlElement, GearMgr* gearMgr) const {
    
    
    TPCParametersImpl* tpcParams = new TPCParametersImpl ;
    
    // first we read the generic parameters
    GearParametersXML::setParametersFromXML( xmlElement, tpcParams  ) ;
    
    tpcParams->setMaxDriftLength(   atof(  getChildElementValue( xmlElement, "maxDriftLength"   ).c_str() )) ;
    tpcParams->setDriftVelocity(    atof(  getChildElementValue( xmlElement, "driftVelocity"    ).c_str() )) ;
    tpcParams->setReadoutFrequency( atof(  getChildElementValue( xmlElement, "readoutFrequency" ).c_str() )) ;
    

    // get the pad layout
    const TiXmlElement* layout = xmlElement->FirstChildElement( "PadRowLayout2D" ) ;

    if( layout == 0 ) {
      
      std::stringstream str ;
      str  << "XMLParser::getChildElementValue missing element \"PadRowLayout2D\" " 
	   << " in element <" << xmlElement->Value() << "/> " ;
      
      throw ParseException( str.str() ) ;
    }

    std::string layoutType = getXMLAttribute( layout , "type" )   ;
  
    PadRowLayout2DXML* layoutXML = PadRowLayout2DXML::getHandler( layoutType ) ;

    if( layoutXML == 0 ) {

      throw ParseException( "TPCParametersXML::fromXML: no handler for " + layoutType + " found !" ) ;
    }

    PadRowLayout2D* dLayout = layoutXML->fromXML( layout ) ;


//     <PadRowLayout2D  type="FixedPadSizeDiskLayout" rMin="320.0" rMax="1680.0" padHeight="6.0" padWidth="2.0" />

//     double rMin      =  atof(  getXMLAttribute( layout , "rMin" ) .c_str() ) ;
//     double rMax      =  atof(  getXMLAttribute( layout , "rMax" ) .c_str() ) ;
//     double padHeight =  atof(  getXMLAttribute( layout , "padHeight" ) .c_str() ) ;
//     double padWidth  =  atof(  getXMLAttribute( layout , "padWidth" ) .c_str() ) ;
//     int    maxRow    =  atoi(  getXMLAttribute( layout , "maxRow" ) .c_str() ) ;
//     double padGap    =  atof(  getXMLAttribute( layout , "padGap" ) .c_str() ) ;

    
//     PadRowLayout2D* dLayout = 
//       new FixedPadSizeDiskLayout( rMin, rMax, padHeight, padWidth, maxRow, padGap ) ;
    
    tpcParams->setPadLayout( dLayout ) ;
    
    
    if( gearMgr != 0 ) {
      
      gearMgr->setTPCParameters( tpcParams ) ;
    }
    
    
    return tpcParams ;
  }
  
  
} // namespace
