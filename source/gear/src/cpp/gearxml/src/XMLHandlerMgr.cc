#include "gearxml/XMLHandlerMgr.h"
#include "gearxml/GearParametersXML.h"
#include "gearxml/SiPlanesParametersXML.h"
#include "gearxml/tinyxml.h"

#include "gear/GEAR.h"


#include <string>
#include <sstream>

namespace gear{

  //  XMLHandlerMgr* XMLHandlerMgr::_me  = 0 ;



  XMLHandlerMgr* XMLHandlerMgr::instance() {

    static XMLHandlerMgr _me ;
      return &_me ;

    //     if( _me == 0 ) {
    //       _me = new XMLHandlerMgr ;
    //     }
    //     return _me ;
  }

  XMLHandlerMgr::~XMLHandlerMgr() {

    XMLHandlerMap::iterator it_end = _map.end() ;

    for( XMLHandlerMap::iterator it = _map.begin() ; it != it_end ; ++ it  ) {

      delete it->second ;
    }
  }


  XMLHandlerMgr::XMLHandlerMgr() {

    _map[ GEAR::GEARPARAMETERS ] = new GearParametersXML ;
    _map[ GEAR::SIPLANESPARAMETERS ] = new SiPlanesParametersXML ;

  }


  const XMLHandler* XMLHandlerMgr::getHandler( const std::string& type ) {

    return _map[ type ] ;
  }


  const void XMLHandlerMgr::setHandler( const std::string& type , XMLHandler* handler ) {

    if( handler == 0 )  // don't allow null pointer
      return  ;

    XMLHandlerMap::iterator it = _map.find( type ) ;

    if( it != _map.end()  ) {

      delete it->second ;

      it->second = handler ;

    } else {

      _map.insert( std::make_pair( type , handler )  ) ;
    }

  }



  // helper function (declared in XMLHandler.h )

  std::string getXMLAttribute(const TiXmlNode* node , const std::string& name ) {

    const TiXmlElement* el = node->ToElement() ;
    if( el == 0 )
      throw ParseException("XMLParser::getAttribute not an XMLElement " ) ;

    const char* at = el->Attribute( name ) ;

    if( at == 0 ){

      std::stringstream str ;
      str  << "XMLParser::getAttribute missing attribute \"" << name
	   << "\" in element <" << el->Value() << "/> " ;
      throw ParseException( str.str() ) ;
    }

    return std::string( at )  ;

  }

  std::string getOptionalXMLAttribute(const  TiXmlNode* node , const std::string& name ,
				      const std::string& defaultValue) {

    std::string result( defaultValue ) ;

    try{

      result = getXMLAttribute( node, name ) ;
    }
    catch( ParseException& p){
    }

    return result ;
  }



  // helper function (declared in XMLHandler.h )

  std::string getChildElementValue(const TiXmlNode* node , const std::string& name ) {

    const TiXmlElement* el = node->ToElement() ;
    if( el == 0 )
      throw ParseException("XMLParser::getAttribute not an XMLElement " ) ;

    const TiXmlElement* cE = el->FirstChildElement( name );

    if( cE == 0 ) {

      std::stringstream str ;
      str  << "XMLParser::getChildElementValue missing element \"" << name
	   << "\" in element <" << el->Value() << "/> " ;

      throw ParseException( str.str() ) ;
    }

    return getXMLAttribute( cE , "value" )   ;

  }

}
