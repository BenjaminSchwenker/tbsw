#include "gearxml/GearParametersXML.h"

#include "gearxml/tinyxml.h"
#include "gearimpl/GearParametersImpl.h"
#include "gear/GearMgr.h"

#include <vector>
#include <algorithm>
#include <sstream>


namespace gear {


  /** helper method */
  int intVal(  const std::string& str ) {
    return atoi( str.c_str() )  ; 
  }
  
  /** helper method */
  double doubleVal( const std::string& str )  {
    return (double) atof( str.c_str() )  ; 
  }

  /** Helper class for XMLParser
   */
  class Tokenizer{
    
    std::vector< std::string >& _tokens ;
    char _del ;
    char _last ;
  public:
    
    Tokenizer( std::vector< std::string >& tokens, char del ) : _tokens(tokens) , _del(del), _last(del) {
    }
    
    
    void operator()(const char& c) { 
      
      if( c != _del  ) {
	
	if( _last == _del  ) {
	  _tokens.push_back("") ; 
	}
	_tokens.back() += c ;
	result() ;
      }
      _last = c ;
      
    } 
    
    ~Tokenizer(){
    }
    
    std::vector<std::string> & result()  { 
      
      return _tokens ; 
      
    }
  };
  



  TiXmlElement GearParametersXML::toXML( const GearParameters & parameters ) const {

    
    TiXmlElement det("detector") ;

    getXMLForParameters( &det , &parameters ) ;
    
    return det ;
    
  }
    
    
  GearParameters* GearParametersXML::fromXML( const TiXmlElement* xmlElement, GearMgr* gearMgr) const {
    

    GearParametersImpl* gearParams = new GearParametersImpl ;


    setParametersFromXML( xmlElement, gearParams  ) ;


    if( gearMgr != 0 ) {

      std::string name = getXMLAttribute( xmlElement , "name" ) ;

      gearMgr->setGearParameters( name , gearParams ) ;

    }


    return gearParams ;
  }

  void GearParametersXML::getXMLForParameters( TiXmlElement* xmlElement, 
					       const GearParameters* gearParams ){

    typedef const std::vector< std::string > KeyVec ;

    if( xmlElement == 0 || gearParams == 0 ){
      return ;
    }
    
    KeyVec& intKeys = gearParams->getIntKeys() ;
    for(unsigned int i=0 ; i < intKeys.size() ; ++i ){
      
      TiXmlElement param("parameter") ;
      
      param.SetAttribute( "name", intKeys[i] ) ;
      param.SetAttribute( "type", "int" ) ;
      param.SetAttribute( "value", gearParams->getIntVal( intKeys[i] )  ) ;
      
      xmlElement->InsertEndChild( param ) ;

    }

    KeyVec& doubleKeys = gearParams->getDoubleKeys() ;

    for(unsigned int i=0 ; i < doubleKeys.size() ; ++i ){
      
      TiXmlElement param("parameter") ;
      
      param.SetAttribute( "name", doubleKeys[i] ) ;
      param.SetAttribute( "type", "double" ) ;
      param.SetDoubleAttribute( "value", gearParams->getDoubleVal( doubleKeys[i] )  ) ;
      
      xmlElement->InsertEndChild( param ) ;

    }

    KeyVec& stringKeys = gearParams->getStringKeys() ;

    for(unsigned int i=0 ; i < stringKeys.size() ; ++i ){
      
      TiXmlElement param("parameter") ;
      
      param.SetAttribute( "name", stringKeys[i] ) ;
      param.SetAttribute( "type", "string" ) ;
      param.SetAttribute( "value", gearParams->getStringVal( stringKeys[i] )  ) ;
      
      xmlElement->InsertEndChild( param ) ;

    }

    KeyVec& intVecKeys = gearParams->getIntVecKeys() ;
    for(unsigned int i=0 ; i < intVecKeys.size() ; ++i ){
      
      TiXmlElement param("parameter") ;
      
      param.SetAttribute( "name", intVecKeys[i] ) ;
      param.SetAttribute( "type", "IntVec" ) ;

      IntVec vec = gearParams->getIntVals( intVecKeys[i]  ) ;
      std::stringstream str ;
      for(unsigned int j=0 ; j < vec.size() ; ++j ){
	if( j!=0 ) 
	  str << " " ; 
	str << vec[j] ;
      }      
      param.SetAttribute( "value", str.str() ) ;
      
      xmlElement->InsertEndChild( param ) ;
      
    }
    
    KeyVec& doubleVecKeys = gearParams->getDoubleVecKeys() ;
    for(unsigned int i=0 ; i < doubleVecKeys.size() ; ++i ){
      
      TiXmlElement param("parameter") ;
      
      param.SetAttribute( "name", doubleVecKeys[i] ) ;
      param.SetAttribute( "type", "DoubleVec" ) ;

      DoubleVec vec = gearParams->getDoubleVals( doubleVecKeys[i]  ) ;
      std::stringstream str ;
      for(unsigned int j=0 ; j < vec.size() ; ++j ){
	if( j!=0 ) 
	  str << " " ; 
	str << vec[j] ;
      }      
      param.SetAttribute( "value", str.str() ) ;
      
      xmlElement->InsertEndChild( param ) ;
      
    }

    KeyVec& stringVecKeys = gearParams->getStringVecKeys() ;
    for(unsigned int i=0 ; i < stringVecKeys.size() ; ++i ){
      
      TiXmlElement param("parameter") ;
      
      param.SetAttribute( "name", stringVecKeys[i] ) ;
      param.SetAttribute( "type", "StringVec" ) ;

      StringVec vec = gearParams->getStringVals( stringVecKeys[i]  ) ;
      std::stringstream str ;
      for(unsigned int j=0 ; j < vec.size() ; ++j ){
	if( j!=0 ) 
	  str << " " ; 
	str << vec[j] ;
      }      
      param.SetAttribute( "value", str.str() ) ;
      
      xmlElement->InsertEndChild( param ) ;
      
    }
    
  }

  void GearParametersXML::setParametersFromXML( const TiXmlElement* xmlElement, 
						GearParametersImpl* gearParams ){


    const TiXmlNode* par = 0 ;
    while( ( par = xmlElement->IterateChildren( "parameter", par ) )  != 0  ){
      
      
      
      // get the parameter value  from attribute or contents of element

      std::string value("") ;
      
      try{  
	
	value = getXMLAttribute( par , "value" )  ; 
      }      
      catch( ParseException ) {
	
	if( par->FirstChild() )
	  value =  par->FirstChild()->Value() ;
      }
      
      
      std::vector<std::string> stringValues ;

      Tokenizer t( stringValues ,' ') ;

      std::for_each( value.begin(), value.end(), t ) ; 


	  // fg: protect against empty value strings  in the gear file
	  if( stringValues.size() == 0 )

		  continue ;

 

      std::string name = getXMLAttribute( par, "name" ) ;
      std::string type = getXMLAttribute( par, "type" ) ;


      if( type == "int" ){

	gearParams->setIntVal( name ,  intVal( stringValues[0] )   )   ;
	
      } else if( type == "double" ){
	
	gearParams->setDoubleVal( name ,  doubleVal( stringValues[0] )   )   ;
	
      } else if( type == "string" ){
	
	gearParams->setStringVal( name ,  stringValues[0]    )   ;
	
      } else if( type == "IntVec" ){

	IntVec v ;
	transform(  stringValues.begin() ,  stringValues.end()  , back_inserter(v) , intVal ) ;
	gearParams->setIntVals( name , v ) ;
	

      } else if( type == "DoubleVec" ){
	
	DoubleVec v ;
	transform(  stringValues.begin() ,  stringValues.end()  , back_inserter(v) , doubleVal ) ;
	gearParams->setDoubleVals( name , v ) ;
	
      } else if( type == "StringVec" ){
	
	gearParams->setStringVals( name , stringValues  ) ;
	
      }
      
    }


  }



//   std::string GearParametersXML::getAttribute(const TiXmlNode* node , const std::string& name ) const {
    
//     const TiXmlElement* el = node->ToElement() ;
//     if( el == 0 ) 
//       throw ParseException("XMLParser::getAttribute not an XMLElement " ) ;
    
//     const char* at = el->Attribute( name ) ;
    
//     if( at == 0 ){
      
//       std::stringstream str ;
//       str  << "XMLParser::getAttribute missing attribute \"" << name 
// 	   << "\" in element <" << el->Value() << "/> " ;
//       throw ParseException( str.str() ) ;
//     }

//     return std::string( at )  ;

//   }


  
  
} // namespace
