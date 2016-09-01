#include "gearxml/ConstantBFieldXML.h"

#include "gearxml/XMLHandlerMgr.h"
#include "gearxml/GearParametersXML.h"

#include "gearxml/tinyxml.h"
#include "gearimpl/ConstantBField.h"
// #include "gear/LayerLayout.h"
#include "gear/GearMgr.h"


namespace gear {


  TiXmlElement ConstantBFieldXML::toXML( const GearParameters & parameters ) const {

    const ConstantBField* param = dynamic_cast<const ConstantBField*> ( &parameters );

    if( param == 0 ) {

      throw Exception( "ConstantBFieldXML::toXML given parmeters not of correct type. "
		       "needs to be gear::ConstantBField." );
    }

    TiXmlElement field("BField");

    field.SetAttribute("type","ConstantBField") ;

    Vector3D b = param->at( Vector3D( 0.,0.,0. ) ) ;

    field.SetDoubleAttribute("x",  b.x()  ) ;
    field.SetDoubleAttribute("y",  b.y()  ) ;
    field.SetDoubleAttribute("z",  b.z()  ) ;

    // Write all other parameters 
    GearParametersXML::getXMLForParameters( &field ,  param ) ;

    return field ;
  }  
  
  GearParameters* ConstantBFieldXML::fromXML( const TiXmlElement* xmlElement, GearMgr* gearMgr) const {
    
    double bx  =  atof(  getXMLAttribute( xmlElement , "x" ) .c_str() ) ;
    double by  =  atof(  getXMLAttribute( xmlElement , "y" ) .c_str() ) ;
    double bz  =  atof(  getXMLAttribute( xmlElement , "z" ) .c_str() ) ;

    ConstantBField* bfield = new ConstantBField(   Vector3D( bx, by, bz)  )  ;
    
    // now we can read the generic parameters
    GearParametersXML::setParametersFromXML( xmlElement, bfield  ) ;    


    //--------- add to proper section in GearMgr ----
    if( gearMgr != 0 ) {

      gearMgr->setBField( bfield ) ;	

    }
      
    return bfield ;
  }



} // namespace
