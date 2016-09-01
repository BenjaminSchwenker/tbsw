#include "gearxml/VXDParametersXML.h"

#include "gearxml/XMLHandlerMgr.h"
#include "gearxml/GearParametersXML.h"

#include "gearxml/tinyxml.h"
#include "gearimpl/VXDParametersImpl.h"

#include "gear/GearMgr.h"

#include <vector>
#include <string>

namespace gear {

  TiXmlElement VXDParametersXML::toXML( const GearParameters & parameters ) const {

    //std::cout << "VXDParameters::toXML called" << std::endl ; //debug

    // check wheterh parameter is valid VXDParameter
    const VXDParameters* param = dynamic_cast<const VXDParameters*> ( &parameters ) ;

    if( param == 0 ) {

      throw Exception( "VXDParametersXML::toXML given parameters not of correct type. "
		       "needs to be gear::VXDParameters." ) ;
    }

    // Set up VXD-Detector as Element
    TiXmlElement det("detector") ;

    //type
    TiXmlElement type( "vxdType" ) ;
    std::string strType ;
    switch( param->getVXDType() ) {
    case ( VXDParameters::CCD ) :
      strType = "CCD" ;
      break ;
    case ( VXDParameters::CMOS ) :
      strType = "CMOS" ;
      break ;
    case ( VXDParameters::HYBRID ) :
      strType = "HYBRID" ;
      break ;
    default :
      strType = "Unknown" ;
    }

    type.SetAttribute( "technology", strType ) ;
    det.InsertEndChild( type ) ;

    //std::cout << "VXDParameters::toXML strType == '" << strType << "'"<< std::endl ; // debug

    // shell
    // std::cout << "VXDParameters::toXML shell      " << std::endl ; // debug
//     std::cout << "                     halflength " << param->getShellHalfLength() << std::endl ; //debug
//     std::cout << "                            gap " << param ->getShellGap() << std::endl ; // debug
//     std::cout << "                         innerR " << param->getShellInnerRadius() << std::endl ; //debug
//     std::cout << "                         outerR " << param->getShellOuterRadius() << std::endl ; //debug
//     std::cout << "                         radLen " << param->getShellRadLength() << std::endl ; //debug
    
    TiXmlElement shell( "shell" ) ;
    shell.SetDoubleAttribute( "halfLength" , param->getShellHalfLength() ) ;
    shell.SetDoubleAttribute( "gap" , param->getShellGap() ) ;
    shell.SetDoubleAttribute( "innerRadius" , param->getShellInnerRadius() ) ;
    shell.SetDoubleAttribute( "outerRadius" , param->getShellOuterRadius() ) ;
    shell.SetDoubleAttribute( "radLength" , param->getShellRadLength() ) ;
    det.InsertEndChild( shell ) ;
    
    // layerLayout
    const VXDLayerLayout& vxdLayers = param->getVXDLayerLayout() ;
    TiXmlElement layers("layers") ;

    for( int i=0 ; i < vxdLayers.getNLayers() ; i++ ) {

      //std::cout << "                     layer #ladders " << vxdLayers.getNLadders( i ) << std::endl ; //debug
      
      TiXmlElement layer("layer" ) ;
      layer.SetAttribute( "nLadders" , vxdLayers.getNLadders( i ) ) ;
      layer.SetDoubleAttribute( "phi0" , vxdLayers.getPhi0( i ) ) ;
      
      TiXmlElement ladder("ladder") ;
      ladder.SetDoubleAttribute( "distance" , vxdLayers.getLadderDistance( i ) ) ;
      ladder.SetDoubleAttribute( "thickness" , vxdLayers.getLadderThickness( i ) ) ;
      ladder.SetDoubleAttribute( "width" , vxdLayers.getLadderWidth( i ) ) ;
      ladder.SetDoubleAttribute( "length" , vxdLayers.getLadderLength( i ) ) ;
      ladder.SetDoubleAttribute( "offset" , vxdLayers.getLadderOffset( i ) ) ;
      ladder.SetDoubleAttribute( "radLength" , vxdLayers.getLadderRadLength( i ) ) ;
      
      TiXmlElement sens("sensitive" ) ;
      sens.SetDoubleAttribute( "distance" , vxdLayers.getSensitiveDistance( i ) ) ;
      sens.SetDoubleAttribute( "thickness" , vxdLayers.getSensitiveThickness( i ) ) ;
      sens.SetDoubleAttribute( "width" , vxdLayers.getSensitiveWidth( i ) ) ;
      sens.SetDoubleAttribute( "length" , vxdLayers.getSensitiveLength( i ) ) ;
      sens.SetDoubleAttribute( "offset" , vxdLayers.getSensitiveOffset( i ) ) ;
      sens.SetDoubleAttribute( "radLength" , vxdLayers.getSensitiveRadLength( i ) ) ;

      // assemble layer
      layer.InsertEndChild(ladder);
      layer.InsertEndChild(sens) ;
      layers.InsertEndChild(layer) ;
      
    }

    // assemble layers
    det.InsertEndChild(layers) ;

    // Assemble Detector
    GearParametersXML::getXMLForParameters( &det , &parameters ) ;

    return det ;

  }

  GearParameters* VXDParametersXML::fromXML( const TiXmlElement* xmlElement, GearMgr* gearMgr) const {

    // type

    const TiXmlElement* vxdType = xmlElement->FirstChildElement( "vxdType" ) ;
    int intType = 0 ;
    
	std::string type =  getXMLAttribute( vxdType , "technology" )  ; 
	//const char* strType = type.c_str() ;

    //std::cout << "VXDParameters::fromXML vxdType == '" << strType << "'" 
	//		<< " string: " << type << std::endl ; // debug

    if(  type == "CCD"  ) {
      intType = VXDParameters::CCD ;
     }
	else if( type =="CMOS"  ) {
      intType = VXDParameters::CMOS ;
    }
	else if( type == "HYBRID"  ) {
      intType = VXDParameters::HYBRID ;
    }
	else{
	  throw Exception( "VXDParametersXML::fromXML technology of vxdType not known: " + type + 
	                   " - Needs to be 'CCD', 'CMOS' or 'HYBRID'." ) ;
    }

    const TiXmlElement* shell = xmlElement->FirstChildElement( "shell" ) ;

    double shellHalfLength = atof( getXMLAttribute( shell , "halfLength" ) .c_str() ) ;
    double shellGap = atof( getXMLAttribute( shell , "gap" ) . c_str() ) ;
    double shellInnerR = atof( getXMLAttribute( shell , "innerRadius" ) .c_str() ) ;
    double shellOuterR = atof( getXMLAttribute( shell , "outerRadius" ) .c_str() ) ;
    double shellRadLen = atof( getXMLAttribute( shell , "radLength" ) .c_str() ) ;

    // create VDXParameters
    VXDParametersImpl* vxdParam = new VXDParametersImpl( intType , shellInnerR , shellOuterR , shellHalfLength , shellGap , shellRadLen ) ;

    // layers
    const TiXmlNode* xmlLayers = xmlElement->FirstChildElement( "layers" ) ;
    const TiXmlNode* xmlLayer = 0 ;
    
    while( ( xmlLayer = xmlLayers->IterateChildren( "layer" , xmlLayer ) ) != 0 ) {
      
      int nLadders = atoi( getXMLAttribute( xmlLayer , "nLadders" ).c_str() ) ;
      double phi0  = atof( getXMLAttribute( xmlLayer , "phi0"     ).c_str() ) ;

      const TiXmlNode* xmlLad = xmlLayer->FirstChildElement( "ladder" ) ;
      const TiXmlNode* xmlSen = xmlLayer->FirstChildElement( "sensitive" ) ;

      double lDist   = atof( getXMLAttribute( xmlLad , "distance"  ).c_str() ) ;
      double lThick  = atof( getXMLAttribute( xmlLad , "thickness" ).c_str() ) ;
      double lWidth  = atof( getXMLAttribute( xmlLad , "width"     ).c_str() ) ;
      double lLength = atof( getXMLAttribute( xmlLad , "length"    ).c_str() ) ;
      double lOffset = atof( getXMLAttribute( xmlLad , "offset"    ).c_str() ) ;
      double lRadLen = atof( getXMLAttribute( xmlLad , "radLength" ).c_str() ) ;
    
      double sDist   = atof( getXMLAttribute( xmlSen , "distance"  ).c_str() ) ;
      double sThick  = atof( getXMLAttribute( xmlSen , "thickness" ).c_str() ) ;
      double sWidth  = atof( getXMLAttribute( xmlSen , "width"     ).c_str() ) ;
      double sLength = atof( getXMLAttribute( xmlSen , "length"    ).c_str() ) ;
      double sOffset = atof( getXMLAttribute( xmlSen , "offset"    ).c_str() ) ;
      double sRadLen = atof( getXMLAttribute( xmlSen , "radLength" ).c_str() ) ;

      vxdParam->addLayer( nLadders, phi0 ,
			  lDist, lOffset, lThick, lLength, lWidth, lRadLen,
			  sDist, sOffset, sThick, sLength, sWidth, sRadLen ) ;

    } // end loop

    // ------- add to GearMgr ----
    if( gearMgr != 0 ) {
      
      gearMgr->setVXDParameters( vxdParam ) ;

    }

    return vxdParam ;

  } // fromXML

} // namespace
  
    
