#include "gearxml/SiPlanesParametersXML.h"

#include "gearxml/XMLHandlerMgr.h"
#include "gearxml/GearParametersXML.h"

#include "gearxml/tinyxml.h"
#include "gearimpl/SiPlanesParametersImpl.h"

#include "gear/GearMgr.h"

#include <vector>
#include <string>

 
#include <cstring>


namespace gear {

  TiXmlElement SiPlanesParametersXML::toXML( const GearParameters & parameters ) const {

    //   std::cout << "SiPlanesParameters::toXML called" << std::endl ; //debug
    
    // check whether parameter is valid SiPlanesParameter
    const SiPlanesParameters* param = dynamic_cast<const SiPlanesParameters*> ( &parameters ) ;

    if( param == 0 ) {

      throw Exception( "SiPlanesParametersXML::toXML given parameters not of correct type. "
		       "needs to be gear::SiPlanesParameters." ) ;
    }

    // Set up Beam TelescopeWithDUT or TelescopeWithoutDUT as Element
    TiXmlElement det("detector") ;

    TiXmlElement setup_id( "siplanesID" ) ;
    setup_id.SetAttribute("ID", param->getSiPlanesID()) ;
    det.InsertEndChild( setup_id ) ;
    
    //type
    TiXmlElement type("siplanesType");
    type.SetAttribute( "type",  "TelescopeWithoutDUT" ) ;
    
    det.InsertEndChild( type ) ;
    
    TiXmlElement nplanes( "siplanesNumber" ) ;
    nplanes.SetAttribute("number", param->getSiPlanesNumber()) ;
    det.InsertEndChild( nplanes ) ;

    // layerLayout
    const SiPlanesLayerLayout& siplanesLayers = param->getSiPlanesLayerLayout() ;
     
    TiXmlElement layers("layers") ;

    for( int i=0 ; i < siplanesLayers.getNLayers() ; i++ ) {

       //      std::cout << " layer " <<  i  << std::endl ; //debug
      
      TiXmlElement layer("layer" ) ;

      TiXmlElement ladder("ladder") ;
      ladder.SetAttribute( "ID" , siplanesLayers.getID( i ) ) ;
      ladder.SetDoubleAttribute( "positionX" , siplanesLayers.getLayerPositionX( i ) ) ;
      ladder.SetDoubleAttribute( "positionY" , siplanesLayers.getLayerPositionY( i ) ) ;
      ladder.SetDoubleAttribute( "positionZ" , siplanesLayers.getLayerPositionZ( i ) ) ;
      ladder.SetDoubleAttribute( "sizeX" , siplanesLayers.getLayerSizeX( i ) ) ;
      ladder.SetDoubleAttribute( "sizeY" , siplanesLayers.getLayerSizeY( i ) ) ;
      ladder.SetDoubleAttribute( "thickness" , siplanesLayers.getLayerThickness( i ) ) ;
      ladder.SetDoubleAttribute( "radLength" , siplanesLayers.getLayerRadLength( i ) ) ;
      
      TiXmlElement sens("sensitive" ) ;
      sens.SetAttribute( "ID" , siplanesLayers.getSensitiveID( i ) ) ;
      sens.SetAttribute( "PixType" , siplanesLayers.getSensitivePixelType( i ) ) ;
      sens.SetDoubleAttribute( "positionX" , siplanesLayers.getSensitivePositionX( i ) ) ;
      sens.SetDoubleAttribute( "positionY" , siplanesLayers.getSensitivePositionY( i ) ) ;
      sens.SetDoubleAttribute( "positionZ" , siplanesLayers.getSensitivePositionZ( i ) ) ;
      sens.SetDoubleAttribute( "sizeX" , siplanesLayers.getSensitiveSizeX( i ) ) ;
      sens.SetDoubleAttribute( "sizeY" , siplanesLayers.getSensitiveSizeY( i ) ) ;
      sens.SetDoubleAttribute( "thickness" , siplanesLayers.getSensitiveThickness( i ) ) ;
      sens.SetAttribute( "npixelX" , siplanesLayers.getSensitiveNpixelX( i ) ) ;
      sens.SetAttribute( "npixelY" , siplanesLayers.getSensitiveNpixelY( i ) ) ;
      sens.SetDoubleAttribute( "pitchX" , siplanesLayers.getSensitivePitchX( i ) ) ;
      sens.SetDoubleAttribute( "pitchY" , siplanesLayers.getSensitivePitchY( i ) ) ;
      sens.SetDoubleAttribute( "resolutionX" , siplanesLayers.getSensitiveResolutionX( i ) ) ;
      sens.SetDoubleAttribute( "resolutionY" , siplanesLayers.getSensitiveResolutionY( i ) ) ;
      sens.SetDoubleAttribute( "alpha" , siplanesLayers.getSensitiveRotationAlpha( i ) ) ;
      sens.SetDoubleAttribute( "beta" , siplanesLayers.getSensitiveRotationBeta( i ) ) ;
      sens.SetDoubleAttribute( "gamma" , siplanesLayers.getSensitiveRotationGamma( i ) ) ;
      sens.SetDoubleAttribute( "rotation1" , siplanesLayers.getSensitiveRotation1( i ) ) ;
      sens.SetDoubleAttribute( "rotation2" , siplanesLayers.getSensitiveRotation2( i ) ) ;
      sens.SetDoubleAttribute( "rotation3" , siplanesLayers.getSensitiveRotation3( i ) ) ;
      sens.SetDoubleAttribute( "rotation4" , siplanesLayers.getSensitiveRotation4( i ) ) ;
      sens.SetDoubleAttribute( "radLength" , siplanesLayers.getSensitiveRadLength( i ) ) ;

      // assemble layer
      layer.InsertEndChild(ladder) ;
      layer.InsertEndChild(sens) ;
      layers.InsertEndChild(layer);
      
    }
    
    det.InsertEndChild(layers) ;

    // Assemble Detector
    GearParametersXML::getXMLForParameters( &det , &parameters ) ;

    return det ;

  }

  GearParameters* SiPlanesParametersXML::fromXML( const TiXmlElement* xmlElement, GearMgr* gearMgr) const {

    // setup ID

    const TiXmlElement* siplanesID = xmlElement->FirstChildElement( "siplanesID" ) ;
    int setupID = atoi( getXMLAttribute( siplanesID , "ID" ).c_str() ) ;
    
    int intType = 0 ;
    
    // number of telescope planes
    
    const TiXmlElement* siplanesNumber = xmlElement->FirstChildElement( "siplanesNumber" ) ;
    int nplanes = atoi( getXMLAttribute( siplanesNumber , "number" ).c_str() ) ;
    
    //    std::cout << "SiPlanesParameters::fromXML siplanesNumber == " << nplanes << std::endl ; // debug
    
    // create SiPlanesParameters
    SiPlanesParametersImpl* siplanesParam = new SiPlanesParametersImpl( setupID, intType , nplanes) ;

    

    // layers
    const TiXmlNode* xmlLayers = xmlElement->FirstChildElement( "layers" ) ;

    const TiXmlNode* xmlLayer = 0 ;
    while( ( xmlLayer = xmlLayers->IterateChildren( "layer" , xmlLayer ) ) != 0 ) {

    const TiXmlNode* xmlLad = xmlLayer->FirstChildElement( "ladder" ) ;
    const TiXmlNode* xmlSen = xmlLayer->FirstChildElement( "sensitive" ) ;
    
    int lID = atoi( getXMLAttribute( xmlLad , "ID" ).c_str() ) ;
    double lPosX   = atof( getXMLAttribute( xmlLad , "positionX" ).c_str() ) ;
    double lPosY   = atof( getXMLAttribute( xmlLad , "positionY" ).c_str() ) ;
    double lPosZ   = atof( getXMLAttribute( xmlLad , "positionZ" ).c_str() ) ;
    double lSizX   = atof( getXMLAttribute( xmlLad , "sizeX" ).c_str() ) ;
    double lSizY   = atof( getXMLAttribute( xmlLad , "sizeY" ).c_str() ) ;
    double lThick   = atof( getXMLAttribute( xmlLad , "thickness" ).c_str() ) ;
    double lRadLen = atof( getXMLAttribute( xmlLad , "radLength" ).c_str() ) ;

    int sID = atoi( getXMLAttribute( xmlSen , "ID" ).c_str() ) ;
    int sPixType = atoi( getXMLAttribute( xmlSen , "PixType" ).c_str() ) ;
    double sPosX   = atof( getXMLAttribute( xmlSen , "positionX" ).c_str() ) ;
    double sPosY   = atof( getXMLAttribute( xmlSen , "positionY" ).c_str() ) ;
    double sPosZ   = atof( getXMLAttribute( xmlSen , "positionZ" ).c_str() ) ;
    double sSizX   = atof( getXMLAttribute( xmlSen , "sizeX" ).c_str() ) ;
    double sSizY   = atof( getXMLAttribute( xmlSen , "sizeY" ).c_str() ) ;
    double sThick   = atof( getXMLAttribute( xmlSen , "thickness" ).c_str() ) ;
    int sNPixX   = atoi(getXMLAttribute( xmlSen , "npixelX" ).c_str() ) ;
    int sNPixY   = atoi(getXMLAttribute( xmlSen , "npixelY" ).c_str() ) ;
    double sPitX   = atof(getXMLAttribute( xmlSen , "pitchX" ).c_str() ) ;
    double sPitY   = atof(getXMLAttribute( xmlSen , "pitchY" ).c_str() ) ;
    double sResolX   = atof(getXMLAttribute( xmlSen , "resolutionX" ).c_str() ) ;
    double sResolY   = atof(getXMLAttribute( xmlSen , "resolutionY" ).c_str() ) ;
    double sAlpha  = atof(getXMLAttribute( xmlSen , "alpha" ).c_str() ) ;
    double sBeta   = atof(getXMLAttribute( xmlSen , "beta" ).c_str() ) ;
    double sGamma  = atof(getXMLAttribute( xmlSen , "gamma" ).c_str() ) ;
    double sRotat1 = atof(getXMLAttribute( xmlSen , "rotation1" ).c_str() ) ;
    double sRotat2 = atof(getXMLAttribute( xmlSen , "rotation2" ).c_str() ) ;
    double sRotat3 = atof(getXMLAttribute( xmlSen , "rotation3" ).c_str() ) ;
    double sRotat4 = atof(getXMLAttribute( xmlSen , "rotation4" ).c_str() ) ;
    double sRadLen = atof( getXMLAttribute( xmlSen , "radLength" ).c_str() ) ;
    
    siplanesParam->addLayer(lID, lPosX, lPosY, lPosZ, lSizX, lSizY, lThick, lRadLen, sID, sPixType, sPosX, sPosY, sPosZ, sSizX, sSizY, sThick, sNPixX, sNPixY, sPitX, sPitY, sResolX, sResolY, sAlpha, sBeta, sGamma, sRotat1, sRotat2,sRotat3,sRotat4,sRadLen) ;

    } // end loop


    // ------- add to GearMgr ----
    if( gearMgr != 0 ) {
      
      gearMgr->setSiPlanesParameters( siplanesParam ) ;

    }

    return siplanesParam ;

  } // fromXML

} // namespace
  
    
