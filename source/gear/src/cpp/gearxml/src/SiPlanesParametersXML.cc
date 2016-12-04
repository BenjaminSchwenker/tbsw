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
    TiXmlElement type( "siplanesType" ) ;
    std::string strType ;
    switch( param->getSiPlanesType() ) {
    case ( SiPlanesParameters::TelescopeWithDUT ) :
      strType = "TelescopeWithDUT" ;
      break ;
    case ( SiPlanesParameters::TelescopeWithoutDUT ) :
      strType = "TelescopeWithoutDUT" ;
      break ;
    default :
      strType = "Unknown" ;
    }

    type.SetAttribute( "type", strType ) ;
    det.InsertEndChild( type ) ;

    //    std::cout << " setup ID = " << param->getSiPlanesID() << std::endl ; //debug

    //    std::cout << "SiPlanesParameters::toXML strType == '" << strType << "'"<< std::endl ; // debug

    //    std::cout << " N of planes = " << param->getSiPlanesNumber() << std::endl ; //debug


    TiXmlElement nplanes( "siplanesNumber" ) ;
    nplanes.SetAttribute("number", param->getSiPlanesNumber()) ;
    det.InsertEndChild( nplanes ) ;

    // layerLayout
    const SiPlanesLayerLayout& siplanesLayers = param->getSiPlanesLayerLayout() ;

    if (strType == "TelescopeWithDUT"){
      
      TiXmlElement dut("dut" ) ;
      
      TiXmlElement dutladder("dut_ladder") ;
      dutladder.SetAttribute( "ID" , siplanesLayers.getDUTID() ) ;
      dutladder.SetDoubleAttribute( "positionX" , siplanesLayers.getDUTPositionX() ) ;
      dutladder.SetDoubleAttribute( "positionY" , siplanesLayers.getDUTPositionY() ) ;
      dutladder.SetDoubleAttribute( "positionZ" , siplanesLayers.getDUTPositionZ() ) ;
      dutladder.SetDoubleAttribute( "sizeX" , siplanesLayers.getDUTSizeX() ) ;
      dutladder.SetDoubleAttribute( "sizeY" , siplanesLayers.getDUTSizeY() ) ;
      dutladder.SetDoubleAttribute( "thickness" , siplanesLayers.getDUTThickness() ) ;
      dutladder.SetDoubleAttribute( "radLength" , siplanesLayers.getDUTRadLength() ) ;
      
      TiXmlElement dutsens("dut_sensitive" ) ;
      dutsens.SetAttribute( "ID" , siplanesLayers.getDUTSensitiveID() ) ;
      dutsens.SetAttribute( "PixType" , siplanesLayers.getDUTSensitivePixelType() ) ;
      dutsens.SetDoubleAttribute( "positionX" , siplanesLayers.getDUTSensitivePositionX() ) ;
      dutsens.SetDoubleAttribute( "positionY" , siplanesLayers.getDUTSensitivePositionY() ) ;
      dutsens.SetDoubleAttribute( "positionZ" , siplanesLayers.getDUTSensitivePositionZ() ) ;
      dutsens.SetDoubleAttribute( "sizeX" , siplanesLayers.getDUTSensitiveSizeX() ) ;
      dutsens.SetDoubleAttribute( "sizeY" , siplanesLayers.getDUTSensitiveSizeY() ) ;
      dutsens.SetDoubleAttribute( "thickness" , siplanesLayers.getDUTSensitiveThickness() ) ;
      dutsens.SetAttribute( "npixelX" , siplanesLayers.getDUTSensitiveNpixelX() ) ;
      dutsens.SetAttribute( "npixelY" , siplanesLayers.getDUTSensitiveNpixelY() ) ;
      dutsens.SetDoubleAttribute( "pitchX" , siplanesLayers.getDUTSensitivePitchX() ) ;
      dutsens.SetDoubleAttribute( "pitchY" , siplanesLayers.getDUTSensitivePitchY() ) ;
      dutsens.SetDoubleAttribute( "resolutionX" , siplanesLayers.getDUTSensitiveResolutionX() ) ;
      dutsens.SetDoubleAttribute( "resolutionY" , siplanesLayers.getDUTSensitiveResolutionY() ) ;
      dutsens.SetDoubleAttribute( "alpha" , siplanesLayers.getDUTSensitiveRotationAlpha() ) ;
      dutsens.SetDoubleAttribute( "beta" , siplanesLayers.getDUTSensitiveRotationBeta() ) ;
      dutsens.SetDoubleAttribute( "gamma" , siplanesLayers.getDUTSensitiveRotationGamma() ) ;
      dutsens.SetDoubleAttribute( "rotation1" , siplanesLayers.getDUTSensitiveRotation1() ) ;
      dutsens.SetDoubleAttribute( "rotation2" , siplanesLayers.getDUTSensitiveRotation2() ) ;
      dutsens.SetDoubleAttribute( "rotation3" , siplanesLayers.getDUTSensitiveRotation3() ) ;
      dutsens.SetDoubleAttribute( "rotation4" , siplanesLayers.getDUTSensitiveRotation4() ) ;
      dutsens.SetDoubleAttribute( "radLength" , siplanesLayers.getDUTSensitiveRadLength() ) ;
      
      // assemble dut
      dut.InsertEndChild(dutladder) ;
      dut.InsertEndChild(dutsens) ;
      det.InsertEndChild(dut) ;
      
    }
    
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

    //    std::cout << "SiPlanesParameters::fromXML siplanesID == " << setupID << std::endl ; // debug

    // type

    const TiXmlElement* siplanesType = xmlElement->FirstChildElement( "siplanesType" ) ;
    int intType = 0 ;
    const char* strType = getXMLAttribute( siplanesType , "type" ) .c_str() ;

    std::cout << "SiPlanesParameters::fromXML siplanesType == '" << strType << "'" << std::endl ; // debug

    if( !std::strcmp( strType , "TelescopeWithDUT" ) ) {
      intType = SiPlanesParameters::TelescopeWithDUT ;
     }
    if( !std::strcmp( strType ,"TelescopeWithoutDUT" ) ) {
      intType = SiPlanesParameters::TelescopeWithoutDUT ;
    }
    if( std::strcmp( strType , "TelescopeWithDUT" ) && std::strcmp( strType ,"TelescopeWithoutDUT" ) ) {
      std::cout << "strType " <<  strType << std::endl ;
      throw Exception( "SiPlanesParametersXML::fromXML type of SiPlanes not known."
		       "Needs to be 'TelescopeWithDUT' or 'TelescopeWithoutDUT'." ) ;
    }

    // number of telescope planes

    const TiXmlElement* siplanesNumber = xmlElement->FirstChildElement( "siplanesNumber" ) ;
    int nplanes = atoi( getXMLAttribute( siplanesNumber , "number" ).c_str() ) ;

    //    std::cout << "SiPlanesParameters::fromXML siplanesNumber == " << nplanes << std::endl ; // debug

    // create SiPlanesParameters
    SiPlanesParametersImpl* siplanesParam = new SiPlanesParametersImpl( setupID, intType , nplanes) ;

    if (intType == SiPlanesParameters::TelescopeWithDUT) {
      const TiXmlNode* xmlDUT = xmlElement->FirstChildElement( "dut" ) ;
      const TiXmlNode* xmlDUTLad = xmlDUT->FirstChildElement( "dut_ladder" ) ;
      const TiXmlNode* xmlDUTSen = xmlDUT->FirstChildElement( "dut_sensitive" ) ;
      
      int ldutID = atoi(getXMLAttribute( xmlDUTLad , "ID" ).c_str() );
      double ldutPosX   = atof(getXMLAttribute( xmlDUTLad , "positionX" ).c_str() ) ;
      double ldutPosY   = atof(getXMLAttribute( xmlDUTLad , "positionY" ).c_str() ) ;
      double ldutPosZ   = atof(getXMLAttribute( xmlDUTLad , "positionZ" ).c_str() ) ;
      double ldutSizX   = atof(getXMLAttribute( xmlDUTLad , "sizeX" ).c_str() ) ;
      double ldutSizY   = atof(getXMLAttribute( xmlDUTLad , "sizeY" ).c_str() ) ;
      double ldutThick   = atof(getXMLAttribute( xmlDUTLad , "thickness" ).c_str() ) ;
      double ldutRadLen = atof(getXMLAttribute( xmlDUTLad , "radLength" ).c_str() ) ;
      
      int sdutID = atoi(getXMLAttribute( xmlDUTSen , "ID" ).c_str() );
      int sdutPixType = atoi(getXMLAttribute( xmlDUTSen , "PixType" ).c_str() );
      double sdutPosX   = atof(getXMLAttribute( xmlDUTSen , "positionX" ).c_str() ) ;
      double sdutPosY   = atof(getXMLAttribute( xmlDUTSen , "positionY" ).c_str() ) ;
      double sdutPosZ   = atof(getXMLAttribute( xmlDUTSen , "positionZ" ).c_str() ) ;
      double sdutSizX   = atof(getXMLAttribute( xmlDUTSen , "sizeX" ).c_str() ) ;
      double sdutSizY   = atof(getXMLAttribute( xmlDUTSen , "sizeY" ).c_str() ) ;
      double sdutThick   = atof(getXMLAttribute( xmlDUTSen , "thickness" ).c_str() ) ;
      int sdutNPixX   = atoi(getXMLAttribute( xmlDUTSen , "npixelX" ).c_str() ) ;
      int sdutNPixY   = atoi(getXMLAttribute( xmlDUTSen , "npixelY" ).c_str() ) ;
      double sdutPitX   = atof(getXMLAttribute( xmlDUTSen , "pitchX" ).c_str() ) ;
      double sdutPitY   = atof(getXMLAttribute( xmlDUTSen , "pitchY" ).c_str() ) ;
      double sdutResolX   = atof(getXMLAttribute( xmlDUTSen , "resolutionX" ).c_str() ) ;
      double sdutResolY   = atof(getXMLAttribute( xmlDUTSen , "resolutionY" ).c_str() ) ;
      double sdutAlpha  =  atof(getXMLAttribute( xmlDUTSen , "alpha" ).c_str() ) ;
      double sdutBeta   =  atof(getXMLAttribute( xmlDUTSen , "beta" ).c_str() ) ;
      double sdutGamma  =  atof(getXMLAttribute( xmlDUTSen , "gamma" ).c_str() ) ;
      double sdutRotat1 =  atof(getXMLAttribute( xmlDUTSen , "rotation1" ).c_str() ) ;
      double sdutRotat2 =  atof(getXMLAttribute( xmlDUTSen , "rotation2" ).c_str() ) ;
      double sdutRotat3 =  atof(getXMLAttribute( xmlDUTSen , "rotation3" ).c_str() ) ;
      double sdutRotat4 =  atof(getXMLAttribute( xmlDUTSen , "rotation4" ).c_str() ) ;
      double sdutRadLen = atof(getXMLAttribute( xmlDUTSen , "radLength" ).c_str() ) ;
    
      siplanesParam->addDUT(ldutID, ldutPosX, ldutPosY, ldutPosZ, ldutSizX, ldutSizY, ldutThick, ldutRadLen, sdutID, sdutPixType, sdutPosX, sdutPosY, sdutPosZ, sdutSizX, sdutSizY, sdutThick, sdutNPixX, sdutNPixY, sdutPitX, sdutPitY, sdutResolX, sdutResolY, sdutAlpha, sdutBeta, sdutGamma, sdutRotat1, sdutRotat2,sdutRotat3,sdutRotat4, sdutRadLen) ;
    }

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
  
    
