#include "gearxml/SiPlanesParametersXML.h"

#include "gearxml/XMLHandlerMgr.h"
#include "gearxml/GearParametersXML.h"

#include "gearxml/tinyxml.h"
#include "gearimpl/SiPlanesParametersImpl.h"

#include "gear/GearMgr.h"

#include <vector>
#include <tuple>
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

    // Set up Beam Telescope as Element
    TiXmlElement det("detector") ;

    TiXmlElement setup_id( "siplanesID" ) ;
    setup_id.SetAttribute("ID", param->getSiPlanesID()) ;
    det.InsertEndChild( setup_id ) ;
        
    TiXmlElement nplanes( "siplanesNumber" ) ;
    nplanes.SetAttribute("number", param->getSiPlanesNumber()) ;
    det.InsertEndChild( nplanes ) ;

    // layerLayout
    const SiPlanesLayerLayout& siplanesLayers = param->getSiPlanesLayerLayout() ;
     
    TiXmlElement layers("layers") ;

    for( int i=0 ; i < siplanesLayers.getNLayers() ; i++ ) {

       //      std::cout << " layer " <<  i  << std::endl ; //debug
      
      TiXmlElement layer("layer" ) ;
       
      TiXmlElement sens("sensitive" ) ;
      sens.SetAttribute( "ID" , siplanesLayers.getSensitiveID( i ) ) ;
      sens.SetDoubleAttribute( "positionX" , siplanesLayers.getSensitivePositionX( i ) ) ;
      sens.SetDoubleAttribute( "positionY" , siplanesLayers.getSensitivePositionY( i ) ) ;
      sens.SetDoubleAttribute( "positionZ" , siplanesLayers.getSensitivePositionZ( i ) ) ;
      sens.SetDoubleAttribute( "thickness" , siplanesLayers.getSensitiveThickness( i ) ) ;
      sens.SetDoubleAttribute( "radLength" , siplanesLayers.getSensitiveRadLength( i ) ) ;
      sens.SetDoubleAttribute( "atomicNumber", siplanesLayers.getSensitiveAtomicNumber( i ) ) ;
      sens.SetDoubleAttribute( "atomicMass", siplanesLayers.getSensitiveAtomicMass( i ) ) ;
      sens.SetDoubleAttribute( "alpha" , siplanesLayers.getSensitiveRotationAlpha( i ) ) ;
      sens.SetDoubleAttribute( "beta" , siplanesLayers.getSensitiveRotationBeta( i ) ) ;
      sens.SetDoubleAttribute( "gamma" , siplanesLayers.getSensitiveRotationGamma( i ) ) ;
      sens.SetDoubleAttribute( "rotation1" , siplanesLayers.getSensitiveRotation1( i ) ) ;
      sens.SetDoubleAttribute( "rotation2" , siplanesLayers.getSensitiveRotation2( i ) ) ;
      sens.SetDoubleAttribute( "rotation3" , siplanesLayers.getSensitiveRotation3( i ) ) ;
      sens.SetDoubleAttribute( "rotation4" , siplanesLayers.getSensitiveRotation4( i ) ) ;
      
      // assemble layer
      layer.InsertEndChild(sens) ;
      
      TiXmlElement ladder("ladder") ;
      ladder.SetDoubleAttribute( "sizeU" , siplanesLayers.getLayerSizeU( i ) ) ;
      ladder.SetDoubleAttribute( "sizeV" , siplanesLayers.getLayerSizeV( i ) ) ;
      ladder.SetDoubleAttribute( "thickness" , siplanesLayers.getLayerThickness( i ) ) ;
      ladder.SetDoubleAttribute( "radLength" , siplanesLayers.getLayerRadLength( i ) ) ;
      ladder.SetDoubleAttribute( "atomicNumber", siplanesLayers.getLayerAtomicNumber( i ) ) ;
      ladder.SetDoubleAttribute( "atomicMass", siplanesLayers.getLayerAtomicMass( i ) ) ;
      
      // assemble layer
      layer.InsertEndChild(ladder) ;
      
      for (auto cellGroup : siplanesLayers.getSensitiveUCells( i ) ) {
        int sMinCell = std::get<0>(cellGroup);
        int sMaxCell = std::get<1>(cellGroup);
        double sPitch = std::get<2>(cellGroup);

        TiXmlElement group("uCellGroup") ;
        group.SetAttribute( "minCell" , sMinCell ) ;
        group.SetAttribute( "maxCell" , sMaxCell ) ;
        group.SetDoubleAttribute( "pitch", sPitch ) ;
        
        // assemble layer
        layer.InsertEndChild(group) ;
      }
  
      for (auto cellGroup : siplanesLayers.getSensitiveVCells( i ) ) {
        int sMinCell = std::get<0>(cellGroup);
        int sMaxCell = std::get<1>(cellGroup);
        double sPitch = std::get<2>(cellGroup);
        
        TiXmlElement group("vCellGroup") ;
        group.SetAttribute( "minCell" , sMinCell ) ;
        group.SetAttribute( "maxCell" , sMaxCell ) ;
        group.SetDoubleAttribute( "pitch", sPitch ) ;
        
        // assemble layer
        layer.InsertEndChild(group) ;
      }
            
      // assemble layers
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
     
    // number of telescope planes
    
    const TiXmlElement* siplanesNumber = xmlElement->FirstChildElement( "siplanesNumber" ) ;
    int nplanes = atoi( getXMLAttribute( siplanesNumber , "number" ).c_str() ) ;
    
    // create SiPlanesParameters
    SiPlanesParametersImpl* siplanesParam = new SiPlanesParametersImpl( setupID, nplanes) ;
    
    // layers
    const TiXmlNode* xmlLayers = xmlElement->FirstChildElement( "layers" ) ;

    const TiXmlNode* xmlLayer = 0 ;
    while( ( xmlLayer = xmlLayers->IterateChildren( "layer" , xmlLayer ) ) != 0 ) {
      
      const TiXmlNode* xmlSensitive = xmlLayer->FirstChildElement( "sensitive" ) ;
      const TiXmlNode* xmlLadder = xmlLayer->FirstChildElement( "ladder" ) ;
      const TiXmlNode* xmlUCellGroup = 0;  
      const TiXmlNode* xmlVCellGroup = 0; 
       
      int sID = atoi( getXMLAttribute( xmlSensitive , "ID" ).c_str() ) ;
      double sPosX   = atof( getXMLAttribute( xmlSensitive , "positionX" ).c_str() ) ;
      double sPosY   = atof( getXMLAttribute( xmlSensitive , "positionY" ).c_str() ) ;
      double sPosZ   = atof( getXMLAttribute( xmlSensitive , "positionZ" ).c_str() ) ;
      double sThick   = atof( getXMLAttribute( xmlSensitive , "thickness" ).c_str() ) ;
      double sRadLen = atof( getXMLAttribute( xmlSensitive , "radLength" ).c_str() ) ;  
      double sAtomicNum = atof( getXMLAttribute( xmlSensitive , "atomicNumber" ).c_str() ) ;  
      double sAtomicMass = atof( getXMLAttribute( xmlSensitive , "atomicMass" ).c_str() ) ;  
      double sAlpha  = atof(getXMLAttribute( xmlSensitive , "alpha" ).c_str() ) ;
      double sBeta   = atof(getXMLAttribute( xmlSensitive , "beta" ).c_str() ) ;
      double sGamma  = atof(getXMLAttribute( xmlSensitive , "gamma" ).c_str() ) ;
      double sRotat1 = atof(getXMLAttribute( xmlSensitive , "rotation1" ).c_str() ) ;
      double sRotat2 = atof(getXMLAttribute( xmlSensitive , "rotation2" ).c_str() ) ;
      double sRotat3 = atof(getXMLAttribute( xmlSensitive , "rotation3" ).c_str() ) ;
      double sRotat4 = atof(getXMLAttribute( xmlSensitive , "rotation4" ).c_str() ) ;
       
      double lSizU   = atof( getXMLAttribute( xmlLadder , "sizeU" ).c_str() ) ;
      double lSizV   = atof( getXMLAttribute( xmlLadder , "sizeV" ).c_str() ) ;
      double lThick   = atof( getXMLAttribute( xmlLadder , "thickness" ).c_str() ) ;
      double lRadLen = atof( getXMLAttribute( xmlLadder , "radLength" ).c_str() ) ;
      double lAtomicNum = atof( getXMLAttribute( xmlLadder , "atomicNumber" ).c_str() ) ;  
      double lAtomicMass  = atof(getXMLAttribute( xmlLadder , "atomicMass" ).c_str() ) ;
      
      std::vector< std::tuple<int,int,double> > uCellGroupVec; 
      while( ( xmlUCellGroup = xmlLayer->IterateChildren( "uCellGroup" , xmlUCellGroup ) ) != 0 ) {
        int sMinCell = atoi( getXMLAttribute( xmlUCellGroup , "minCell" ).c_str() ) ;
        int sMaxCell = atoi( getXMLAttribute( xmlUCellGroup , "maxCell" ).c_str() ) ;
        double sPitch = atof(getXMLAttribute( xmlUCellGroup , "pitch" ).c_str() ) ;
        uCellGroupVec.push_back( std::tuple<int, int, double>(sMinCell, sMaxCell, sPitch) );   
      }

      std::vector< std::tuple<int,int,double> > vCellGroupVec; 
      while( ( xmlVCellGroup = xmlLayer->IterateChildren( "vCellGroup" , xmlVCellGroup ) ) != 0 ) {
        int sMinCell = atoi( getXMLAttribute( xmlVCellGroup , "minCell" ).c_str() ) ;
        int sMaxCell = atoi( getXMLAttribute( xmlVCellGroup , "maxCell" ).c_str() ) ;
        double sPitch = atof(getXMLAttribute( xmlVCellGroup , "pitch" ).c_str() ) ;
        vCellGroupVec.push_back( std::tuple<int, int, double>(sMinCell, sMaxCell, sPitch) );
      }      
      
      siplanesParam->addLayer(sID, sPosX, sPosY, sPosZ, sThick, sRadLen, sAtomicNum, sAtomicMass, sAlpha, sBeta, sGamma, sRotat1, sRotat2, sRotat3, sRotat4, uCellGroupVec, vCellGroupVec, lSizU, lSizV, lThick, lRadLen, lAtomicNum, lAtomicMass) ;
     
    } // end loop


    // ------- add to GearMgr ----
    if( gearMgr != 0 ) {
      
      gearMgr->setSiPlanesParameters( siplanesParam ) ;

    }

    return siplanesParam ;

  } // fromXML

} // namespace
  
    
