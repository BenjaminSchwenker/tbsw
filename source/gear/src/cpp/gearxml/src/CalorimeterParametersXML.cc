#include "gearxml/CalorimeterParametersXML.h"

#include "gearxml/XMLHandlerMgr.h"
#include "gearxml/GearParametersXML.h"

#include "gearxml/tinyxml.h"
#include "gearimpl/CalorimeterParametersImpl.h"
// #include "gear/LayerLayout.h"
#include "gear/GearMgr.h"

#include "gearimpl/FixedPadSizeDiskLayout.h"

#include <vector>
#include <algorithm>
#include <sstream>
#include <math.h>
#include <string>

#define _EPSILON 0.0001


namespace gear {


  TiXmlElement CalorimeterParametersXML::toXML( const GearParameters & parameters ) const {

    //check wheter parameter is valid CalorimeterParameter
    const CalorimeterParameters* param = dynamic_cast<const CalorimeterParameters*> ( &parameters );

    if( param == 0 ) {

      throw Exception( "CaloriemterParametersXML::toXML given parmeters not of correct type. "
		       "needs to be gear::CaloriemeterParameters." );
    }

    //Set up Caloriemter-Detector as Element
    TiXmlElement det("detector");

    //Dimension
    const std::vector<double> vecDim=param->getExtent();
    TiXmlElement dimensions("dimensions");

    //Layout
    TiXmlElement layout("layout");

    switch( param->getLayoutType() ) {
      
    case CalorimeterParameters::BARREL : 
      
      layout.SetAttribute("type","Barrel");
      dimensions.SetDoubleAttribute("inner_r",vecDim[0]);
      dimensions.SetDoubleAttribute("outer_z",vecDim[3]);
      break;
      
    case CalorimeterParameters::ENDCAP : 

      layout.SetAttribute("type","Endcap");
      dimensions.SetDoubleAttribute("inner_r",vecDim[0]);
      dimensions.SetDoubleAttribute("outer_r",vecDim[1]);
      dimensions.SetDoubleAttribute("inner_z",vecDim[2]);
      break;
      
    default: 

      layout.SetAttribute("type","Unknown");
      dimensions.SetDoubleAttribute("inner_r",vecDim[0]);
      dimensions.SetDoubleAttribute("outer_r",vecDim[1]);
      dimensions.SetDoubleAttribute("inner_z",vecDim[2]);
      dimensions.SetDoubleAttribute("outer_z",vecDim[3]);
      break;
    }

    layout.SetAttribute("symmetry",param->getSymmetryOrder());
    layout.SetDoubleAttribute("phi0",param->getPhi0());


    //Assemble first parts
    det.InsertEndChild(layout);
    det.InsertEndChild(dimensions);

    //LayerLayout
    const LayerLayout& layerLayout = param->getLayerLayout();
    
    //<------Go through layers to find out how to repeat
    // prepare currentValues for each repeater
    double curThickness = 0.0 ;
    double curAbsorber = 0.0 ;
    double curCell0 = 0.0 ;
    double curCell1 = 0.0 ;

    int repN = 0 ;
    
    for(int i=0; i<layerLayout.getNLayers();i++) {

      // check if anything has changed

      bool layerChanged = 	
	( ! isEqual( curThickness, layerLayout.getThickness(i) )         ) ||
	( ! isEqual( curAbsorber,  layerLayout.getAbsorberThickness(i) ) ) ||
	( ! isEqual( curCell0,     layerLayout.getCellSize0(i) )         ) ||
	( ! isEqual( curCell1,     layerLayout.getCellSize1(i) )         ) ;
      
      if( layerChanged && repN > 0 ) {

	//Write parameters

	TiXmlElement layer("layer");

	layer.SetAttribute("repeat",repN);

	layer.SetDoubleAttribute("thickness",curThickness);
	layer.SetDoubleAttribute("absorberThickness",curAbsorber);
	layer.SetDoubleAttribute("cellSize0",curCell0);
	layer.SetDoubleAttribute("cellSize1",curCell1);

	det.InsertEndChild(layer);

	repN = 1;

      } else {

	repN += 1;
      }

      // take on new values
      curThickness = layerLayout.getThickness(i);
      curAbsorber = layerLayout.getAbsorberThickness(i);
      curCell0 = layerLayout.getCellSize0(i);
      curCell1 = layerLayout.getCellSize1(i);	  
    }



    //Assemble last Repeater
    if ( repN > 0 ) {

      //Write parameters
      TiXmlElement layer("layer");
      layer.SetAttribute("repeat",repN);
      layer.SetDoubleAttribute("thickness",curThickness);
      layer.SetDoubleAttribute("absorberThickness",curAbsorber);
      layer.SetDoubleAttribute("cellSize0",curCell0);
      layer.SetDoubleAttribute("cellSize1",curCell1);
      det.InsertEndChild(layer);
    }
    //Assemble Detector
    
    GearParametersXML::getXMLForParameters( &det , &parameters ) ;

    return det ;
  }
  
  
  GearParameters* CalorimeterParametersXML::fromXML( const TiXmlElement* xmlElement, GearMgr* gearMgr) const {
    
    
    // ---- layout -----

    const TiXmlElement* layout = xmlElement->FirstChildElement( "layout" ) ;
 
    if( layout == 0 ) {
      
      std::stringstream str ;
      str  << "XMLParser::getChildElementValue missing element \"layout\" " 
	   << " in element <" << xmlElement->Value() << "/> " ;
      
      throw ParseException( str.str() ) ;
    }

    std::string layoutType = getXMLAttribute( layout , "type" )   ;

    std::transform( layoutType.begin() , layoutType.end() , layoutType.begin(),  tolower ) ;
  
    if( layoutType != "barrel" &&  layoutType != "endcap" ) {

      throw ParseException( "CalorimeterParametersXML::toXML: currently only barrel and endcap type layouts implemented ! " ) ;
    }


    // ---- dimensions -----
    const TiXmlElement* dimensions = xmlElement->FirstChildElement( "dimensions" ) ;
 
    if( dimensions == 0 ) {
      
      std::stringstream str ;
      str  << "XMLParser::getChildElementValue missing element \"dimensions\" " 
	   << " in element <" << xmlElement->Value() << "/> " ;
      
      throw ParseException( str.str() ) ;
    }

    int    symmetry   =  atoi(  getXMLAttribute( layout , "symmetry" ) .c_str() ) ;
    double phi0       =  atof(  getXMLAttribute( layout , "phi0" ) .c_str() ) ;

    CalorimeterParametersImpl* caloParams = 0 ;
    
    double firstLayerDistance = 0. ;
    bool isBarrelCalorimeter = false ;

    if( layoutType == "barrel" ) {

      isBarrelCalorimeter = true ;

      // barrel dimensions defined through inner_r and outer_z
      double rMin = atof(  getXMLAttribute( dimensions , "inner_r" ) .c_str() ) ;
      double zMax = atof(  getXMLAttribute( dimensions , "outer_z" ) .c_str() ) ;
      
      caloParams = new CalorimeterParametersImpl( rMin, zMax,  symmetry , phi0 ) ;

      firstLayerDistance = rMin ;

    } else if( layoutType == "endcap" ) { 
      
      // endcap dimensions defined through inner_r, outer_r  and inner_z
      double rMin = atof(  getXMLAttribute( dimensions , "inner_r" ) .c_str() ) ;
      double rMax = atof(  getXMLAttribute( dimensions , "outer_r" ) .c_str() ) ;
      double zMin = atof(  getXMLAttribute( dimensions , "inner_z" ) .c_str() ) ;
      
      caloParams = new CalorimeterParametersImpl( rMin, rMax, zMin,  symmetry , phi0 ) ;

      firstLayerDistance = zMin ;

    }

    // now we can read the generic parameters
    GearParametersXML::setParametersFromXML( xmlElement, caloParams  ) ;    


    // ----- layers ------------
    int layerId = 0 ;

    const TiXmlNode* layer = 0 ;
    while( ( layer = xmlElement->IterateChildren( "layer", layer ) )  != 0  ){

      int    nLayer            = atoi(  getXMLAttribute( layer , "repeat" ) .c_str() ) ;
      double thickness         = atof(  getXMLAttribute( layer , "thickness" ) .c_str() ) ;
      double absorberThickness = atof(  getXMLAttribute( layer , "absorberThickness" ) .c_str() ) ;
      double cellSize0         = atof(  getXMLAttribute( layer , "cellSize0" ) .c_str() ) ;
      double cellSize1         = atof(  getXMLAttribute( layer , "cellSize1" ) .c_str() ) ;
      
      for(int i=0 ; i < nLayer ; ++i ) {

	if( layerId == 0 )
	  caloParams->layerLayout().positionLayer( firstLayerDistance , thickness, 
						   cellSize0, cellSize1, absorberThickness ) ;
	else
	  caloParams->layerLayout().addLayer( thickness, 
					      cellSize0, cellSize1, absorberThickness ) ;

	++ layerId ;
      }
    }



    //--------- add to proper section in GearMgr ----
    if( gearMgr != 0 ) {

      std::string caloName = getXMLAttribute( xmlElement , "name" ) ;

//       bool isEcal = ( caloName[0] == 'E' || caloName[0] == 'e' ) ; 

      std::transform( caloName.begin() , caloName.end() , 
		      caloName.begin(), tolower ) ;

      
      if( caloName.find("ecal") != std::string::npos ){
	
	if( caloName.find("plug" ) != std::string::npos ){
	  
	  gearMgr->setEcalPlugParameters( caloParams ) ;	  
	}
	else if( isBarrelCalorimeter )
	  gearMgr->setEcalBarrelParameters( caloParams ) ;
	else
	  gearMgr->setEcalEndcapParameters( caloParams ) ;
	
      } else if( caloName.find("yoke") != std::string::npos ){
	
	if( caloName.find("plug" ) != std::string::npos ){
	  
	  gearMgr->setYokePlugParameters( caloParams ) ;	  
	}
	else if( isBarrelCalorimeter )
	  gearMgr->setYokeBarrelParameters( caloParams ) ;
	else
	  gearMgr->setYokeEndcapParameters( caloParams ) ;
	
      }  else if( caloName.find("lcal" ) != std::string::npos ) {
	
	gearMgr->setLcalParameters( caloParams ) ;	
	
      } else if( caloName.find("lhcal" ) != std::string::npos ) {
	
        gearMgr->setLHcalParameters( caloParams ) ;	
	
      } else if( caloName.find("beamcal" ) != std::string::npos ) {
	
	gearMgr->setBeamCalParameters( caloParams ) ;	

      } else if( caloName.find("hcal" ) != std::string::npos ) {
	
	// the additional hcal ring might be endcap or barrel like..
	if( caloName.find("ring" ) != std::string::npos ){

	  gearMgr->setHcalRingParameters( caloParams ) ;	  
	}
	else if( isBarrelCalorimeter )
	  gearMgr->setHcalBarrelParameters( caloParams ) ;
	else
	  gearMgr->setHcalEndcapParameters( caloParams ) ;
	
      }
    }
    
    return caloParams ;
  }


  bool CalorimeterParametersXML::isEqual( const double valueOne, const double valueTwo ) const {

     // save calculating time if equal
    if ( valueOne == valueTwo ) return true ;
    
    // get avg value and calculate maximum allowed delta
    double maxDelta = fabs( ( valueOne + valueTwo ) / 2 * _EPSILON ) ;
    
    // check if values differ by less than maximal delta
    bool differsLess = ( fabs( valueOne - valueTwo ) < maxDelta ) ;

    // return 
    return differsLess ;
  }

} // namespace
