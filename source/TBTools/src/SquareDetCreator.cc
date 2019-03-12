#include "SquareDetCreator.h"
#include "DetCreatorFactory.h"

#include "ThreeDModel.h"
#include "SquareDet.h"
#include "TBDetector.h"

#include <streamlog/streamlog.h>
#include <marlin/tinyxml.h>
#include <marlin/Exceptions.h>
#include <Eigen/LU>

using namespace depfet;
using namespace std;


namespace {
  std::string getXMLAttribute(const TiXmlNode* node , const std::string& name ) {
    
    const TiXmlElement* el = node->ToElement() ;
    if( el == 0 )
      throw marlin::ParseException("XMLParser::getAttribute not an XMLElement " ) ;
    
    const std::string*  at = el->Attribute( name ) ;
    
    if( at == 0 ){
      std::stringstream str ;
      str  << "XMLParser::getAttribute missing attribute \"" << name
	       << "\" in element <" << el->Value() << "/> " ;
      throw marlin::ParseException( str.str() ) ;
    }
     
    return std::string( *at )  ;
  }
} 

/** Register the creator */
// FIXME Change name from SiPlanesParameters to SquareDet, but this requires to changes xml files too. do it later.
DetCreatorFactory<SquareDetCreator> SquareDetFactory("SiPlanesParameters");

void SquareDetCreator::create(const TiXmlNode* content, std::vector<Det*>& dets) {
     
  // layers
  const TiXmlNode* xmlLayers = content->FirstChildElement( "layers" ) ;
      
  const TiXmlNode* xmlLayer = 0 ;  
  while( ( xmlLayer = xmlLayers->IterateChildren( "layer" , xmlLayer ) ) != 0 ) {
        
    const TiXmlNode* xmlSensitive = xmlLayer->FirstChildElement( "sensitive" ) ;
    const TiXmlNode* xmlLadder = xmlLayer->FirstChildElement( "ladder" ) ;
    const TiXmlNode* xmlUCellGroup = 0;  
    const TiXmlNode* xmlVCellGroup = 0; 
           
    int sensorID = atoi( getXMLAttribute( xmlSensitive , "ID" ).c_str() ) ; 
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
    int sRotat1 = atoi(getXMLAttribute( xmlSensitive , "rotation1" ).c_str() ) ;
    int sRotat2 = atoi(getXMLAttribute( xmlSensitive , "rotation2" ).c_str() ) ;
    int sRotat3 = atoi(getXMLAttribute( xmlSensitive , "rotation3" ).c_str() ) ;
    int sRotat4 = atoi(getXMLAttribute( xmlSensitive , "rotation4" ).c_str() ) ;
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
        
    // Construct discrete rotation from mounting frame to global frame.  
    Matrix3d DiscreteRotation(Matrix3d::Identity()); 
    DiscreteRotation(0,0) = sRotat1; 
    DiscreteRotation(0,1) = sRotat2; 
    DiscreteRotation(1,0) = sRotat3; 
    DiscreteRotation(1,1) = sRotat4;
         
    if ( sRotat1*sRotat4 - sRotat2*sRotat3 == 1 ) {     
      DiscreteRotation(2,2) = 1;
    } else if ( sRotat1*sRotat4 - sRotat2*sRotat3 == -1 ) {
      DiscreteRotation(2,2) = -1;  
    } else {
      streamlog_out(MESSAGE3) << std::endl << "Rotation parameters in geometry xml file wrong" << std::endl;  
    }
        
    if ( std::abs( DiscreteRotation.determinant() - 1 ) ==  1.e-5 )  
      streamlog_out(MESSAGE3) << "Rotation matrix BUG. Discrete matrix determinant is " << DiscreteRotation.determinant() << std::endl;      
            
    // Construct discrete reference frame 
    ReferenceFrame discrete;
    discrete.SetRotation(DiscreteRotation);   
        
    // Read Euler rotation from local to mounting frame
     
    const double MYPI = std::atan(1.0)*4;  
    // XML file stores angles in degree 
    double alpha = sAlpha*MYPI/180.; 
    double beta = sBeta*MYPI/180.; 
    double gamma = sGamma*MYPI/180.; 
    Matrix3d EulerRotation;
    FillRotMatrixKarimaki(EulerRotation, 
                              alpha, 
                              beta,
                              gamma);
        
    if ( std::abs( EulerRotation.determinant() - 1 ) ==  1.e-5 )  
          streamlog_out(MESSAGE3) << "Rotation matrix BUG. Euler rotation matrix determinant is " << EulerRotation.determinant() << std::endl;   
         
    // Combine the two factors in proper order
    Matrix3d NominalRotation = EulerRotation*DiscreteRotation;
        
    // Read position of origin of local uvw frame in global coordinates 
    Vector3d NominalPosition;
    NominalPosition << sPosX, sPosY, sPosZ;
        
    // Construct nominal reference frame 
    ReferenceFrame nominal;
    nominal.SetPosition(NominalPosition);
    nominal.SetRotation(NominalRotation);     
        
    // Add new Det object to detector, ownership goes with detector
    // The plane number is set to -1 because we do not know the ordering of Dets yet 
    dets.push_back( new SquareDet( "SquareDet",
                                  sensorID,
                                  -1,
                                  sThick, 
                                  sRadLen,
                                  sAtomicNum,
                                  sAtomicMass, 
                                  lThick, 
                                  lRadLen, 
                                  lAtomicNum, 
                                  lAtomicMass,
                                  lSizU, 
                                  lSizV,  
                                  uCellGroupVec, 
                                  vCellGroupVec, 
                                  discrete, 
                                  nominal
                                 ));
     
        
  } // end loop SiPlanes 

}


