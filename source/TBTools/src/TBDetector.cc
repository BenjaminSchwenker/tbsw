// Local include files
#include "TBDetector.h"
#include "ThreeDModel.h"
#include "DetCreatorManager.h"
#include "DetCreatorBase.h"

// Include basic C header files
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <memory>

// Include ROOT classes
#include <TH1F.h>
#include <TFile.h>

// Include Gear header files
#include <gear/GEAR.h>


// Include Marlin
#include <marlin/Global.h>
#include <marlin/tinyxml.h>
#include <streamlog/streamlog.h>

#include <Eigen/LU>

// Namespaces
using namespace marlin;


namespace depfet {

// Helper function copied over from gear (this all that is actually needed)
namespace {
  
  std::string getXMLAttribute(const TiXmlNode* node , const std::string& name ) {
    
    const TiXmlElement* el = node->ToElement() ;
    if( el == 0 )
      throw gear::ParseException("XMLParser::getAttribute not an XMLElement " ) ;
    
    const std::string*  at = el->Attribute( name ) ;
    
    if( at == 0 ){
      std::stringstream str ;
      str  << "XMLParser::getAttribute missing attribute \"" << name
	   << "\" in element <" << el->Value() << "/> " ;
      throw gear::ParseException( str.str() ) ;
    }
     
    return std::string( *at )  ;
  }
 
} // end anonymouse namespace


TBDetector& TBDetector::GetInstance()
{
  static std::unique_ptr<TBDetector> instance(new TBDetector());
  return *instance;
}

TBDetector::TBDetector( ) 
{
  // This is a setup with no detectors 
  m_alignmentDBFilePath = "";
  m_numberOfSensors  = 0; 

  m_Bx = 0;
  m_By = 0;
  m_Bz = 0;
}

TBDetector::~TBDetector()
{
  for (Det* aDet : m_Dets) delete aDet;
  m_Dets.clear();
}

void TBDetector::ReadGearConfiguration( const std::string & geometryXMLFile )
{
  
  streamlog_out ( MESSAGE3) << "Construct Test Beam Detector NEW WAY" << std::endl;
   
  TiXmlDocument doc ;
  bool loadOkay = doc.LoadFile( geometryXMLFile  ) ;
  
  if( !loadOkay ) {
      
    std::stringstream str ;
    str  << "Error reading geometry XML file [" << geometryXMLFile
	     << ", row: " << doc.ErrorRow() << ", col: " << doc.ErrorCol() << "] : "
	     << doc.ErrorDesc() ;
      
    throw gear::ParseException( str.str() ) ;
  }
  
  TiXmlElement* root = doc.RootElement();
   
  if( root == 0 ){
    throw gear::ParseException( std::string( "No root tag found in ") + geometryXMLFile  ) ;
  }
   
  TiXmlNode* global = root->FirstChild("global") ;
  if( global != 0 ){
    std::string detName  =  getXMLAttribute( global, "detectorName" )  ;
    // TODO put the detName somewhere
    streamlog_out ( MESSAGE3) << "detectorName " << detName  << std::endl;
  }
     
  // 
  // Read data about constant magnetic field vector 
   
  TiXmlNode* field = root->FirstChild("BField")  ;
  if( field != 0 ){
    m_Bx  =  atof(  getXMLAttribute( field , "x" ) .c_str() ) ;
    m_By  =  atof(  getXMLAttribute( field , "y" ) .c_str() ) ;
    m_Bz  =  atof(  getXMLAttribute( field , "z" ) .c_str() ) ;
  }
  
  //
  // Read data about detectors planes 
   
  TiXmlNode* detectors = root->FirstChild("detectors")  ;
  if( detectors == 0 ){
    throw gear::ParseException( std::string( "No detectors tag found in  ") + geometryXMLFile  ) ;
  }
   
  TiXmlNode* det = 0 ;
  while( ( det = detectors->IterateChildren( "detector", det ) )  != 0  ){
    
    std::string name  =  getXMLAttribute( det, "name" );  
    std::string creatorName("UNKOWN") ;
    
    try {       
	  creatorName =  getXMLAttribute( det, "geartype" )  ;
       
      //std::cout << "Reading detector " << name
      //          << " with \"geartype\" " << creatorName << std::endl ;
    } catch( gear::ParseException& e){
      
	  streamlog_out ( MESSAGE3) << "Igoring detector " << name
                                << " with missing attribute \"geartype\" " << std::endl ; 
	  continue ;
    }
     
    DetCreatorBase* creator = DetCreatorManager::getCreator(creatorName);
    if (creator) {
      creator->create(det, m_Dets);
    } else {
      streamlog_out (ERROR3) << "Could not find creator " << creatorName << std::endl;
    }
    
  } // end loop detector components
    
  // Set total number of sensors read from XML file
  m_numberOfSensors = int(m_Dets.size());
      
  // Order the sensors according the the z position along 
  // the beam line.
  std::sort(std::begin(m_Dets), std::end(m_Dets), [](auto const *t1, auto const *t2) {
                   return t1->GetNominal().GetZPosition() < t2->GetNominal().GetZPosition(); });
      
  int planeNumber = 0;
  for (Det* aDet : m_Dets) {
    aDet->SetPlaneNumber(planeNumber);
    m_indexMap[aDet->GetSensorID()] = planeNumber;
    planeNumber++;  
  }
}



void TBDetector::SetAlignmentDBPath( std::string FilePath )
{
  // Store name of alignment data base
  m_alignmentDBFilePath = FilePath;
}


void TBDetector::ApplyAlignmentDB(  )
{
    
  // Open alignment data base file
  TFile * rootFile = new TFile(m_alignmentDBFilePath.c_str(), "READ");
  
  std::map< std::string, TH1F *> histoMap;
 
 
  if ( (TH1F *) rootFile->Get("hSensorID") != nullptr) {
    histoMap["hSensorID"] = (TH1F *) rootFile->Get("hSensorID");  
  } else { 
    streamlog_out ( WARNING ) <<  "Alignment DB empty!!" << std::endl;
    return; 
  } 
  
  if ( (TH1F *) rootFile->Get("hPositionX") != nullptr) {
    histoMap["hPositionX"] = (TH1F *) rootFile->Get("hPositionX");  
  } else { 
    streamlog_out ( WARNING ) <<  "Alignment DB empty!!" << std::endl;
    return; 
  } 

  if ( (TH1F *) rootFile->Get("hPositionY") != nullptr) {
    histoMap["hPositionY"] = (TH1F *) rootFile->Get("hPositionY");  
  } else { 
    streamlog_out ( WARNING ) <<  "Alignment DB empty!!" << std::endl;
    return; 
  } 

  if ( (TH1F *) rootFile->Get("hPositionZ") != nullptr) {
    histoMap["hPositionZ"] = (TH1F *) rootFile->Get("hPositionZ");  
  } else {
    streamlog_out ( WARNING ) <<  "Alignment DB empty!!" << std::endl; 
    return; 
  } 

  if ( (TH1F *) rootFile->Get("hRotationAlpha") != nullptr) {
    histoMap["hRotationAlpha"] = (TH1F *) rootFile->Get("hRotationAlpha");  
  } else {
    streamlog_out ( WARNING ) <<  "Alignment DB empty!!" << std::endl; 
    return; 
  } 

  if ( (TH1F *) rootFile->Get("hRotationBeta") != nullptr) {
    histoMap["hRotationBeta"] = (TH1F *) rootFile->Get("hRotationBeta");  
  } else {
    streamlog_out ( WARNING ) <<  "Alignment DB empty!!" << std::endl; 
    return; 
  } 

  if ( (TH1F *) rootFile->Get("hRotationGamma") != nullptr) {
    histoMap["hRotationGamma"] = (TH1F *) rootFile->Get("hRotationGamma");  
  } else {
    streamlog_out ( WARNING ) <<  "Alignment DB empty!!" << std::endl; 
    return; 
  } 

  for ( int ipl = 0; ipl < m_numberOfSensors; ipl++ ) {
    
    Det & adet = GetDet(ipl);
    
    int bin = ipl + 1; 

    // Check consistency with gear file 
    if ( histoMap["hSensorID"]->GetBinContent(bin) != adet.GetSensorID() ) {
      streamlog_out ( WARNING ) <<  "Alignment DB inconsistent to Gear file!!" << std::endl;   
    }
    
    // Now, read nominal position from DB
    ReferenceFrame nominal;
    Vector3d NominalPosition;
    NominalPosition << histoMap["hPositionX"]->GetBinContent(bin), histoMap["hPositionY"]->GetBinContent(bin), histoMap["hPositionZ"]->GetBinContent(bin);
    nominal.SetPosition(NominalPosition);
      
    // AlignmentDB stores Euler angles in rad 
    Matrix3d EulerRotation;
    double alpha = histoMap["hRotationAlpha"]->GetBinContent(bin);   
    double beta  = histoMap["hRotationBeta"]->GetBinContent(bin);  
    double gamma = histoMap["hRotationGamma"]->GetBinContent(bin);  
    FillRotMatrixKarimaki(EulerRotation, alpha, beta, gamma);
      
    // Combine the two factors in proper order
    Matrix3d DiscreteRotation = adet.GetDiscrete().GetRotation(); 
    Matrix3d NominalRotation = EulerRotation*DiscreteRotation;
    nominal.SetRotation(NominalRotation);  
      
    nominal.SetRotation(NominalRotation);
    adet.SetNominalFrame(nominal); 
    
  }
   
  // Close root  file
  rootFile->Close();
  delete rootFile;
  
  // Print detector data sheets  
  Print();    
}


void TBDetector::WriteAlignmentDB( )
{
  
  streamlog_out(MESSAGE3) << std::endl << "Write alignment DB file " << m_alignmentDBFilePath << std::endl << std::endl;
  
  TFile * rootFile = new TFile( m_alignmentDBFilePath.c_str(),"recreate");
  rootFile->cd("");
  
  std::map< std::string, TH1F *> _histoMap;
   
  _histoMap["hSensorID"] = new TH1F("hSensorID", "", m_numberOfSensors,0,m_numberOfSensors); 
  _histoMap["hSensorID"]->SetStats( false );
  _histoMap["hSensorID"]->SetYTitle("sensor id"); 
  _histoMap["hSensorID"]->SetXTitle("plane");

  _histoMap["hPositionX"] = new TH1F("hPositionX", "", m_numberOfSensors,0,m_numberOfSensors); 
  _histoMap["hPositionX"]->SetStats( false );
  _histoMap["hPositionX"]->SetYTitle("position x [mm]"); 
  _histoMap["hPositionX"]->SetXTitle("plane"); 

  _histoMap["hPositionY"] = new TH1F("hPositionY", "", m_numberOfSensors,0,m_numberOfSensors); 
  _histoMap["hPositionY"]->SetStats( false );
  _histoMap["hPositionY"]->SetYTitle("position y [mm]"); 
  _histoMap["hPositionY"]->SetXTitle("plane"); 

  _histoMap["hPositionZ"] = new TH1F("hPositionZ", "", m_numberOfSensors,0,m_numberOfSensors); 
  _histoMap["hPositionZ"]->SetStats( false );
  _histoMap["hPositionZ"]->SetYTitle("position z [mm]"); 
  _histoMap["hPositionZ"]->SetXTitle("plane"); 

  _histoMap["hRotationAlpha"] = new TH1F("hRotationAlpha", "", m_numberOfSensors,0,m_numberOfSensors); 
  _histoMap["hRotationAlpha"]->SetStats( false );
  _histoMap["hRotationAlpha"]->SetYTitle("rotation alpha [rad]"); 
  _histoMap["hRotationAlpha"]->SetXTitle("plane"); 

  _histoMap["hRotationBeta"] = new TH1F("hRotationBeta", "", m_numberOfSensors,0,m_numberOfSensors); 
  _histoMap["hRotationBeta"]->SetStats( false );
  _histoMap["hRotationBeta"]->SetYTitle("rotation beta [rad]"); 
  _histoMap["hRotationBeta"]->SetXTitle("plane"); 

  _histoMap["hRotationGamma"] = new TH1F("hRotationGamma", "", m_numberOfSensors,0,m_numberOfSensors); 
  _histoMap["hRotationGamma"]->SetStats( false );
  _histoMap["hRotationGamma"]->SetYTitle("rotation gamma [rad]"); 
  _histoMap["hRotationGamma"]->SetXTitle("plane");

  for ( int ipl = 0; ipl < m_numberOfSensors; ipl++ ) {
    
    // Load current pixel module 
    Det & adet = GetDet(ipl);    
    
    // Now, write nominal tracker geometry parameters  
    auto Origin = adet.GetNominal().GetPosition(); 
    auto NominalRotation = adet.GetNominal().GetRotation(); 
      
    // Factor out the discrete rotation 
    auto DiscreteRotation = adet.GetDiscrete().GetRotation(); 
    auto EulerRotation = NominalRotation*DiscreteRotation.transpose();
       
    double alpha, beta, gamma; 
    GetAnglesKarimaki(EulerRotation, alpha, beta, gamma);   
    
    int bin = ipl+1;
    
    _histoMap["hSensorID"]->SetBinContent( bin, adet.GetSensorID() );
    _histoMap["hSensorID"]->SetBinError( bin, 0 );
    
    _histoMap["hPositionX"]->SetBinContent( bin, Origin[0] );
    _histoMap["hPositionX"]->SetBinError( bin, 0 );
     
    _histoMap["hPositionY"]->SetBinContent( bin, Origin[1] );
    _histoMap["hPositionY"]->SetBinError( bin, 0 ); 

    _histoMap["hPositionZ"]->SetBinContent( bin, Origin[2] );
    _histoMap["hPositionZ"]->SetBinError( bin, 0 ); 
      
    _histoMap["hRotationAlpha"]->SetBinContent( bin, alpha );
    _histoMap["hRotationAlpha"]->SetBinError( bin, 0 );  

    _histoMap["hRotationBeta"]->SetBinContent( bin, beta );
    _histoMap["hRotationBeta"]->SetBinError( bin, 0 );   

    _histoMap["hRotationGamma"]->SetBinContent( bin, gamma );
    _histoMap["hRotationGamma"]->SetBinError( bin, 0 );    
     
  }
  
  // Close root  file
  rootFile->Write();
  rootFile->Close();
  delete rootFile;   

  streamlog_out(MESSAGE3) << std::endl << "Finish alignment DB file " << m_alignmentDBFilePath << std::endl << std::endl;
}


Det & TBDetector::GetDet( int ipl ) 
{
  return *m_Dets[ipl];
}

const Det & TBDetector::GetDet( int ipl ) const 
{
  return *m_Dets[ipl];
}



int TBDetector::GetPlaneNumber(int sensorID) const 
{
  if (m_indexMap.count(sensorID) > 0 ) 
    return m_indexMap.at(sensorID);
  else 
    return -99;  
}


void TBDetector::Print( ) 
{
  streamlog_out(MESSAGE3) << std::endl
                          << " "
                          << "Pixel Tracking Telescope Parameters: "
                          << " "
                          << std::endl  << std::endl;
  
  
  for (int ipl=0; ipl< GetNSensors(); ++ipl){
    Det & adet = GetDet(ipl);
     
    // Print general detector data 
    adet.Print();
    
    // Print nominal reference frame
    adet.GetNominal().PrintHepMatrix(); 
  }  
  
}

} // Namespace;

