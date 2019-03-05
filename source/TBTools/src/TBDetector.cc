// Local include files
#include "TBDetector.h"
#include "ThreeDModel.h"
#include "SquareDet.h"

// Include basic C header files
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cassert> 

// Include ROOT classes
#include <TH1F.h>
#include <TFile.h>

// Include Gear header files
#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>
#include <gear/BField.h>

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

#include <Eigen/LU>

// Namespaces
using namespace marlin;


namespace depfet {


TBDetector::TBDetector( ) 
{
  // This is a setup with no detectors 
  m_alignmentDBFileName = "xxx";
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



void TBDetector::ReadGearConfiguration( )
{
  
  streamlog_out ( MESSAGE3) << "Construct Test Beam Detector" << std::endl;
   
  /*
  //Now create all subcomponents
  for (const GearDir& component : detectorDir.getNodes("DetectorComponent")) {
    string name;
    string creatorName;
    try {
      name        = component.getString("@name");
      creatorName = component.getString("Creator");
    } catch (gearbox::PathEmptyError& e) {
      B2ERROR("Could not find required element Name or Creator for " << component.getPath());
      continue;
    }
    
    if (!m_components.empty() && m_components.count(name) == 0) {
      B2DEBUG(50, "DetectorComponent " << name << " not in list of components, skipping");
      continue;
    }
    
    string libraryName = component.getString("Creator/@library", "");
    if (!iov.empty()) {
      CreatorBase* creator = CreatorManager::getCreator(creatorName, libraryName);
      if (creator) {
        creator->createPayloads(GearDir(component, "Content"), iov);
      } else {
        B2ERROR("Could not load creator " << creatorName << " from " << libraryName);
      }
    }
    config.addComponent({name, creatorName, libraryName});
  }
  */

  
  // Check iff gear file is available  
  if ( Global::GEAR == 0x0 ) {
    streamlog_out ( ERROR4 ) <<  "The GearMgr is not available, for an unknown reason." << std::endl;
    exit(-1);
  }
  
  // 
  // Read data about constant magnetic field vector 
   
  gear::BField * bField = 
                  const_cast<gear::BField* > (&(Global::GEAR->getBField()));  

  m_Bx = bField->at(gear::Vector3D())[0];
  m_By = bField->at(gear::Vector3D())[1];
  m_Bz = bField->at(gear::Vector3D())[2];
  
  //
  // Read data about detectors planes 

  gear::SiPlanesParameters * siPlanesParameters = 
                  const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters())); 
  
  gear::SiPlanesLayerLayout * siPlanesLayerLayout = 
                  const_cast<gear::SiPlanesLayerLayout*> ( &(siPlanesParameters->getSiPlanesLayerLayout() ));
     
  // 
  // Sorting all sensors according to position along beam 
  
  m_numberOfSensors = siPlanesLayerLayout->getNLayers();    
  std::vector<int> planeSort;  
  std::vector<double> planeZPosition;
   
  for(int ipl=0; ipl < m_numberOfSensors; ipl++) {
    planeZPosition.push_back( siPlanesLayerLayout->getSensitivePositionZ(ipl) );
    planeSort.push_back( ipl ); 
  }
     
  bool sorted;
  do {
     sorted=false;
     for(int iz=0; iz<m_numberOfSensors-1 ; iz++) {
       if( planeZPosition[iz]> planeZPosition[iz+1])   {
         double posZ = planeZPosition[iz];
         planeZPosition[iz] = planeZPosition[iz+1];
         planeZPosition[iz+1] = posZ;
          
         int idZ = planeSort[iz];
         planeSort[iz] = planeSort[iz+1];
         planeSort[iz+1] = idZ;
         
         sorted=true;
       }
     }
  } while(sorted);
  
  
   
  // Loop over all pixel modules  
  for (int ipl=0; ipl < m_numberOfSensors ; ipl++) {
    
    
    // TODO: all code in this loop is currently SquareDet specific
    //       and needs to be move into a SquareDetCreator class

    // Gear layer number is needed to read the layout data  
    int ilayer = planeSort[ipl];
     
    // Set sensorID   
    int sensorID = siPlanesLayerLayout->getSensitiveID(ilayer);
    m_indexMap[sensorID] = ipl;

     
    // Read discrete rotation from mounting frame to global frame.  
    
    Matrix3d DiscreteRotation(Matrix3d::Identity()); 
    int r1 = siPlanesLayerLayout->getSensitiveRotation1(ilayer); 
    int r2 = siPlanesLayerLayout->getSensitiveRotation2(ilayer); 
    int r3 = siPlanesLayerLayout->getSensitiveRotation3(ilayer);  
    int r4 = siPlanesLayerLayout->getSensitiveRotation4(ilayer);
    
    DiscreteRotation(0,0) = r1; 
    DiscreteRotation(0,1) = r2; 
    DiscreteRotation(1,0) = r3; 
    DiscreteRotation(1,1) = r4;
    
    if ( r1*r4 - r2*r3 == 1 ) {     
      DiscreteRotation(2,2) = 1;
    } else if ( r1*r4 - r2*r3 == -1 ) {
      DiscreteRotation(2,2) = -1;  
    } else {
      streamlog_out(MESSAGE3) << std::endl << "Rotation parameters in gear file wrong" << std::endl;  
    }
    
    if ( std::abs( DiscreteRotation.determinant() - 1 ) ==  1.e-5 )  
      streamlog_out(MESSAGE3) << "Rotation matrix BUG. Discrete matrix determinant is " << DiscreteRotation.determinant() << std::endl;      
    
    // Construct discrete reference frame 
    ReferenceFrame discrete;
    discrete.SetRotation(DiscreteRotation);   
    
    // Read Euler rotation from local to mounting frame
     
    const double MYPI = std::atan(1.0)*4;  
    // Gear file stores angles in degree 
    double alpha = siPlanesLayerLayout->getSensitiveRotationAlpha(ilayer)*MYPI/180.; 
    double beta = siPlanesLayerLayout->getSensitiveRotationBeta(ilayer)*MYPI/180.; 
    double gamma = siPlanesLayerLayout->getSensitiveRotationGamma(ilayer)*MYPI/180.; 
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
    NominalPosition << siPlanesLayerLayout->getSensitivePositionX(ilayer),  siPlanesLayerLayout->getSensitivePositionY(ilayer), siPlanesLayerLayout->getSensitivePositionZ(ilayer);
    
    // Construct nominal reference frame 
    ReferenceFrame nominal;
    nominal.SetPosition(NominalPosition);
    nominal.SetRotation(NominalRotation);   
    
    // Create a new Det object, ownership goes with TBDetector 
    m_Dets.push_back( new SquareDet( "SquareDet",
                                  sensorID,
                                  ipl,
                                  siPlanesLayerLayout->getSensitiveThickness(ilayer), 
                                  siPlanesLayerLayout->getSensitiveRadLength(ilayer),
                                  siPlanesLayerLayout->getSensitiveAtomicNumber(ilayer),
                                  siPlanesLayerLayout->getSensitiveAtomicMass(ilayer), 
                                  siPlanesLayerLayout->getLayerThickness(ilayer), 
                                  siPlanesLayerLayout->getLayerRadLength(ilayer), 
                                  siPlanesLayerLayout->getLayerAtomicNumber(ilayer), 
                                  siPlanesLayerLayout->getLayerAtomicMass(ilayer),
                                  siPlanesLayerLayout->getLayerSizeU(ilayer), 
                                  siPlanesLayerLayout->getLayerSizeV(ilayer),  
                                  siPlanesLayerLayout->getSensitiveUCells( ilayer ), 
                                  siPlanesLayerLayout->getSensitiveVCells( ilayer ), 
                                  discrete, 
                                  nominal
                                 ));
  }
  
    
}


void TBDetector::SetAlignmentDBName( std::string FileName )
{
     
  // Store name of alignment data base
  m_alignmentDBFileName = FileName;
}


void TBDetector::ReadAlignmentDB( std::string FileName )
{
     
  // Store name of alignment data base
  m_alignmentDBFileName = FileName;
  
  // Open alignment data base
  TFile * rootFile = new TFile(m_alignmentDBFileName.c_str(), "READ");
  
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
  
  streamlog_out(MESSAGE3) << std::endl << "Write alignment DB file " << m_alignmentDBFileName << std::endl << std::endl;
  
  TFile * rootFile = new TFile( m_alignmentDBFileName.c_str(),"recreate");
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

  streamlog_out(MESSAGE3) << std::endl << "Finish alignment DB file " << m_alignmentDBFileName << std::endl << std::endl;
}


Det & TBDetector::GetDet( int ipl ) 
{
  return reinterpret_cast<Det&>(m_Dets[ipl]);
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

