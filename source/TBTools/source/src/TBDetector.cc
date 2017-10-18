// Local include files
#include "TBDetector.h"
#include "ThreeDModel.h"

// Include basic C header files
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cassert> 

// Include ROOT classes
#include <TH1F.h>
#include <TFile.h>

// Include CLHEP header files
#include <CLHEP/Matrix/DiagMatrix.h>

// Include Gear header files
#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>
#include "gear/BField.h"

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

// Namespaces
using namespace CLHEP;
using namespace marlin;


namespace depfet {

//
// Constructor
//
TBDetector::TBDetector( ) 
{
  // This is a setup with no detectors 
  _alignmentDBFileName = "xxx";
  _numberOfSensors  = 0; 

  _Bx = 0;
  _By = 0;
  _Bz = 0;
  
}

//
// Destructor - Delete the detectors from memory
//
TBDetector::~TBDetector()
{
}


//
//  Build detector from gear file
//
void TBDetector::ReadGearConfiguration( )
{
  
  streamlog_out ( MESSAGE3) << "Construct Test Beam Detector" << std::endl;
    
  // Check iff gear file is available  
  if ( Global::GEAR == 0x0 ) {
    streamlog_out ( ERROR4 ) <<  "The GearMgr is not available, for an unknown reason." << std::endl;
    exit(-1);
  }
  
  // 
  // Read data about constant magnetic field vector 
   
  gear::BField * bField = 
                  const_cast<gear::BField* > (&(Global::GEAR->getBField()));  

  _Bx = bField->at(gear::Vector3D())[0];
  _By = bField->at(gear::Vector3D())[1];
  _Bz = bField->at(gear::Vector3D())[2];
  
  //
  // Read data about detectors planes 

  gear::SiPlanesParameters * siPlanesParameters = 
                  const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters())); 
  
  gear::SiPlanesLayerLayout * siPlanesLayerLayout = 
                  const_cast<gear::SiPlanesLayerLayout*> ( &(siPlanesParameters->getSiPlanesLayerLayout() ));
     
  // 
  // Sorting all sensors according to position along beam 
  
  _numberOfSensors = siPlanesLayerLayout->getNLayers();    
  int * planeSort = new int[_numberOfSensors];
  double * planeZPosition = new double[_numberOfSensors];
   
  for(int ipl=0; ipl < _numberOfSensors; ipl++) {
    planeZPosition[ipl]=  siPlanesLayerLayout->getSensitivePositionZ(ipl);
    planeSort[ipl]=ipl; 
  }
   
  bool sorted;
  do {
     sorted=false;
     for(int iz=0; iz<_numberOfSensors-1 ; iz++) {
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
  
  // Resize space for pixel modules (Det's) 
  _DetVec.resize( _numberOfSensors );
   
  // Loop over all pixel modules  
  for (int ipl=0; ipl < _numberOfSensors ; ipl++) {
    
    // Create a new Det object  
    Det adet;

    // Set plane number 
    adet.SetPlaneNumber(ipl);

    // Gear layer number is needed to read the layout data  
    int ilayer = planeSort[ipl];
    
    // Set device type 
    int Type = siPlanesLayerLayout->getSensitivePixelType(ilayer);  
    adet.SetDeviceType(Type);
     
    // Set DAQ ID   
    int ID = siPlanesLayerLayout->getSensitiveID(ilayer);
    adet.SetDAQID(ID);
    _indexMap[ID] = ipl;
      
    // Set other stuff  
    int NColumns = siPlanesLayerLayout->getSensitiveNpixelX(ilayer);
    adet.SetNColumns(NColumns); 
    
    int NRows = siPlanesLayerLayout->getSensitiveNpixelY(ilayer); 
    adet.SetNRows(NRows);
      
    double PitchU = siPlanesLayerLayout->getSensitivePitchX(ilayer);
    adet.SetPitchU(PitchU);
    
    double PitchV = siPlanesLayerLayout->getSensitivePitchY(ilayer);
    adet.SetPitchV(PitchV);
    
    double ResoU = siPlanesLayerLayout->getSensitiveResolutionX(ilayer);
    adet.SetResolutionU(ResoU);
    
    double ResoV = siPlanesLayerLayout->getSensitiveResolutionY(ilayer);
    adet.SetResolutionV(ResoV);
     
    double SensThick = siPlanesLayerLayout->getSensitiveThickness(ilayer);
    adet.SetSensitiveThickness(SensThick);
    
    double LadderThick = siPlanesLayerLayout->getLayerThickness(ilayer);
    adet.SetLadderThickness(LadderThick);
        
    double SensRadLength = siPlanesLayerLayout->getSensitiveRadLength(ilayer);
    adet.SetSensitiveRadLength(SensRadLength);
       
    double LadderRadLength = siPlanesLayerLayout->getLayerRadLength(ilayer);
    adet.SetLadderRadLength(LadderRadLength);
      
    double LadderSizeX= siPlanesLayerLayout->getLayerSizeX(ilayer);
    adet.SetModuleBoxSizeU(LadderSizeX); 
         
    double LadderSizeY = siPlanesLayerLayout->getLayerSizeY(ilayer); 
    adet.SetModuleBoxSizeV(LadderSizeY); 
    
    double SensSizeX = siPlanesLayerLayout->getSensitiveSizeX(ilayer); 
    adet.SetSensitiveSizeU(SensSizeX); 
    
    double SensSizeY = siPlanesLayerLayout->getSensitiveSizeY(ilayer); 
    adet.SetSensitiveSizeV(SensSizeY); 
       
    // Discrete rotation frame 
    ReferenceFrame discrete;
    
    // Construct a nominal local <-> global rotation matrix
    // from gear file. 
    
    // Read discrete rotation from local to global coord.  
    HepMatrix DiscreteRotation(3, 3, 1); 
    int r1 = siPlanesLayerLayout->getSensitiveRotation1(ilayer); 
    int r2 = siPlanesLayerLayout->getSensitiveRotation2(ilayer); 
    int r3 = siPlanesLayerLayout->getSensitiveRotation3(ilayer);  
    int r4 = siPlanesLayerLayout->getSensitiveRotation4(ilayer);
    
    DiscreteRotation[0][0] = r1; 
    DiscreteRotation[0][1] = r2; 
    DiscreteRotation[1][0] = r3; 
    DiscreteRotation[1][1] = r4;
    
    if ( r1*r4 - r2*r3 == 1 ) {     
      DiscreteRotation[2][2] = 1;
    } else if ( r1*r4 - r2*r3 == -1 ) {
      DiscreteRotation[2][2] = -1;  
    } else {
      streamlog_out(MESSAGE3) << std::endl << "Rotation parameters in gear file wrong" << std::endl;  
    }
    
    if ( std::abs( DiscreteRotation.determinant() - 1 ) ==  1.e-5 )  
      streamlog_out(MESSAGE3) << "Rotation matrix BUG. Discrete matrix determinant is " << DiscreteRotation.determinant() << std::endl; 
    
    discrete.SetRotation(DiscreteRotation);    
    
    // Set discrete frame 
    adet.SetDiscreteFrame(discrete);
    
    // Nominal reference frame 
    ReferenceFrame nominal;
    
    HepVector NominalPosition(3);
    NominalPosition[0] = siPlanesLayerLayout->getSensitivePositionX(ilayer); 
    NominalPosition[1] = siPlanesLayerLayout->getSensitivePositionY(ilayer);
    NominalPosition[2] = siPlanesLayerLayout->getSensitivePositionZ(ilayer);
    nominal.SetPosition(NominalPosition);


    // Read Euler rotation from local to global frame
    HepMatrix EulerRotation;
    const double MYPI = std::atan(1.0)*4;  
    // Gear file stores angles in degree 
    double alpha = siPlanesLayerLayout->getSensitiveRotationAlpha(ilayer)*MYPI/180.; 
    double beta = siPlanesLayerLayout->getSensitiveRotationBeta(ilayer)*MYPI/180.; 
    double gamma = siPlanesLayerLayout->getSensitiveRotationGamma(ilayer)*MYPI/180.; 
    FillRotMatrixKarimaki(EulerRotation, 
            alpha, 
            beta,
            gamma);
    
    if ( std::abs( EulerRotation.determinant() - 1 ) ==  1.e-5 )  
      streamlog_out(MESSAGE3) << "Rotation matrix BUG. Euler rotation matrix determinant is " << EulerRotation.determinant() << std::endl; 
    
    // Combine the two factors in proper order
    HepMatrix NominalRotation = EulerRotation*DiscreteRotation;
    nominal.SetRotation(NominalRotation);    
    
    // Set nominal frame - initial guess where detector is in space 
    adet.SetNominalFrame(nominal);
    
    // Save detector in collection
    _DetVec[ipl] = adet;
    
  }
  
  // Clean up 
  delete [] planeSort ; 
  delete [] planeZPosition ;   
  
}

// 
// Read name of alignment data base file 
//   
void TBDetector::SetAlignmentDBName( std::string FileName )
{
     
  // Store name of alignment data base
  _alignmentDBFileName = FileName;
}

// 
// Read alignment data base file 
//   
void TBDetector::ReadAlignmentDB( std::string FileName )
{
     
  // Store name of alignment data base
  _alignmentDBFileName = FileName;
  
  // Open alignment data base
  TFile * rootFile = new TFile(_alignmentDBFileName.c_str(), "READ");
  
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

  for ( int ipl = 0; ipl < _numberOfSensors; ipl++ ) {
    
    Det & adet = GetDet(ipl);
    
    int bin = ipl + 1; 

    // Check consistency with gear file 
    if ( histoMap["hSensorID"]->GetBinContent(bin) != adet.GetDAQID() ) {
      streamlog_out ( WARNING ) <<  "Alignment DB inconsistent to Gear file!!" << std::endl;   
    }
    
    // Now, read nominal position from DB
    ReferenceFrame nominal;
    HepVector NominalPosition(3);
    NominalPosition[0] = histoMap["hPositionX"]->GetBinContent(bin);
    NominalPosition[1] = histoMap["hPositionY"]->GetBinContent(bin);
    NominalPosition[2] = histoMap["hPositionZ"]->GetBinContent(bin);
    nominal.SetPosition(NominalPosition);
      
    // AlignmentDB stores Euler angles in rad 
    HepMatrix EulerRotation;
    double alpha = histoMap["hRotationAlpha"]->GetBinContent(bin);   
    double beta  = histoMap["hRotationBeta"]->GetBinContent(bin);  
    double gamma = histoMap["hRotationGamma"]->GetBinContent(bin);  
    FillRotMatrixKarimaki(EulerRotation, alpha, beta, gamma);
      
    // Combine the two factors in proper order
    HepMatrix DiscreteRotation = adet.GetDiscrete().GetRotation(); 
    HepMatrix NominalRotation = EulerRotation*DiscreteRotation;
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

// 
// Write alignment data base - overwrites old DB file   
//   
void TBDetector::WriteAlignmentDB( )
{
  
  streamlog_out(MESSAGE3) << std::endl << "Write alignment DB file " << _alignmentDBFileName << std::endl << std::endl;
  
  TFile * rootFile = new TFile( _alignmentDBFileName.c_str(),"recreate");
  rootFile->cd("");
  
  std::map< std::string, TH1F *> _histoMap;
   
  _histoMap["hSensorID"] = new TH1F("hSensorID", "", _numberOfSensors,0,_numberOfSensors); 
  _histoMap["hSensorID"]->SetStats( false );
  _histoMap["hSensorID"]->SetYTitle("sensor id"); 
  _histoMap["hSensorID"]->SetXTitle("plane");

  _histoMap["hPositionX"] = new TH1F("hPositionX", "", _numberOfSensors,0,_numberOfSensors); 
  _histoMap["hPositionX"]->SetStats( false );
  _histoMap["hPositionX"]->SetYTitle("position x [mm]"); 
  _histoMap["hPositionX"]->SetXTitle("plane"); 

  _histoMap["hPositionY"] = new TH1F("hPositionY", "", _numberOfSensors,0,_numberOfSensors); 
  _histoMap["hPositionY"]->SetStats( false );
  _histoMap["hPositionY"]->SetYTitle("position y [mm]"); 
  _histoMap["hPositionY"]->SetXTitle("plane"); 

  _histoMap["hPositionZ"] = new TH1F("hPositionZ", "", _numberOfSensors,0,_numberOfSensors); 
  _histoMap["hPositionZ"]->SetStats( false );
  _histoMap["hPositionZ"]->SetYTitle("position z [mm]"); 
  _histoMap["hPositionZ"]->SetXTitle("plane"); 

  _histoMap["hRotationAlpha"] = new TH1F("hRotationAlpha", "", _numberOfSensors,0,_numberOfSensors); 
  _histoMap["hRotationAlpha"]->SetStats( false );
  _histoMap["hRotationAlpha"]->SetYTitle("rotation alpha [rad]"); 
  _histoMap["hRotationAlpha"]->SetXTitle("plane"); 

  _histoMap["hRotationBeta"] = new TH1F("hRotationBeta", "", _numberOfSensors,0,_numberOfSensors); 
  _histoMap["hRotationBeta"]->SetStats( false );
  _histoMap["hRotationBeta"]->SetYTitle("rotation beta [rad]"); 
  _histoMap["hRotationBeta"]->SetXTitle("plane"); 

  _histoMap["hRotationGamma"] = new TH1F("hRotationGamma", "", _numberOfSensors,0,_numberOfSensors); 
  _histoMap["hRotationGamma"]->SetStats( false );
  _histoMap["hRotationGamma"]->SetYTitle("rotation gamma [rad]"); 
  _histoMap["hRotationGamma"]->SetXTitle("plane");

  for ( int ipl = 0; ipl < _numberOfSensors; ipl++ ) {
    
    // Load current pixel module 
    Det & adet = GetDet(ipl);    
    
    // Now, write nominal tracker geometry parameters  
    HepVector Origin = adet.GetNominal().GetPosition(); 
    HepMatrix NominalRotation = adet.GetNominal().GetRotation(); 
      
    // Factor out the discrete rotation 
    HepMatrix DiscreteRotation = adet.GetDiscrete().GetRotation(); 
    HepMatrix EulerRotation = NominalRotation*DiscreteRotation.T();
       
    double alpha, beta, gamma; 
    GetAnglesKarimaki(EulerRotation, alpha, beta, gamma);   
    
    int bin = ipl+1;
    
    _histoMap["hSensorID"]->SetBinContent( bin, adet.GetDAQID() );
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

  streamlog_out(MESSAGE3) << std::endl << "Finish alignment DB file " << _alignmentDBFileName << std::endl << std::endl;
}

/** Get data handle for pixel module at position ipl
 */
Det & TBDetector::GetDet( int ipl ) 
{
  return _DetVec[ipl];
}


/** Get plane number from sensorID 
 */  
int TBDetector::GetPlaneNumber(int ID) 
{
  if (_indexMap.count(ID) > 0 ) 
    return _indexMap[ID];
  else 
    return -99;  
}

/** Method printing detector
 */
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

