// Local include files
#include "TBDetector.h"
#include "ThreeDModel.h"

// Include basic C header files
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cassert> 

// Include CLHEP header files
#include <CLHEP/Matrix/DiagMatrix.h>

// Include Gear header files
#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>
#include "gear/BField.h"


// Include LCIO header files 
#include "lcio.h"
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCFlagImpl.h>

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

// Namespaces
using namespace CLHEP;
using namespace marlin;
using namespace lcio;

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
  
  // Check iff external db file is available 
  LCReader * lcReader = LCFactory::getInstance()->createLCReader();  
  
  try {
    
    lcReader->open( _alignmentDBFileName );
    LCEvent* alignEvt = lcReader->readNextEvent();
     
    LCCollectionVec * alignVec = static_cast<  LCCollectionVec * >  ( alignEvt->getCollection( "alignDB" )); 
          
    for ( int ipl = 0; ipl < _numberOfSensors; ipl++ ) 
    {
      
      // Load alignment constants from DB 
      AlignmentConstant * align = static_cast< AlignmentConstant * > ( alignVec->getElementAt( ipl ) );
      
      // Load current pixel module    
      Det & adet = GetDet(ipl);
      
      // Check consistency with gear file 
      if ( align->getSensorID() != adet.GetDAQID() ) {
        streamlog_out ( ERROR4 ) <<  "Alignment DB inconsistent to Gear file!!" << std::endl;
        exit(-1);
      }
      
      // Now, read nominal position from DB
      ReferenceFrame nominal;
      HepVector NominalPosition(3);
      NominalPosition[0] = align->getXOffset() ;
      NominalPosition[1] = align->getYOffset() ;
      NominalPosition[2] = align->getZOffset() ;
      nominal.SetPosition(NominalPosition);
      
      // AlignmentDB stores Euler angles in rad 
      HepMatrix EulerRotation;
      double alpha = align->getAlpha() ;
      double beta  = align->getBeta() ;
      double gamma = align->getGamma() ;
      FillRotMatrixKarimaki(EulerRotation, alpha, beta, gamma);
      
      // Combine the two factors in proper order
      HepMatrix DiscreteRotation = adet.GetDiscrete().GetRotation(); 
      HepMatrix NominalRotation = EulerRotation*DiscreteRotation;
      nominal.SetRotation(NominalRotation);  
      
      nominal.SetRotation(NominalRotation);
      adet.SetNominalFrame(nominal); 
     
    }
    
    // Close GeoDB file 
    lcReader->close();
        
  } catch (IOException& e) {
    // Print error message 
    std::cerr << e.what() << std::endl; 
  }
  
  // Print detector data sheets  
  Print();   
   
}

// 
// Write alignment data base - overwrites old DB file   
//   
void TBDetector::WriteAlignmentDB( )
{
  
  streamlog_out(MESSAGE3) << std::endl << "Write alignment DB file " << _alignmentDBFileName << std::endl << std::endl;
    
  try 
  {
    // Reopen the LCIO file this time in new mode
    LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
    
    lcWriter->open( _alignmentDBFileName, LCIO::WRITE_NEW );
    
    // Write an almost empty run header
    LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl;
    lcHeader->setRunNumber( 0 );
        
    lcWriter->writeRunHeader(lcHeader);
    delete lcHeader;
      
    LCEventImpl * event = new LCEventImpl;
    event->setRunNumber( 0 );
    event->setEventNumber( 0 );
       
    LCTime * now = new LCTime;
    event->setTimeStamp( now->timeStamp() );
    delete now;
    
    LCCollectionVec * constantsCollection = new LCCollectionVec( LCIO::LCGENERICOBJECT );
    
    for ( int ipl = 0; ipl < _numberOfSensors; ipl++ ) 
    {
      
      // Load current pixel module 
      Det & adet = GetDet(ipl);      
      
      // Stores all alignment data for this module
      AlignmentConstant * alignConst = new AlignmentConstant;
      alignConst->setSensorID( adet.GetDAQID() );
      
      // Now, write nominal tracker geometry parameters  
      HepVector Origin = adet.GetNominal().GetPosition(); 
      HepMatrix NominalRotation = adet.GetNominal().GetRotation(); 
      
      // Factor out the discrete rotation 
      HepMatrix DiscreteRotation = adet.GetDiscrete().GetRotation(); 
      HepMatrix EulerRotation = NominalRotation*DiscreteRotation.T();
       
      double alpha, beta, gamma; 
      GetAnglesKarimaki(EulerRotation, alpha, beta, gamma); 
      
      alignConst->setXOffset( Origin[0] );
      alignConst->setYOffset( Origin[1] );
      alignConst->setZOffset( Origin[2] );
      alignConst->setAlpha( alpha );
      alignConst->setBeta( beta );
      alignConst->setGamma( gamma );
       
      constantsCollection->push_back( alignConst ); 
    }                     
    event->addCollection( constantsCollection,  "alignDB" );
    lcWriter->writeEvent( event );
    delete event;
         
    lcWriter->close(); 
    streamlog_out(MESSAGE3) << std::endl << "Finish alignment DB file " << _alignmentDBFileName << std::endl << std::endl;
    
  }
  catch ( IOException& e ) 
  {
    streamlog_out ( ERROR4 ) << e.what() << std::endl
                             << "Problem creating alignmentDB. Sorry for quitting. " << std::endl;
    exit(-1);
  }
          
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

