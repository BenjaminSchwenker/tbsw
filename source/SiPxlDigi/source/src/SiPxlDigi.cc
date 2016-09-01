// History:
//- original version - ILCsoft VTXDigitizer by A. Raspereza
//- added procesor parameters: electronic noise, Z. Drasal Jan 2009
//- completely changed code structure, Z. Drasal Mar 2009
//- renamed to SiPxlGeom & created as standalone pckg, Z. Drasal Jul 2009
//- added geometry interface (SiPxlGeom), Z. Drasal Aug 2009
//- added parameter for cut on time to emulate the integration time, K. Prothmann Dec. 2009
//- new proper clustering implemented & noise generation updated, Z. Drasal Jan 2010
//- use Geant 4 steps as SimTrackerHits in input collection, B. Schwenker Jan 2010
//- no fluctuate() in ProduceSignalPoints function, B. Schwenker Jan 2010  
//- validation againts TB data (see Depfet Prag meeting 2010), B. Schwenker Jan 2010
//- added zero supression for hit pixel output;threshold is SNthr*elNoise, B. Schwenker May 2010
//- implemented sideward depletion model for charge carrier drift, B. Schwenker Sept. 2010 
//- added 'diffusive' model for charge collection into internal gates, B. Schwenker Sept. 2010 
//- added option to produce local truth hits from geant4 steps, B. Schwenker Oct. 2010
//- update of validation againts TB data (see Depfet Valencia meeting 2010), B. Schwenker Jan 2010
//- added option to specify noise seperately for all sensors, B. Schwenker Oct. 2010
//- added optional simulation of ADC, B. Schwenker Oct. 2010
//- fixed problems with function ProduceTruthHits(), B. Schwenker Nov. 2010
//- deposit photon Edep energy at step (bug fix) , B. Schwenker Feb. 2011
//- validation with S3B Cd109 source data, B. Schwenker Feb. 2011 
//- truth hits now stored in local sensor coordinates, B. Schwenker July. 2011 
//- changed decoding of zs pixels to DEPFET format, B. Schwenker Nov. 2011 
//- removed dependence on EUTELESCOPE (now DEPFETTrackTools), B. Schwenker Jan. 2012 
//- changed drift model for PXD6 sensors, B. Schwenker Feb. 2012
//- produce full analog matrix for all sensors, always, B. Schwenker May. 2013 

#include "SiPxlDigi.h"
#include "Colours.h"
#include "PhysicalConstants.h"

// Include DEPFETTrackTools 
#include "DEPFET.h" 
#include "MatrixDecoder.h"

// Include basic C
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <string>
#include <map>

// Include CLHEP classes
#include <CLHEP/Vector/Rotation.h>


// Include Gear classes
#include <gear/GEAR.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>


// Include LCIO classes
#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// Used namespaces
using namespace depfet;
using namespace CLHEP;
using namespace lcio ;
using namespace marlin ;

namespace sipxl {

//
// Instantiate this object
//
SiPxlDigi aSiPxlDigi ;



//
// Constructor
//
SiPxlDigi::SiPxlDigi() : Processor("SiPxlDigi")
{

// Processor description
   _description = "SiPxlDigi: Marlin processor creating Digits from SimTrackerHits" ;

//
// Processor parameters

   // Define compulsory parameters
   registerInputCollection(    LCIO::SIMTRACKERHIT,
                               "InputCollectionName" ,
                               "Name of input SimTrackerHit collection"  ,
                               _inColName ,
                               std::string("PXDCollection") ) ;

   registerOutputCollection(   LCIO::TRACKERDATA, 
                               "SparseDataCollectionName",
                               "Name of output sparse data collection",
                               _sparseDataCollectionName, 
                               std::string("zsdata"));

   registerProcessorParameter( "Bricked",
                               "Use bricked structure?",
                               _bricked,
                               bool(false));
   
   registerProcessorParameter( "DoublePixel",
                               "Use double pixel structure?",
                               _doublePixel,
                               bool(true));
   
   
   registerProcessorParameter( "FullMatrix",
                               "Create full analog output matrix collection?",
                               _createFullMatrix,
                               bool(false));
   
   registerProcessorParameter( "TruthHits",
                               "Create truth hit collection?",
                               _createTruthHits,
                               bool(true));
   
   registerProcessorParameter( "SiBulkDoping",
                               "Effective bulk doping concentration, in um^-3",
                               _bulkDoping,
                               double(10));
   
   registerProcessorParameter( "Uback",
                               "Back contact voltage wrt. source, in Volts",
                               _Uback,
                               double(-30));
   
   registerProcessorParameter( "Utop",
                               "Top plane voltage wrt. source, in Volts",
                               _Utop,
                               double(-5));
   
   const size_t nDetectorExample = 1;
     
   FloatVec sourceBorderExample(nDetectorExample, 7.);
   registerProcessorParameter( "SourceBorderLength",
                               "Source border length in um (one value for detector)",
                               _sourceBorderLength, 
                               sourceBorderExample );
    
   FloatVec drainBorderExample(nDetectorExample, 9.);
   registerProcessorParameter( "DrainBorderLength",
                               "Drain border length in um (one value for detector)",
                               _drainBorderLength, 
                               drainBorderExample );
   
   FloatVec clearBorderExample(nDetectorExample, 10.);
   registerProcessorParameter( "ClearBorderLength",
                               "Clear border length in um (one value for each detector)",
                               _clearBorderLength, 
                               clearBorderExample );
    
   registerProcessorParameter( "ElectronicEffects",
                               "Apply electronic effects?",
                               _electronicEffects,
                               bool(true));
   
   FloatVec elNoiseExample(nDetectorExample, 120.);
   registerProcessorParameter( "ElectronicNoise",
                               "Noise added by the electronics, set in ENC (one value for each detector)",
                               _elNoise,
                               elNoiseExample);
   
   registerProcessorParameter( "IntegrationWindow",
                               "Use integration window?",
                               _integrationWindow,
                               bool(false));
   
   registerProcessorParameter( "PoissonSmearing",
                               "Apply Poisson smearing of electrons collected on pixels?",
                               _PoissonSmearing,
                               bool(true));
   
   registerProcessorParameter( "MaxSegmentLength",
                               "Maximum distance between ionization points (in mm)",
                               _maxSegmentLength,
                               double(0.005));
   
   registerProcessorParameter( "ElectronGroupSize",
                               "Split Signalpoints in smaller groups of N electrons (in e)",
                               _eGroupSize,
                               double(100));
   
   registerProcessorParameter( "ElectronStepTime",
                               "Time step for tracking electron groups in readout plane (in ns)",
                               _eStepTime,
                               double(0.3));
   
   registerProcessorParameter( "ZSThreshold",
                               "Keep only pixel hits >= threshold in ADC units",
                               _ZSthr,
                               float(4.0) );
   
   registerProcessorParameter( "NoiseFraction",
                               "Set fraction of noise hits per readout frame",
                               _Fraction,
                               float(0.001) );
   
   registerProcessorParameter( "TanLorentz",
                               "Tangent of Lorentz angle",
                               _tanLorentzAngle,
                               double(0.0));

   registerProcessorParameter( "StartIntegration",
                               "Only Simulated hits after the StartIntegration time in ns will be digitized",
                               _startIntegration,
                               double(-10000.0));
   
   registerProcessorParameter( "StopIntegration",
                               "Only Simulated hits before the StopIntegration time in ns will be digitized",
                               _stopIntegration,
                               double(10000.0));
     
   registerProcessorParameter( "ADC",
                               "Simulate ADC?",
                               _ADC,
                               bool(false));
   
   registerProcessorParameter( "ADCRange",
                               "Set analog-to-digital converter range 0 - ? (in e)",
                               _ADCRange,
                               int(50000));
   
   registerProcessorParameter( "ADCBits",
                               "Set how many bits the ADC uses",
                               _ADCBits,
                               int(8));
   
   
}

//
// Method called at the beginning of data processing
//
void SiPxlDigi::init() {

// Random number generator - Mersenne Twister -> high quality random numbers 
   myRng = new TRandom3();

// Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;

   _currentLayerID  =  0 ;
   _currentLadderID =  0 ;
   _currentSensorID =  0 ;

   _sensorThick     = -1.;
   _sensorWidth     = -1.;
   _sensorLength    = -1.;

// Set variables in appropriate physical units
   _ADCRange             *= e;
   _bulkDoping *= 1/um/um/um; 
   _Utop *= V; 
   _Uback *= V; 
   _maxSegmentLength     *= mm;
   _eGroupSize           *= e; 
   _eStepTime            *= ns;

   _startIntegration *= ns;
   _stopIntegration  *= ns;
  

// Get geometry parameters from Gear xml file
   _geometry = new SiPxlGeom();
   _geometry->initGearParams();
   _geometry->printGearParams();

// Cross check sensor geometry data
   unsigned int numberOfSensors = 0; 
   for (short int iLayer=0; iLayer<_geometry->getNPXDLayers(); iLayer++) {
     for (short int iLadder=0; iLadder<_geometry->getNLadders(iLayer); iLadder++) {
       for (short int iSensor=0; iSensor<_geometry->getNSensors(iLayer); iSensor++) {

         numberOfSensors++;
           	          
       }
     }
   }  

   if ( numberOfSensors != _sourceBorderLength.size() ) {
     streamlog_out( WARNING3 ) << "The number of values in the sourceBorderLength vector does not match the number of sensors\n"
                               << "Resize vector consequently." << std::endl;
     _sourceBorderLength.resize(numberOfSensors, _sourceBorderLength.back());
   }

   if ( numberOfSensors != _drainBorderLength.size() ) {
     streamlog_out( WARNING3 ) << "The number of values in the drainBorderLength vector does not match the number of sensors\n"
                               << "Resize vector consequently." << std::endl;
     _drainBorderLength.resize(numberOfSensors, _drainBorderLength.back());
   }

   if ( numberOfSensors != _clearBorderLength.size() ) {
     streamlog_out( WARNING3 ) << "The number of values in the clearBorderLength vector does not match the number of sensors\n"
                               << "Resize vector consequently." << std::endl;
     _clearBorderLength.resize(numberOfSensors, _clearBorderLength.back());
   }

   if ( numberOfSensors != _elNoise.size() ) {
     streamlog_out( WARNING3 ) << "The number of values in the elNoise vector does not match the number of sensors\n"
                               << "Resize vector consequently." << std::endl;
     _elNoise.resize(numberOfSensors, _elNoise.back());
   }

   // Set border length variables to microns
   for (unsigned int iSensor = 0; iSensor < numberOfSensors; iSensor++) {
     _sourceBorderLength[iSensor] *=um;
     _drainBorderLength[iSensor] *=um;
     _clearBorderLength[iSensor] *=um;
     _elNoise[iSensor] *= e; 
   }
   
// Print set parameters
   printProcessorParams();

// CPU time start
   _timeCPU = clock()*ms/1000;

//
// ROOT variables
//
#ifdef ROOT_OUTPUT
   _rootFile = new TFile("BelleII_PXD_Hits.root","recreate");
   _rootFile->cd("");

// Declare Tree
   _rootTreeHits = new TTree("Hits","Hit info");

   _rootTreeHits->Branch("Layer"         ,&_rootLayerID     ,"Layer/I"         );
   _rootTreeHits->Branch("Ladder"        ,&_rootLadderID    ,"Ladder/I"        );
   _rootTreeHits->Branch("det"           ,&_rootSensorID    ,"det/I"           );
   _rootTreeHits->Branch("pdg"           ,&_rootPrimaryPDG  ,"pdg/I"           );
   _rootTreeHits->Branch("x_cog"         ,&_rootX           ,"x_cog/D"         );
   _rootTreeHits->Branch("y_cog"         ,&_rootY           ,"y_cog/D"         );
   _rootTreeHits->Branch("z_cog"         ,&_rootZ           ,"z_cog/D"         );
   _rootTreeHits->Branch("energy"        ,&_rootTOT         ,"energy/D" );
   _rootTreeHits->Branch("iEvt"          ,&_nEvt            ,"iEvt/I"          );   

 
// Declare SimHit Tree
   _rootTreeSimHits = new TTree("SimHits","SimHit info");
   _rootTreeSimHits->Branch("Layer"      ,&_rootLayerID         ,"Layer/I"      );
   _rootTreeSimHits->Branch("Ladder"     ,&_rootLadderID        ,"Ladder/I"     );
   _rootTreeSimHits->Branch("det"        ,&_rootSensorID        ,"det/I"        );
   _rootTreeSimHits->Branch("pdg"        ,&_rootPDG             ,"pdg/I"        );
   _rootTreeSimHits->Branch("length"     ,&_rootTrackLength     ,"length/D"     );
   _rootTreeSimHits->Branch("momentum"   ,&_rootMomentum        ,"momentum/D"   );
   _rootTreeSimHits->Branch("dEdx"       ,&_rootdEdx            ,"dEdx/D"       );
   _rootTreeSimHits->Branch("iEvt"       ,&_nEvt                ,"iEvt/I"       );
   _rootTreeSimHits->Branch("x_pos"      ,&_rootXPos            ,"x_pos/D"         );
   _rootTreeSimHits->Branch("y_pos"      ,&_rootYPos            ,"y_pos/D"         );
   _rootTreeSimHits->Branch("z_pos"      ,&_rootZPos            ,"z_pos/D"         );
 
   
// Declare charge collection time histogram, in ns
   _rootChargeCollectionTime = new TH1D("collectionTime", "Charge Collection Time", 300, 0, 300);
   
// Declare number of random walk steps histogram
   _rootRandomWalkSteps = new TH1D("steps", "Random walk steps", 200, 0, 200);
   
// Declare killed random walk histogram 
   _rootKilledRandomWalks = new TH1D("killed", "Killed random walks", 2, 0, 2);
   
#endif

}

//
// Method called for each run
//
void SiPxlDigi::processRunHeader(LCRunHeader * run)
{

// Print run number
   streamlog_out(MESSAGE3) << DGREEN
                           << " Processing run: "
                           << ENDCOLOR
                           << (run->getRunNumber())
                           << std::endl << std::endl;

   _nRun++ ;
}

//
// Method called for each event
//
void SiPxlDigi::processEvent(LCEvent * evt)
{
	//local error count variable
	static int hits_cutted = 0;

// Print event number
   if ((evt->getEventNumber()+1)%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                   << (evt->getEventNumber()+1)
                                                                   << std::endl << std::endl;

   //streamlog_out(MESSAGE3) << "Events processed: " << (evt->getEventNumber()+1) << std::endl << std::endl;


//
// Open collections
   try {

      // Open SimTrackerHit collection
      LCCollection * STHcol = evt->getCollection( _inColName );

      // Set collection decoder
      CellIDDecoder<SimTrackerHit> cellIDDec(STHcol);

      // Number of SimTrackerHits
      int nSTH = STHcol->getNumberOfElements();

      // Initialize PXD map of all digits -> needed for proper clustering
      DigitsMap digitsMap;

      //
      // Loop over SimTracker hits;
      streamlog_out(MESSAGE2) << " Producing all digits ..." << std::endl;

      for (int i=0; i<nSTH; ++i) {

         SimTrackerHit * simTrkHit = dynamic_cast<SimTrackerHit*>(STHcol->getElementAt(i));

         // Set current - layer ID, ladder ID and sensor ID
         _currentLayerID   = _geometry->getLayerIDCTypeNo(cellIDDec(simTrkHit)["layer"]);
         _currentLadderID  = cellIDDec(simTrkHit)["ladder"];
         _currentSensorID  = cellIDDec(simTrkHit)["sensor"];


         // Cut on simHit creation time --> simulate integration time of a sensor (if option switched on))
	 if ((simTrkHit != 0) && (_integrationWindow)) {
         
	   if (simTrkHit->getTime()*ns < _startIntegration || simTrkHit->getTime()*ns > _stopIntegration) {
	     hits_cutted++;
	     continue;
	   }
				
	}


         // Set current - sensor thickness, width, length
         _sensorThick      = _geometry->getSensorThick( _currentLayerID, _currentSensorID);
         _sensorWidth      = _geometry->getSensorWidth( _currentLayerID, _currentSensorID);
         _sensorLength     = _geometry->getSensorLength(_currentLayerID, _currentSensorID);

         // Print
         streamlog_out(MESSAGE2) << " LayerID: "               << _currentLayerID
                                 << ", LadderID: "             << _currentLadderID
                                 << ", SensorID: "             << _currentSensorID << std::endl;
         streamlog_out(MESSAGE2) << std::setiosflags(std::ios::fixed | std::ios::internal )
                                 << std::setprecision(3)
                                 << " Simulated hit global[mm]: (" << simTrkHit->getPosition()[0]
                                 << ", "                           << simTrkHit->getPosition()[1]
                                 << ", "                           << simTrkHit->getPosition()[2] << ")"
                                 << std::setprecision(0)
                                 << std::endl;

         // Check if digitizing pixels
         if (_geometry->getLayerType(_currentLayerID) != pixel ) {

            streamlog_out(ERROR) << "SiPxlDigi::processEvent - sensor is not of pixel type!!!" << std::endl;
            exit(1);
         }

         //
         // Produce ionisation points along the track
         IonisationPointVec ionisationPoints;
         ProduceIonisationPoints(simTrkHit, ionisationPoints);

         //
         // Produce signal points
         streamlog_out(MESSAGE1) << " Producing signal points ..." << std::endl;
         SignalPointVec signalPoints;
         ProduceSignalPoints(ionisationPoints, signalPoints);

         // Release memory - clear ionisationPoints
         IonisationPointVec::iterator iterIPVec;

         for (iterIPVec=ionisationPoints.begin(); iterIPVec!=ionisationPoints.end(); iterIPVec++) {

            IonisationPoint * iPoint = (*iterIPVec);
            delete iPoint;
            iPoint = 0;
         }
         ionisationPoints.clear();

         //
         // Produce digits
         streamlog_out(MESSAGE1) << " Producing digits ..." << std::endl;
         DigitVec digits;
         ProduceDigits(signalPoints, digits);

         // Release memory - clear signalPoints
         SignalPointVec::iterator iterSPVec;

         for (iterSPVec=signalPoints.begin(); iterSPVec!=signalPoints.end(); iterSPVec++) {

            SignalPoint * sPoint = (*iterSPVec);
            delete sPoint;
            sPoint = 0;
         }
         signalPoints.clear();

         //
         // Update PXD map (of all digits) with new digits -> necessary to perform correct clustering
         UpdateDigitsMap(digitsMap, digits);

         // Clear vector of digits (all digits will be deleted after clustering
         digits.clear();

	
      } // Loop over simTrkHits

      if ( _createFullMatrix  ) {
	
	// 
        // Produce calibrated full frame data, mimics TB data  
	streamlog_out(MESSAGE1) << " Writing Full Matrix ..." << std::endl;

	// Create TrackerData collection (full frame Digits)
      	LCCollectionVec * TMatrixCol = new LCCollectionVec(LCIO::TRACKERDATA);

	// create full anolog output for all PXD sensor matrices
	ProduceFullMatrix(TMatrixCol, digitsMap); 

	// strore stuff in LCIO file 
	evt->addCollection(TMatrixCol,"FullAnalogPXD");
	
      } 

     
      if (_electronicEffects) {

        // Produce noise digits and update PXD map
        streamlog_out(MESSAGE2) << " Producing noisy digits ..." << std::endl;
	
	// add electronic noise to truth digits
	ProduceNoiseEffects(digitsMap);
 
	// create pure noise digits
	ProduceNoiseDigits(digitsMap);

      }
      
      // 
      // Produce sparse pixel output collection  
      streamlog_out(MESSAGE1) << " Writing Sparse Pixels ..." << std::endl;
	
      // Create sparse data collection
      LCCollectionVec * sparseDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);       
        
      ProduceSparsePixels( sparseDataCollection  , digitsMap); 
	
      // Store stuff in LCIO file 
      evt->addCollection(sparseDataCollection, _sparseDataCollectionName);
        
      
      if ( _createTruthHits  ) {

        // Produce truth hits from simtrakerhits
	streamlog_out(MESSAGE1) << " Writing Truth Hits ..." << std::endl;
        
	// Create truth hits collection
      	LCCollectionVec * TTruthHitCol = new LCCollectionVec(LCIO::TRACKERHIT);   

	// Fill truth hits collection from simtrackerhits 
	ProduceTruthHits(TTruthHitCol, STHcol ); 
        
	// Strore stuff in LCIO file 
	evt->addCollection(TTruthHitCol,"TruthHit");     
      }

      // Release memory - clear digits
      DigitsMap::const_iterator iterDigitsMap;
      DigitMap::const_iterator iterSensorMap;

      for (iterDigitsMap=digitsMap.begin(); iterDigitsMap!=digitsMap.end(); iterDigitsMap++) {

         // Denote PXD sensor map as
         DigitMap digitsSensorMap = iterDigitsMap->second;

         for (iterSensorMap=digitsSensorMap.begin(); iterSensorMap!=digitsSensorMap.end(); iterSensorMap++) {

            Digit * digit = iterSensorMap->second;

            if (digit!=0) delete digit;
            digit = 0;
         }
         digitsSensorMap.clear();
      }
      digitsMap.clear();


   }
   catch(DataNotAvailableException &e){}
   if(hits_cutted != 0) {
	   streamlog_out(DEBUG4) << "SimTrackerHit outside " << hits_cutted << " times outside of integration time of sensor" << std::endl;
   }
   hits_cutted = 0;

   _nEvt ++ ;
}

//
// Method called after each event to check the data processed
//
void SiPxlDigi::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void SiPxlDigi::end()
{

   // CPU time end
   _timeCPU = clock()*ms/1000 - _timeCPU;

   // clean up
   delete  myRng; 

   // Print message
   streamlog_out(MESSAGE3) << std::endl
                           << " "
                           << "Time per event: "
                           << std::setiosflags(std::ios::fixed | std::ios::internal )
                           << std::setprecision(3)
                           << _timeCPU/_nEvt/ms
                           << " ms"
                           << std::endl
                           << std::setprecision(3)
                           << std::endl
                           << " "
                           << "Processor succesfully finished!"
                           << std::endl;

  


#ifdef ROOT_OUTPUT

   // Close file
   _rootFile->cd("");
   _rootFile->Write();
   _rootFile->Close();

#endif
}

//
// Method transforming hit to local coordinates
//
void SiPxlDigi::TransformToLocal(const SimTrackerHit * simTrkHit, SpacePoint & hitLocal)
{
   Hep3Vector position(simTrkHit->getPosition()[0]*mm ,simTrkHit->getPosition()[1]*mm ,simTrkHit->getPosition()[2]*mm );
   Hep3Vector momentum(simTrkHit->getMomentum()[0]*GeV,simTrkHit->getMomentum()[1]*GeV,simTrkHit->getMomentum()[2]*GeV);

   position = _geometry->transformPointToLocal(_currentLayerID, _currentLadderID, _currentSensorID, position);
   momentum = _geometry->transformVecToLocal(  _currentLayerID, _currentLadderID, _currentSensorID, momentum);

   // Save final results
   hitLocal.position  = position;
  
   if ( momentum.mag() != 0 ) {

     hitLocal.direction = momentum/momentum.mag();

   } else {
     
     hitLocal.direction = momentum; 
   }
}

//
// Method transforming hit to global coordinates
//
void SiPxlDigi::TransformToGlobal(const SpacePoint & hitLocal, SpacePoint & hitGlobal)
{
   Hep3Vector position(hitLocal.position);
   Hep3Vector direction(hitLocal.direction);

   position  = _geometry->transformPointToGlobal(_currentLayerID, _currentLadderID, _currentSensorID, position );
   direction = _geometry->transformVecToGlobal(  _currentLayerID, _currentLadderID, _currentSensorID, direction);

   // Save final results
   hitGlobal.position  = position;
   hitGlobal.direction = direction;

}

//
// Method producing ionisation points along the track
//
void SiPxlDigi::ProduceIonisationPoints(const SimTrackerHit * simTrkHit, IonisationPointVec & ionisationPoints)
{
    
   
   // Space point
   SpacePoint hitLocal;
   
   // Transform geant4 step to local sensor coordinates
   TransformToLocal(simTrkHit, hitLocal);

   /* delete benni 
   if (_currentSensorID == 2) {

     streamlog_out(MESSAGE3) << std::setiosflags(std::ios::fixed | std::ios::internal )
                           << std::setprecision(4)
                           << " Simulated hit global[mm]: (" << hitLocal.position.getX()
                           << ", "                           << hitLocal.position.getY()
                           << ", "                           << hitLocal.position.getZ() << ")"
                           << std::resetiosflags(std::ios::showpos)
                           << std::setprecision(0)
                           << std::endl << std::endl;
     
     Hep3Vector pos(simTrkHit->getPosition()[0]*mm ,simTrkHit->getPosition()[1]*mm ,simTrkHit->getPosition()[2]*mm );
     streamlog_out(MESSAGE3) << std::setiosflags(std::ios::fixed | std::ios::internal )
                           << std::setprecision(4)
                           << " Simulated hit local[mm]:  (" << pos.getX()
                           << ", "                           << pos.getY()
                           << ", "                           << pos.getZ() << ")"
                           << std::resetiosflags(std::ios::showpos)
                           << std::setprecision(0)
                           << std::endl << std::endl;
     
   }
   */
   
   streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                           << std::setprecision(3)
                           << " Simulated hit local[mm]:  (" << hitLocal.position.getX()
                           << ", "                           << hitLocal.position.getY()
                           << ", "                           << hitLocal.position.getZ() << ")"
                           << std::resetiosflags(std::ios::showpos)
                           << std::setprecision(0)
                           << std::endl << std::endl;
   
   
   // Print info
   streamlog_out(MESSAGE1) << " Producing ionisation points ..." << std::endl;
   
   // Calculate entryPoint and exitPoint (in mm)
   double trackLength  = simTrkHit->getPathLength()*mm;
   
   // Energy deposit along G4 step 
   double Edep = (simTrkHit->getdEdx()*GeV);

   Hep3Vector entryPoint = hitLocal.position - hitLocal.direction*trackLength/2.;
   Hep3Vector exitPoint  = hitLocal.position + hitLocal.direction*trackLength/2.;

   // Get particle ID 
   int mcPDG = simTrkHit->getMCParticle()->getPDG(); 
   
   // Calculate ionisation points in sensitive silicon 
   // depending in particle type     
   if ( mcPDG == 22 ) {
       
     // For photons: energy deposit due to photoelectric or 
     // compton effect. Energy deposit is atom relaxation 
     // energy at end of photon step.
     
     // Create single ionisation point at end of step  
     ionisationPoints.resize(1);
     
     IonisationPoint * iPoint = new IonisationPoint;
       
     iPoint->position = exitPoint; 
     iPoint->eLoss    = Edep; 
     ionisationPoints[0] = iPoint;
     
     // Check if iPoint is within sensor boundaries
     if (_geometry->isPointOutOfSensor(_currentLayerID, _currentSensorID, iPoint->position)) {
          
       streamlog_out(ERROR) << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(3)
                              << "SiPxlDigi::ProduceIonisationPoints - ionPoint: " << iPoint->position/mm << " out of sensor!!!"
                              << std::setprecision(0) << std::endl
                              << "SensorID is " << _currentSensorID
                              << std::endl;
       exit(1);
     }
     
     // Print
     streamlog_out(MESSAGE1) << "  Hit local ionPoints (photoeffect): " << std::endl;
     streamlog_out(MESSAGE2) << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(3)
                              << "   Pos [mm]: ( " << iPoint->position.getX()/mm << ", " << iPoint->position.getY()/mm << ", " << iPoint->position.getZ()/mm << " )"
                              << " , dE [keV]: "   << iPoint->eLoss/keV
                              << std::setprecision(0)
                              << std::endl;
     
   } else {
     
     // For charged particles: energy is smeared along the step in 
     // small segments. 
     
     int numberOfSegments = int(trackLength/_maxSegmentLength) + 1;
     
     // Calculate mean energy loss in each segment
     double dEMean = Edep/((double)numberOfSegments);
     
     ionisationPoints.resize(numberOfSegments);
          
     // Set ionisation points
     for (int i=0; i<numberOfSegments; i++) {
       
       IonisationPoint * iPoint = new IonisationPoint;
       
       iPoint->position = entryPoint + hitLocal.direction*trackLength/numberOfSegments*(i+0.5);
       iPoint->eLoss    = dEMean; 
       ionisationPoints[i] = iPoint;
       
       // Check if iPoint is within sensor boundaries
       if (_geometry->isPointOutOfSensor(_currentLayerID, _currentSensorID, iPoint->position)) {
         
         streamlog_out(ERROR) << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(3)
                              << "SiPxlDigi::ProduceIonisationPoints - ionPoint: " << iPoint->position/mm << " out of sensor!!!"
                              << std::setprecision(0) << std::endl
                              << "SensorID is " << _currentSensorID
                              << std::endl;
         exit(1);
       }
       
       // Print
       streamlog_out(MESSAGE1) << "  Hit local ionPoints (ionisation): " << std::endl;
       streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(3)
                              << "   Pos [mm]: ( " << iPoint->position.getX()/mm << ", " << iPoint->position.getY()/mm << ", " << iPoint->position.getZ()/mm << " )"
                              << " , dE [keV]: "   << iPoint->eLoss/keV
                              << std::setprecision(0)
                              << std::endl;
     
     } // end segment loop 
     
   } // endif particle type 
   
}

//
// Method producing signal points on each sensor
//
void SiPxlDigi::ProduceSignalPoints(const IonisationPointVec & ionisationPoints, SignalPointVec & signalPoints)
{
   // Calculate number of ionisation points
   int numberOfIonPoints = ionisationPoints.size();

   // Run over all ionisation points and create signal points
   signalPoints.clear();
   signalPoints.resize(numberOfIonPoints);
   
   // Print
   streamlog_out(MESSAGE1) << "  Hit local signalPoints: " << std::endl;
   
   for (int i=0; i<numberOfIonPoints; i++) {
      
      // Get current ionisation point
      IonisationPoint * iPoint = ionisationPoints[i];
      
      // Final point
      Hep3Vector finalPos(iPoint->position);
      
      //  Charge cloud created at distance to top plane 
      double d = _sensorThick - iPoint->position.getX();
      d = sqrt( d * d );
      
      //std::cout << "Depth of charge deposit [mm]: " << d/mm << std::endl;
      
      // Potential valley is at distance to top plane
      double x0 = _sensorThick/2 + Perm_Si*(_Uback-_Utop)/(e*_bulkDoping*_sensorThick);
      double dx = 0.003;  
      
      //std::cout << "Depth of potential vallay [mm]: " << x0/mm << std::endl;
      
      // Drift time into potential valley -  
      double td = 1.4*1.4*um*um/( 2* Utherm * e_mobility);
        
      if ( d > x0 + dx ) {
        td += Perm_Si * ::log( (d - x0)/dx ) / ( e_mobility * e * _bulkDoping);
      } else if ( d < x0 - dx )  {
        td += Perm_Si * ::log( (x0 - d)/dx ) / ( e_mobility * e * _bulkDoping);
      }
      
      //std::cout << "Vertical drift time [ns]: " << td/ns << std::endl; 
      
      // Diffusive spread in lateral plane
      double sigmaDiffus = sqrt( 2 * Utherm * e_mobility * td ); 
      
      //std::cout << "Diffusive spread [mm]: " << sigmaDiffus/mm << std::endl; 
      
      double sigmaY = sigmaDiffus;
      double sigmaZ = sigmaDiffus; 
      
      //  After Lorentz shift
      double onPlaneY = iPoint->position.getY() + _tanLorentzAngle * (d - x0);
      double onPlaneZ = iPoint->position.getZ();
      
      finalPos.setY(onPlaneY);
      finalPos.setZ(onPlaneZ);
      finalPos.setX(x0);
       
      // Save info in signal point
      SignalPoint * sPoint = new SignalPoint;
      
      sPoint->position.setX(finalPos.getX());
      sPoint->position.setY(finalPos.getY());
      sPoint->position.setZ(finalPos.getZ());
      
      sPoint->sigma.setX(_sensorThick-x0);
      sPoint->sigma.setY(sigmaY);
      sPoint->sigma.setZ(sigmaZ);
      
      // Charge in electrons
      sPoint->charge = iPoint->eLoss/Eeh * e;
         
      // Check if sinal point not out of sensor - if yes set border of the sensor
      _geometry->correctPointOutOfSensor(_currentLayerID, _currentSensorID, sPoint->position);
      
      // Save signal point
      signalPoints[i] = sPoint;
       
      // Print
      streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(3)
                              << "   Pos (X/Y/Z) [mm]:( " << sPoint->position.getX()/mm << ", " << sPoint->position.getY()/mm << ", " << sPoint->position.getZ()/mm << " )"
                              << " , q [fC]: "    << sPoint->charge/fC
                              << std::setprecision(0)
                              << std::endl;
       
      streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(3)
                              << "   Sigma (X/Y/Z) [um]:( " << sPoint->sigma.getX()/um << ", " << sPoint->sigma.getY()/um << ", " << sPoint->sigma.getZ()/um << " )"
                              << std::setprecision(0)
                              << std::endl;
      
   }
   
   
}

//
// Method producing digits from signal points
//
void SiPxlDigi::ProduceDigits(const SignalPointVec & signalPoints, DigitVec & digits)
{
   
   // Calculate number of signal points
   int numberOfSigPoints = signalPoints.size(); 
   
   // Run over all signal points and create digits
   digits.clear();
   
   // Print
   streamlog_out(MESSAGE1) << "Producing hit digits ... " << std::endl;
   
   for (int i=0; i<numberOfSigPoints; ++i) {
      
      // Get current signal point
      SignalPoint * sPoint = signalPoints[i];
      
      // Calculate centre of gaussian charge cloud
      double centreZ    = sPoint->position.getZ();
      double centreRPhi = sPoint->position.getY();
      
      // Calculate width of charge cloud
      double sigmaZ     = sPoint->sigma.getZ();
      double sigmaRPhi  = sPoint->sigma.getY();
      
      // Get number of electrons in cloud
      double clusterCharge = sPoint->charge;
      
      streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(3) << std::endl;
      
      streamlog_out(MESSAGE1) << "Signal point on sensor " << _currentSensorID <<  std::endl;
      streamlog_out(MESSAGE1) << "point charge [e-] " << clusterCharge  <<  std::endl; 
      streamlog_out(MESSAGE1) << "point position [mm] " <<  centreZ/mm << ", " << centreRPhi/mm <<  std::endl;
      streamlog_out(MESSAGE1) << "point sigma [um]    " << sigmaZ/um << ", " << sigmaRPhi/um <<  std::endl;  
      
      
      // Now, the signal point is split into groups of electons spread around the 
      // signal point center. Each group is tracked into an internal gate.
      int numberOfGroups = int(clusterCharge/_eGroupSize) + 1;
      double groupCharge = clusterCharge/((double)numberOfGroups);
      
      // Random walk step length 
      double sigmaDiffus = sqrt( 2*Utherm*e_mobility*_eStepTime ) ;
      
      // N.B. diffusion step length must be smaller than any border length (drain/source/clear)
      // to simulate random walk correctly  
      //std::cout << "sigmaDiffus [um]: " <<    sigmaDiffus/um << std::endl;    
      
      
      for (int iGroup=0; iGroup<numberOfGroups; ++iGroup)  {         
        
        // Initial group position   
        double groupPosZ = centreZ + myRng->Gaus(0, sigmaZ);
        double groupPosRPhi = centreRPhi + myRng->Gaus(0, sigmaRPhi);
            
        // Pixel cellID in Z 
        int iZ = _geometry->getPixelIDInZ(_currentLayerID, _currentSensorID, groupPosZ);
        double pixelPosZ = _geometry->getPixelPosInZ(_currentLayerID, _currentSensorID, iZ);
            
        // Pixel cellID in RPhi
        int iRPhi = _geometry->getPixelIDInRPhi(_currentLayerID, _currentSensorID, _bricked, groupPosRPhi, groupPosZ);
        double pixelPosRPhi = _geometry->getPixelPosInRPhi(_currentLayerID, _currentSensorID, _bricked, iRPhi, iZ); 
           
        // control variables 
        double collectionTime = 0;
        bool insideIG = false; 
        
        for (int iStep = 0; iStep < 400; ++iStep) {
          
          //std::cout << "  iStep " <<  iStep << " at time (ns) " << collectionTime/ns <<std::endl;
          
          // Calculate border of internal gate region
          double pitchZ    = _geometry->getSensorPitchInZ(_currentLayerID, _currentSensorID);
          double pitchRPhi = _geometry->getSensorPitchInRPhi(_currentLayerID, _currentSensorID);           
          double lowerZ      = pitchZ/2. - _drainBorderLength[_currentSensorID]; 
          double upperZ      = pitchZ/2. - _sourceBorderLength[_currentSensorID];          
          double halfwidth   = pitchRPhi/2. - _clearBorderLength[_currentSensorID];
          double deltaZ = groupPosZ - pixelPosZ; 
          double deltaRPhi = groupPosRPhi - pixelPosRPhi;  
          
          // Check if sensor uses double pixel structure
          if ( _doublePixel ) {
            if(iZ % 2 != 0)  {
	      // N.B. odd pixel rows have source at lower edge
              lowerZ  = pitchZ/2. - _sourceBorderLength[_currentSensorID]; 
              upperZ  = pitchZ/2. - _drainBorderLength[_currentSensorID]; 
            } 
          }
          // debug
          //std::cout << "  SensorID " << _currentSensorID <<std::endl;
          //std::cout << "  in pixel cellID: " << iZ << ", " << iRPhi  <<std::endl;
          //std::cout << "  pix  centre (mm): " << pixelPosZ/mm << ", " << pixelPosRPhi/mm << std::endl;
          //std::cout << "  group pos   (mm): " << groupPosZ/mm << ", " << groupPosRPhi/mm << std::endl;
          //std::cout << "  delta       (mm): " << deltaZ/mm    << ", " << deltaRPhi/mm    << std::endl;
          
          // Check if charge is inside IG region
          if (   deltaRPhi > - halfwidth &&  deltaRPhi <  halfwidth &&
                 deltaZ    > - lowerZ   &&  deltaZ    < upperZ        )
          {
            insideIG = true;
            //std::cout  << " inside!! " << std::endl;
            #ifdef ROOT_OUTPUT
            _rootRandomWalkSteps->Fill(iStep);
            #endif 
            break;     
          }
            
  	  // Update position of group, possibly leaving current pixel 
          collectionTime += _eStepTime;
          groupPosZ += myRng->Gaus(0, sigmaDiffus);
          groupPosRPhi += myRng->Gaus(0, sigmaDiffus); 
              
          // Update current cellID in Z  
          iZ = _geometry->getPixelIDInZ(_currentLayerID, _currentSensorID, groupPosZ);
          pixelPosZ = _geometry->getPixelPosInZ(_currentLayerID, _currentSensorID, iZ);
                       
          // Update current cellID in RPHI  
          iRPhi = _geometry->getPixelIDInRPhi(_currentLayerID, _currentSensorID, _bricked, groupPosRPhi, groupPosZ);
          pixelPosRPhi = _geometry->getPixelPosInRPhi(_currentLayerID, _currentSensorID, _bricked, iRPhi, iZ);
               
        } // end of group tracking
         
        #ifdef ROOT_OUTPUT 
        _rootChargeCollectionTime->Fill(collectionTime/ns);  
        if ( !insideIG ) _rootKilledRandomWalks->Fill(0); 
  	else _rootKilledRandomWalks->Fill(1); 
        #endif
        
        // Now, pixel cell [iZ,iRPhi] collects group of electrons         
        // Update digit info if exists, otherwise create new one
        DigitVec::iterator iterDVec;
        bool digitFound = false;
               
        // Update digit info
        for (iterDVec=digits.begin(); iterDVec!=digits.end(); iterDVec++) {
              
          Digit * digit = (*iterDVec);
                 
          if ( (digit->cellIDRPhi == iRPhi) && (digit->cellIDZ == iZ) ) {
                   
            digit->charge += groupCharge;
            digitFound = true;
            break;
          }
        }
                 
        // Such a digit not found - create new digit
        if (digitFound == false) {
               
          Digit * digit = new Digit;
                
          digit->cellIDZ    = iZ;
          digit->cellIDRPhi = iRPhi;
          digit->cellPosZ   = pixelPosZ;
          digit->cellPosRPhi= pixelPosRPhi;
          digit->charge     = groupCharge;                
          digits.push_back(digit);
        } 
         
     } // For groups  
  
   } // For signalPoints
   
}

//
// Method producing full analog output matrix from all digits (noise + signal)
//
void SiPxlDigi::ProduceFullMatrix(LCCollectionVec * matrix , const DigitsMap & digitsMap)
{
  
  // Prepare an encoder also for the matrix output
  CellIDEncoder<TrackerDataImpl> idMatrixEncoder(DEPFET::MATRIXDEFAULTENCODING, matrix);
   
  // Prepare navigation in digits map 
  DigitsMap::const_iterator iterDigitsMap;
  DigitMap::const_iterator iterSensorMap;
  
  // Go through all sensors & generate digits
  for (short int iLayer=0; iLayer<_geometry->getNPXDLayers(); iLayer++) {
    for (short int iLadder=0; iLadder<_geometry->getNLadders(iLayer); iLadder++) {
      for (short int iSensor=0; iSensor<_geometry->getNSensors(iLayer); iSensor++) {
         
        // Describe sensor by unique ID
        int uniqSensorID = _geometry->encodeSensorID(iLayer, iLadder, iSensor);
         
        // Now that we know which is the sensorID, we can ask to GEAR
        // about the number of pixels.
        // note: X == RPhi and Y == Z
        
        int noOfXPixels = _geometry->getSensorNPixelsInRPhi(iLayer,iSensor); 
        int noOfYPixels = _geometry->getSensorNPixelsInZ(iLayer,iSensor); 
         
        // Internal matrix decoding, must be sync. to reconstruction algorithms!!
        MatrixDecoder matrixDecoder( noOfXPixels, noOfYPixels); 
        
        // Store full analog frame 
        std::vector<float> chargeVec( noOfXPixels*noOfYPixels , 0.0  );     
        
        // Check if there is signal       
        iterDigitsMap = digitsMap.find(uniqSensorID);
        
        // IF COPY SIGNAL   
        if (iterDigitsMap!=digitsMap.end()) {
          
          // These are the signal digits on current sensor
          DigitMap digitsSensorMap = iterDigitsMap->second;
          
          // Copy digits into local table
          for (iterSensorMap=digitsSensorMap.begin(); iterSensorMap!=digitsSensorMap.end(); iterSensorMap++) {
                       
            Digit * currentDigit = iterSensorMap->second;
            
            // At this point, digits contain only true signal charges.
            // That are electrons from ionisations processes collected 
            // in the internal gate a DEPFET pixel cell. 
            double chargeValue = currentDigit->charge;
            int xPixel = currentDigit->cellIDRPhi; 
            int yPixel = currentDigit->cellIDZ; 
            
            // Add Poisson smearing - only for ionisation charge 
            if (_PoissonSmearing) {
                           
              // For big charge assume Gaussian distr.
              if (chargeValue > (1000.*e)) {   
                double sigma  = sqrt(chargeValue);
                chargeValue = myRng->Gaus(chargeValue,sigma);
              }
              // Otherwise Poisson distr.
              else {
                chargeValue = myRng->Poisson(chargeValue);
              }
            
            } // Poisson smearing
             
            // Internal matrix decoding 
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  
            chargeVec[iPixel] = chargeValue;
            
          }
          
        } // ENDIF COPY SIGNAL
        
        // IF ADD ELETRONIC NOISE
        if (_electronicEffects) { 
          
          double elNoise = _elNoise[iSensor]; 
          
          // Transistor and readout noise + ADC 
          for ( int iPixel = 0; iPixel < (int) chargeVec.size(); iPixel++ ) {
               
            // Assume gaussian readout noise   
            double noise   = myRng->Gaus(0.,elNoise);
            chargeVec[iPixel] += noise;
            
            // now, simulate ADC
            if (_ADC) chargeVec[iPixel] = (double) getInADCUnits(chargeVec[iPixel]);
        
          }
        } // ENDIF ELETRONIC NOISE
        
        // Store final results
        TrackerDataImpl  * sensorMatrix = new TrackerDataImpl;
        sensorMatrix->setChargeValues(chargeVec);
           
        idMatrixEncoder["sensorID"] = iSensor;
        idMatrixEncoder["xMin"]     = 0;
        idMatrixEncoder["xMax"]     = noOfXPixels-1;
        idMatrixEncoder["yMin"]     = 0;
        idMatrixEncoder["yMax"]     = noOfYPixels-1;
        idMatrixEncoder.setCellID(sensorMatrix);
        matrix->push_back(sensorMatrix);
           
      }
    }
  } // 3x SENSOR LOOP 
          
}

//
// Method producing sparsified pixels output from all digits (noise + signal)
//
void SiPxlDigi::ProduceSparsePixels(LCCollectionVec * sparseDataCollection , const DigitsMap & digitsMap)
{
  
  // Print
  streamlog_out(MESSAGE2) << " Writing zero suppressed data " << std::endl;
  
  // Prepare an encoder for zero suppressed pixels
  CellIDEncoder<TrackerDataImpl> sparseDataEncoder( "sensorID:6,sparsePixelType:5" , sparseDataCollection ); 
  
  // Prepare navigation in digits map 
  DigitsMap::const_iterator iterDigitsMap;
  DigitMap::const_iterator iterSensorMap;
  
  // Go through all sensors
  for (iterDigitsMap=digitsMap.begin(); iterDigitsMap!=digitsMap.end(); iterDigitsMap++) {
    
    // Get layerID, ladderID & sensorID corresponding to current sensor
    _geometry->decodeSensorID(_currentLayerID, _currentLadderID, _currentSensorID, iterDigitsMap->first);
    
    // Now that we know which is the sensorID, we can ask to GEAR
    // which are the minX, minY, maxX and maxY.
    // note: X == RPhi and Y == Z
      
    
    
    // reset the pixel counter for the PixelID
    int pixelCounter = 0;
    
    // Denote PXD sensor map as
    DigitMap digitsSensorMap = iterDigitsMap->second;
     
    // Prepare a TrackerData to store zs pixels
    TrackerDataImpl* zspixels = new TrackerDataImpl;
     
    // Set description for zspixels 
    sparseDataEncoder["sensorID"] = _currentSensorID;
    sparseDataEncoder["sparsePixelType"] = 0;
    sparseDataEncoder.setCellID( zspixels );
    
    // Loop over all digits on this sensor 
    for (iterSensorMap=digitsSensorMap.begin(); iterSensorMap!=digitsSensorMap.end(); iterSensorMap++) {
       
      Digit * currentDigit = iterSensorMap->second;
       
      double signal  = currentDigit->charge;
       
      // Note: zero suppressed pixels have signal == 0, just skip 
      if (signal == 0) continue;  
      
      int yPixel = currentDigit->cellIDZ;
      int xPixel = currentDigit->cellIDRPhi;
      
      // Store pixel data int DEPFET cluster format 
      zspixels->chargeValues().push_back( xPixel );
      zspixels->chargeValues().push_back( yPixel );
      zspixels->chargeValues().push_back( signal );    
      
      ++pixelCounter;
        
      // Print detailed pixel summary 
      streamlog_out(MESSAGE2) << "Found pixel Nr. " << pixelCounter << " on sensor " << _currentSensorID << std::endl;  
      streamlog_out(MESSAGE2) << "   x:" << xPixel << ", y:" << yPixel << ", charge:" << signal << std::endl;
      
    
    }
    
    // Add zs pixels to collection 
    sparseDataCollection->push_back( zspixels );
        
  } // end for: loop over sensors
  
}

// 
// Method producing truth hits from simhits
//  
void SiPxlDigi::ProduceTruthHits( LCCollectionVec * truthhit, LCCollection * simTrkHits)
{
   
   
   // Number of SimTrackerHits
   int nHits = simTrkHits->getNumberOfElements();
   
   // Main loop over simhits == geant4 steps 
   // Create a tracker hit for each particle entering an active sensor
   for (int iHit=0; iHit<nHits; ++iHit) {   

     // Set collection decoder
     CellIDDecoder<SimTrackerHit> stepDecoder(simTrkHits );
         
     SimTrackerHit * currentSimTrkHit = dynamic_cast<SimTrackerHit*>(simTrkHits->getElementAt(iHit));
     
     // Set current - layer ID, ladder ID and sensor ID
     _currentLayerID   = _geometry->getLayerIDCTypeNo(stepDecoder(currentSimTrkHit)["layer"]);
     _currentLadderID  = stepDecoder(currentSimTrkHit)["ladder"];
     _currentSensorID  = stepDecoder(currentSimTrkHit)["sensor"];       
     
     
     #ifdef ROOT_OUTPUT
     
     // fill simhit tree     
     _rootLayerID   = _geometry->getLayerIDCTypeNo(stepDecoder(currentSimTrkHit)["layer"]);
     _rootLadderID  = stepDecoder(currentSimTrkHit)["ladder"];
     _rootSensorID  = stepDecoder(currentSimTrkHit)["sensor"];
     
     Hep3Vector partMomentum(currentSimTrkHit->getMomentum()[0]*GeV, currentSimTrkHit->getMomentum()[1]*GeV, currentSimTrkHit->getMomentum()[2]*GeV);
     double partMomentumMag = partMomentum.mag();
     
     _rootPDG = currentSimTrkHit->getMCParticle()->getPDG() ; 
     _rootTrackLength = currentSimTrkHit->getPathLength()/um;
     _rootMomentum = partMomentumMag/keV;
     _rootdEdx = currentSimTrkHit->getdEdx()*GeV/keV;
     
      
     _rootXPos =  currentSimTrkHit->getPosition()[0]; 
     _rootYPos = currentSimTrkHit->getPosition()[1];
     _rootZPos = currentSimTrkHit->getPosition()[2];
       
     // Fill simhit tree
     _rootFile->cd("");
     _rootTreeSimHits->Fill();  
     
     #endif
     
     // Print - for debugging only 
     streamlog_out(MESSAGE2) << std::endl
                             << "SimHit on LayerID: "     << _currentLayerID
                             << ", LadderID: "            << _currentLadderID
                             << ", SensorID: "            << _currentSensorID << std::endl;
   
     streamlog_out(MESSAGE2) << " Particle type is "      << currentSimTrkHit->getMCParticle()->getPDG()
                             << ", entry: "               << stepDecoder(currentSimTrkHit)["isEntry"]  
                             << ", exit: "                << stepDecoder(currentSimTrkHit)["isExit"]  << std::endl 
                             << std::setprecision(3)
                             << " Edep is "               << currentSimTrkHit->getdEdx()*GeV/keV << " keV" << std::endl
                             << " Length is "             <<  currentSimTrkHit->getPathLength()/um << " um"          
                             << std::setprecision(0)
                             << std::endl;
      
     
     // The isEntry bit flags the first geant4 after a particle entered a sensitive medium 
     if ( stepDecoder(currentSimTrkHit)["isEntry"] == 1 ) {
       
       // Get particle entering sensor
       MCParticle * currentMCPart = dynamic_cast<MCParticle *>(currentSimTrkHit->getMCParticle());  
       
       // True ionization along the path of the primary particle 
       double trueEDep = 0; 
       
       _rootPrimaryPDG =  currentMCPart->getPDG() ; 
       
       // Find last step of the current track in current sensor; there are only two possibilities: 
       // a) current track is stopped in current sensor -> isExit never appears, a new track starts    
       // b) current track leaves current sensor -> check isExit flag for steps of current track 
       // note: all secondary particles created inside the sensor are just skipped, not taken into account!!!        
        
       SimTrackerHit * lastSimTrkHit = currentSimTrkHit; 
       
       
       // Find last steps of current particle in this detector
       // note: SimTrackerHits are expected to be ordered into particle tracks  
       for (int eHit=iHit; eHit<nHits; ++eHit) {   
          
         SimTrackerHit * nextSimTrkHit = dynamic_cast<SimTrackerHit*>(simTrkHits->getElementAt(eHit));
         MCParticle * nextMCPart = dynamic_cast<MCParticle *>(nextSimTrkHit->getMCParticle());          
        
         // check if current particle is stopped   
         if ( currentMCPart != nextMCPart ) {  
            
           // new track started, last step of current particle is the previous step   
           lastSimTrkHit = dynamic_cast<SimTrackerHit*>(simTrkHits->getElementAt(eHit-1));
           streamlog_out(MESSAGE2) << " Track is stopped in sensor!! " << std::endl;
           break; 
            
         } else {
           
           // Ok, nextSimTrkHit belongs to current particle
           trueEDep += nextSimTrkHit->getdEdx()*GeV; 
           
           // Check if particle leaves the sensor.   
           if (stepDecoder(nextSimTrkHit)["isExit"] == 1) {
             lastSimTrkHit = nextSimTrkHit;
             streamlog_out(MESSAGE2) << " Track leaves sensor!! " << std::endl; 
             break; 
           }
               
         } 
       }
        
       // Entry step into current sensor 
       SpacePoint entryStep;
        
       // Transform step to local sensor coordinates
       TransformToLocal(currentSimTrkHit, entryStep);
       double entryStepLength  = currentSimTrkHit->getPathLength()*mm;
       
       // Calculate entry point 
       Hep3Vector entryPoint = entryStep.position - entryStep.direction*entryStepLength/2.;
       
       // Last step in current sensor
       SpacePoint lastStep;
       
       // Transform geant4 step to local sensor coordinates
       TransformToLocal(lastSimTrkHit, lastStep);
       double lastStepLength  = lastSimTrkHit->getPathLength()*mm;
        
       // Calculate entry point 
       Hep3Vector lastPoint  = lastStep.position + lastStep.direction*lastStepLength/2.;
             
       // 
       // Finally, get truht hit position as mean of entry position and last position
       Hep3Vector meanPosition = 0.5*(lastPoint + entryPoint);
       
       // 
       // Finally finally, translate in UV to the centre of active area 
       meanPosition[1] -= _geometry->getSensorWidth(_currentLayerID, _currentSensorID)/2.;
       meanPosition[2] -= _geometry->getSensorLength(_currentLayerID, _currentSensorID)/2.;
       
       streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                               << std::setprecision(3)
                               << " Truth hit global[mm]: (" << meanPosition.getX()/mm
                               << ", "                       << meanPosition.getY()/mm
                               << ", "                       << meanPosition.getZ()/mm << ")"
                               << std::setprecision(0)
                               << std::endl << std::endl;
          
       #ifdef ROOT_OUTPUT      
       
       // fill hit tree variables             
       _rootX = meanPosition.getX()/mm;               
       _rootY = meanPosition.getY()/mm;  
       _rootZ = meanPosition.getZ()/mm;              
       _rootTOT = currentMCPart->getEnergy()*GeV/keV;                     
       
 
       _rootFile->cd("");
       _rootTreeHits->Fill();
       
       #endif
       

          
       //Create truth hit  
       streamlog_out(MESSAGE2) << " Producing a truth hit ... " << std::endl; 
       TrackerHitImpl * trkHit = new TrackerHitImpl();
        
       // Final output uses TB axis convention: 
       //  X -> (TB) Z , Z -> (TB) Y and Y == (TB) X     
       
       double pos[3] = {0.,0.,0.};
       pos[0] = meanPosition.getY()/mm;   // TB X
       pos[1] = meanPosition.getZ()/mm;   // TB Y
       pos[2] = _currentSensorID;         // plane number 
       trkHit->setPosition( &pos[0] );
       
       // Set covariance matrix as (c_xx, c_yy, c_xy=0)
       float cov[3] = {0.,0.,0.}; 
       cov[0] = 1; 
       cov[1] = 1;   
       trkHit->setCovMatrix ( &cov[0] );
       
       // Set quality  
       trkHit->setType(0);    
       
       // debug 
       streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                               << std::setprecision(3)
                               << " Local truth sensor " << _currentSensorID << " hit[mm]: (" << pos[0]/mm
                               << ", "                                        << pos[1]/mm
                               << ", "                                        << pos[2]/mm << ")"
                               << std::setprecision(0)
                               << std::endl << std::endl;
        
        
       // Eloss
       trkHit->setdEdx( trueEDep/keV );
        
       // Save only the entry simTrackerHit
       trkHit->rawHits().push_back(dynamic_cast<SimTrackerHit*>(currentSimTrkHit));
       
       // Add tracker hit to truth hit collection
       truthhit->push_back(trkHit);
       
     } // end of isEntry 
     
   } // simhit loop
    
}

//
// Method updating PXD map (containing all digits) using given vector of digits
//
void SiPxlDigi::UpdateDigitsMap(DigitsMap & digitsMap, DigitVec & digits)
{
   // Go through all digits in given vector
   DigitVec::iterator  iterDigitVec;

   for (iterDigitVec=digits.begin(); iterDigitVec!=digits.end(); iterDigitVec++) {

      Digit * digit = *iterDigitVec;

      // Describe sensor by unique ID & pixel by unique ID
      int uniqSensorID = _geometry->encodeSensorID(_currentLayerID, _currentLadderID, _currentSensorID);
      int uniqPixelID  = _geometry->encodePixelID(_currentLayerID, _currentSensorID, digit->cellIDRPhi, digit->cellIDZ);

      // Find if sensor already has some signal
      if (digitsMap.find(uniqSensorID)!=digitsMap.end()) {

         // Find if pixel already has some signal --> update
         if(digitsMap[uniqSensorID].find(uniqPixelID)!=digitsMap[uniqSensorID].end()) {

            digitsMap[uniqSensorID][uniqPixelID]->charge += digit->charge;
            
            // Release memory
            delete digit;
         }
         // Record info
         else {

            digitsMap[uniqSensorID][uniqPixelID] = digit;
         }
      }
      // Record info
      else {

         digitsMap[uniqSensorID][uniqPixelID] = digit;
      }
   }
}

//
// Method for calculating electronic effects on signal digits
//
void SiPxlDigi::ProduceNoiseEffects(DigitsMap & digitsMap)
{
  
  // prepare navigation in digits map 
  DigitsMap::const_iterator iterDigitsMap;
  DigitMap::const_iterator iterSensorMap;

  // Go through all sensors
  for (iterDigitsMap=digitsMap.begin(); iterDigitsMap!=digitsMap.end(); iterDigitsMap++) {
      
    // Get layerID, ladderID & sensorID corresponding to current sensor
    _geometry->decodeSensorID(_currentLayerID, _currentLadderID, _currentSensorID, iterDigitsMap->first);
     
    // Get sensor noise level
    double elNoise = _elNoise[_currentSensorID];  

    // Denote PXD sensor map as
    DigitMap digitsSensorMap = iterDigitsMap->second;

    // Go through all signal digits and add Poisson smearing and electronics effects
    for (iterSensorMap=digitsSensorMap.begin(); iterSensorMap!=digitsSensorMap.end(); iterSensorMap++) {
      
      Digit * digit = iterSensorMap->second;

      // Print
      streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                               << std::setprecision(3)
                               << "   Pos [mm]: ( "     << _sensorThick/2./mm  << ", " << digit->cellPosRPhi/mm << ", " << digit->cellPosZ/mm
                               << " , q: " << digit->charge    << " )"
                               << std::setprecision(0)
                               << std::endl;	
      
      // Add Poisson smearing
      if (_PoissonSmearing) {
	 
        // Current digit charge
        double charge = digit->charge;

	// For big charge assume Gaussian distr.
	if (charge > (1000.*e)) {
	   
	  double sigma  = sqrt(charge);
	  digit->charge = myRng->Gaus(charge,sigma);     
	}
	// Otherwise Poisson distr.
	else {
	  digit->charge = myRng->Poisson(charge);    
	}

      } // Poisson smearing

      // Add electronics effects
      if (_electronicEffects) {


	// assume gaussian readout noise
	double noise   = myRng->Gaus(0,elNoise);     
	digit->charge += noise;

        // now, simulate ADC
        if (_ADC) digit->charge = (double) getInADCUnits(digit->charge);

        // now, zero suppression
        double thresholdCut = _ZSthr;
        if (thresholdCut==0) thresholdCut++;
        
        if ( digit->charge < thresholdCut ) digit->charge = 0; 
	
      } // Electronics effects

      // Print
      streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                               << std::setprecision(3)
                               << "   Pos [mm]: ( "     << _sensorThick/2./mm  << ", " << digit->cellPosRPhi/mm << ", " << digit->cellPosZ/mm
                               << " , q: " << digit->charge    << " )"
                               << std::setprecision(0)
                               << std::endl;

       
    } // end loop: sensor map 

  } // end loop: PXD map

}


//
// Method producing noise digits and updating PXD map containing all digits
//
void SiPxlDigi::ProduceNoiseDigits(DigitsMap & digitsMap)
{
   // Initialize
   Digit * digit = 0; 
   
   // Go through all sensors & generate noise digits
   for (short int iLayer=0; iLayer<_geometry->getNPXDLayers(); iLayer++) {
      for (short int iLadder=0; iLadder<_geometry->getNLadders(iLayer); iLadder++) {
         for (short int iSensor=0; iSensor<_geometry->getNSensors(iLayer); iSensor++) {
            
            // Describe sensor by unique ID
            int uniqSensorID = _geometry->encodeSensorID(iLayer, iLadder, iSensor);
            
           
            
            // Average number of noise pixels
            double meanNoisePixels = _Fraction*_geometry->getSensorNPixelsInRPhi(iLayer,iSensor)*_geometry->getSensorNPixelsInZ(iLayer,iSensor);
            
            // Total number of noise pixels has Poison distribution
            int fractionPixels = myRng->Poisson(meanNoisePixels);  

            // Generate noise digits
            for (int iNoisePixel=0; iNoisePixel<fractionPixels; iNoisePixel++) {
                    
              int iPixelRPhi  = int(myRng->Uniform(_geometry->getSensorNPixelsInRPhi(iLayer, iSensor)));
              int iPixelZ     = int(myRng->Uniform(_geometry->getSensorNPixelsInZ(iLayer, iSensor)));
              
              if (iPixelRPhi==_geometry->getSensorNPixelsInRPhi(iLayer, iSensor)) iPixelRPhi--;
              if (iPixelZ   ==_geometry->getSensorNPixelsInZ(   iLayer, iSensor)) iPixelZ--;
                   
              // Describe pixel by unique ID
              int uniqPixelID  = _geometry->encodePixelID(iLayer, iSensor, iPixelRPhi, iPixelZ);
                
              // Find if pixel doesn't already have some signal+noise or just noise
              if( !(digitsMap[uniqSensorID].find(uniqPixelID)!=digitsMap[uniqSensorID].end()) ) {
                     
                // Create a noise charge  
                double charge = _ZSthr; 
                     
                // Create new noise digit      
                digit = new Digit;        		       	
                digit->cellIDRPhi  = iPixelRPhi;
                digit->cellIDZ     = iPixelZ;
                         
                digit->cellPosRPhi = _geometry->getPixelPosInRPhi(iLayer, iSensor, _bricked, iPixelRPhi, iPixelZ);
                digit->cellPosZ    = _geometry->getPixelPosInZ(iLayer, iSensor, iPixelZ);
                       
                digit->charge      = charge;
                            
                // Record it
                digitsMap[uniqSensorID][uniqPixelID] = digit;
                  
              }
            } // For noise pixels
            
         }
      }
   } // 3xFor: layers, ladders, sensors

}

//
// Method returning collected charge in ADC units
//

int SiPxlDigi::getInADCUnits(double charge) {

   static double ADCunit = _ADCRange/short(pow(2, _ADCBits));

   return int(charge/ADCunit);
}


//
// Method providing mathematical round
//

int SiPxlDigi::round(double num) {

   return (int)(num+0.5);
}

//
// Method printing processor parameters
//
void SiPxlDigi::printProcessorParams() const
{
   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << DUNDERL
                            << DBLUE
                            << "SiPxlDigi parameters:"
                            << ENDCOLOR
                            << " "
                            << std::endl  << std::endl;

   streamlog_out(MESSAGE3)  << std::setiosflags(std::ios::fixed | std::ios::internal )
                            << std::setprecision(2)                        << std::endl;
   if (_bricked) {
   streamlog_out(MESSAGE3)  << "  Bricked structure used?:           " << "   yes"     << std::endl;
   }
   else {
   streamlog_out(MESSAGE3)  << "  Bricked structure used?:           " << "    no"     << std::endl;
   }
   if (_electronicEffects) {
   //streamlog_out(MESSAGE3)  << "  El. noise [fC]:                    " << std::setw(6) << _elNoise/fC              << std::endl;
   }
   if (_PoissonSmearing) {
   streamlog_out(MESSAGE3)  << "  Electronics - Poisson smearing?:   " << "   yes"     << std::endl;
   }
   else {
   streamlog_out(MESSAGE3)  << "  Electronics - Poisson smearing?:   " << "    no"     << std::endl;
   }
   if (_integrationWindow) {
   streamlog_out(MESSAGE3)  << "  Integration window [us]:           " << std::setw(6) << _startIntegration/us     << " ; " << _stopIntegration/us << std::endl;
   }
   streamlog_out(MESSAGE3)  << "  Space precision [um]:              " << std::setw(6) << _maxSegmentLength/um        << std::endl
                            << "  ZS threshold for zs:              " << std::setw(6) << _ZSthr                   << std::endl
                            << "  Tangent of Lorentz angle:          " << std::setw(6) << _tanLorentzAngle         << std::endl
                            << std::resetiosflags(std::ios::showpos)
                            << std::setprecision(0)
                            << std::endl;
}

} // Namespace






