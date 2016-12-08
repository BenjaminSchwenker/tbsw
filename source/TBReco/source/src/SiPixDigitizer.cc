// SiPixDigitizer
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "SiPixDigitizer.h"


// Include DEPFETTrackTools 
#include "DEPFET.h" 
#include "MatrixDecoder.h"
#include "PhysicalConstants.h"

// Include basic C
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <string>
#include <map>

// Include CLHEP classes
#include <CLHEP/Vector/Rotation.h>

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
using namespace CLHEP;
using namespace lcio ;
using namespace marlin ;

namespace depfet {

  //
  // Instantiate this object
  //
  SiPixDigitizer aSiPixDigitizer ;
  
  //
  // Constructor
  //
  SiPixDigitizer::SiPixDigitizer() : Processor("SiPixDigitizer")
  {
    
    // Processor description
    _description = "SiPixDigitizer: A rather generic silicon pixel detector digitizer" ;
    
    //
    // Input collections  
    registerInputCollection (LCIO::SIMTRACKERHIT, "SimTrackerHitCollectionName",
                             "Collection name for SimTrackerHits",
                             m_SimTrackerHitCollectionName, string ("SimTrackerHits") );
      
    registerOutputCollection(LCIO::TRACKERDATA, "DigitCollectionName",
                             "Collection name for Digits",
                              m_digitCollectionName, std::string("zsdata"));
    
    std::vector<int> initFilterIDs;
    registerProcessorParameter ("FilterIDs",
                                "Apply digitization only to SimTrackerHits for sensors having DAQ IDs in this list",
                                _filterIDs, initFilterIDs);
      
    registerProcessorParameter( "SiBulkDoping",
                                "Effective bulk doping concentration, in um^-3",
                                m_bulkDoping,
                                double(10));
    
    registerProcessorParameter( "BackPlaneVoltage",
                               "Back contact voltage wrt. source, in Volts",
                               m_backVoltage,
                               double(-30));
   
    registerProcessorParameter( "TopPlaneVoltage",
                                "Top contact voltage wrt. source, in Volts",
                                m_topVoltage,
                                double(-5));
    
    registerProcessorParameter( "uSideBorderLength",
                                "Border length along u side, in um",
                                m_uSideBorderLength,
                                double(3));
    
    registerProcessorParameter( "vSideBorderLength",
                                "Border length along v side, in um",
                                m_vSideBorderLength,
                                double(3));  
       
    registerProcessorParameter( "MaxSegmentLength",
                               "Maximum distance between ionization points (in mm)",
                               m_maxSegmentLength,
                               double(0.005));
    
    registerProcessorParameter( "ElectronGroupSize",
                                "Split Signalpoints in smaller groups of N electrons (in e)",
                                m_eGroupSize,
                                double(100));
    
    registerProcessorParameter( "ElectronStepTime",
                                "Time step for tracking electron groups in readout plane (in ns)",
                                m_eStepTime,
                                double(0.3));
    
    registerProcessorParameter( "TanLorentz",
                                "Tangent of Lorentz angle",
                                m_tanLorentzAngle,
                                double(0.0));
    
    registerProcessorParameter( "IntegrationWindow",
                                "Use integration window?",
                                m_integrationWindow,
                                bool(true));
     
    registerProcessorParameter( "StartIntegration",
                                "Only Simulated hits after the StartIntegration time in ns will be digitized",
                                m_startIntegration,
                                double(0.0));
   
    registerProcessorParameter( "StopIntegration",
                                "Only Simulated hits before the StopIntegration time in ns will be digitized",
                                m_stopIntegration,
                                double(20000.0));
    
    registerProcessorParameter( "ElectronicEffects",
                                "Apply electronic effects?",
                                m_electronicEffects,
                                bool(true));
    
    registerProcessorParameter( "ElectronicNoise",
                                "Noise added by the electronics, set in ENC (one value for each detector)",
                                m_elNoise,
                                double(100) );
     
    registerProcessorParameter( "NoiseFraction",
                                "Set fraction of noise hits per sensor per readout frame",
                                m_noiseFraction,
                                float(0.001) );
    
    registerProcessorParameter( "ZSThreshold",
                                "Keep only hits >= threshold in ADC units",
                                m_zsThreshold,
                                float(4.0) );
   
    registerProcessorParameter( "ADC",
                                "Simulate ADC?",
                                m_useADC,
                                bool(false));
    
    registerProcessorParameter( "ADCRange",
                                "Set analog-to-digital converter range 0 - ? (in e)",
                                m_ADCRange,
                                int(50000));
    
    registerProcessorParameter( "ADCBits",
                                "Set how many bits the ADC uses",
                                _ADCBits,
                                int(8));
    
   
  }

  //
  // Method called at the beginning of data processing
  //
  void SiPixDigitizer::init() {
    
    // Initialize variables
    m_nRun = 0 ;
    m_nEvt = 0 ;
    
    // Set variables in appropriate physical units
    m_ADCRange             *= e;
    m_bulkDoping           *= 1/um/um/um; 
    m_topVoltage           *= V; 
    m_backVoltage          *= V; 
    m_maxSegmentLength     *= mm;
    m_eGroupSize           *= e; 
    m_eStepTime            *= ns; 
    m_startIntegration     *= ns;
    m_stopIntegration      *= ns;
    m_uSideBorderLength    *=um;
    m_uSideBorderLength    *=um;
    m_elNoise              *= e; 
    
    // Read detector constants from gear file
    m_detector.ReadGearConfiguration();      
    
    // Print set parameters
    printProcessorParams();
    
    // CPU time start
    m_timeCPU = clock()*ms/1000;
  }

  //
  // Method called for each run
  //
  void SiPixDigitizer::processRunHeader(LCRunHeader * run)
  {
  
    // Print run number
    streamlog_out(MESSAGE3) << "Processing run: "
                            << (run->getRunNumber())
                            << std::endl << std::endl;
    
    m_nRun++ ;
  }
  
  //
  // Method called for each event
  //
  void SiPixDigitizer::processEvent(LCEvent * evt)
  {
   
    //
    // Open collections
    try {
      
      // Open SimTrackerHit collection
      LCCollection * simHitCol = evt->getCollection( m_SimTrackerHitCollectionName );
      
      // Set collection decoder
      CellIDDecoder<SimTrackerHit> cellIDDec(simHitCol);
      
      // Number of SimTrackerHits
      int nSimHits = simHitCol->getNumberOfElements();
      
      // Initialize map of all digits -> needed for proper clustering
      DigitsMap digitsMap;
      
      //
      // Loop over SimTracker hits;
      streamlog_out(MESSAGE2) << " Producing all digits ..." << std::endl;
      
      for (int i=0; i<nSimHits; ++i) {
        
        SimTrackerHit * simTrkHit = dynamic_cast<SimTrackerHit*>(simHitCol->getElementAt(i));
        
        // Set current - layer ID, ladder ID and sensor ID
        m_ipl = cellIDDec(simTrkHit)["layer"];
         
        // Cut on simHit creation time --> simulate integration time of a sensor (if option switched on))
	    if ((simTrkHit != 0) && (m_integrationWindow)) {
	      if (simTrkHit->getTime()*ns < _startIntegration || simTrkHit->getTime()*ns > _stopIntegration) {
	        continue;
	      }		
	    }
        
        // Print
        streamlog_out(MESSAGE2) << " Sensor plane: "  << m_ipl << std::endl;
        streamlog_out(MESSAGE2) << std::setiosflags(std::ios::fixed | std::ios::internal )
                                << std::setprecision(3)
                                << " Simulated hit global[mm]: (" << simTrkHit->getPosition()[0]
                                << ", "                           << simTrkHit->getPosition()[1]
                                << ", "                           << simTrkHit->getPosition()[2] << ")"
                                << std::setprecision(0)
                                << std::endl;
         
        // Only filter SimTrackerHits from sensor in filter list
        /*
        if (_geometry->getLayerType(_currentLayerID) != pixel ) {
          
          continue;
        }
        */
        
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
      
      if (_electronicEffects) {
        
        // Produce noise digits and update PXD map
        streamlog_out(MESSAGE2) << " Producing noisy digits ..." << std::endl;
	    
	    // add electronic noise to digits
	    ProduceNoiseEffects(digitsMap);
        
	    // create pure noise digits
	    ProduceNoiseDigits(digitsMap);
      }
      
      // 
      // Produce sparse pixel output collection  
      streamlog_out(MESSAGE1) << " Writing digits ..." << std::endl;
	  
      // Create sparse data collection
      LCCollectionVec * sparseDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);       
        
      ProduceSparsePixels( sparseDataCollection  , digitsMap); 
	  
      // Store stuff in LCIO file 
      evt->addCollection(sparseDataCollection, _sparseDataCollectionName);
        
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
    
    _nEvt ++ ;
  }

  //
  // Method called after each event to check the data processed
  //
  void SiPixDigitizer::check( LCEvent * evt ) { }
  
  //
  // Method called after all data processing
  //
  void SiPixDigitizer::end()
  {
    // CPU time end
    m_timeCPU = clock()*ms/1000 - m_timeCPU;
    
    // Print message
    streamlog_out(MESSAGE3) << std::endl
                            << " "
                            << "Time per event: "
                            << std::setiosflags(std::ios::fixed | std::ios::internal )
                            << std::setprecision(3)
                            << m_timeCPU/m_nEvt/ms
                            << " ms"
                            << std::endl
                            << std::setprecision(3)
                            << std::endl
                            << " "
                            << "Processor succesfully finished!"
                            << std::endl;
  
  }

  //
  // Method transforming hit to local coordinates
  //
  void SiPixDigitizer::TransformToLocal(const SimTrackerHit * simTrkHit, SpacePoint & hitLocal)
  {
    Hep3Vector position(simTrkHit->getPosition()[0]*mm ,simTrkHit->getPosition()[1]*mm ,simTrkHit->getPosition()[2]*mm );
    Hep3Vector momentum(simTrkHit->getMomentum()[0]*GeV,simTrkHit->getMomentum()[1]*GeV,simTrkHit->getMomentum()[2]*GeV);
     
    // Save final results
    hitLocal.position  = position;
      
    if ( momentum.mag() != 0 ) {
       hitLocal.direction = momentum/momentum.mag();
    } else {
       hitLocal.direction = momentum; 
    }
  }
  
  //
  // Method producing ionisation points along the track
  //
  void SiPixDigitizer::ProduceIonisationPoints(const SimTrackerHit * simTrkHit, IonisationPointVec & ionisationPoints)
  {
    // Space point
    SpacePoint hitLocal;
    
    // Transform geant4 step to local sensor coordinates
    TransformToLocal(simTrkHit, hitLocal);
    
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
      /* FIXME
      if (_geometry->isPointOutOfSensor(_currentLayerID, _currentSensorID, iPoint->position)) {
         
        streamlog_out(ERROR) << std::setiosflags(std::ios::fixed | std::ios::internal )
                             << std::setprecision(3)
                             << "SiPixDigitizer::ProduceIonisationPoints - ionPoint: " << iPoint->position/mm << " out of sensor!!!"
                             << std::setprecision(0) << std::endl
                             << "SensorID is " << _currentSensorID
                             << std::endl;
       
      }
      */ 
      // Print
      streamlog_out(MESSAGE1) << "  Hit local ionPoints (ionisation): " << std::endl;
      streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(3)
                              << "   Pos [mm]: ( " << iPoint->position.getX()/mm << ", " << iPoint->position.getY()/mm << ", " << iPoint->position.getZ()/mm << " )"
                              << " , dE [keV]: "   << iPoint->eLoss/keV
                              << std::setprecision(0)
                              << std::endl;
     
    } // end segment loop  
  }
  
  //
  // Method producing signal points on each sensor
  //
  void SiPixDigitizer::ProduceSignalPoints(const IonisationPointVec & ionisationPoints, SignalPointVec & signalPoints)
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
  void SiPixDigitizer::ProduceDigits(const SignalPointVec & signalPoints, DigitVec & digits) 
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
  // Method producing sparsified pixels output from all digits (noise + signal)
  //
  void SiPixDigitizer::ProduceSparsePixels(LCCollectionVec * sparseDataCollection , const DigitsMap & digitsMap) 
  {
  
    // Print
    streamlog_out(MESSAGE2) << " Writing digits " << std::endl;
    
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
        
    } // end for: loop over sensor
  }

  //
  // Method updating PXD map (containing all digits) using given vector of digits
  //
  void SiPixDigitizer::UpdateDigitsMap(DigitsMap & digitsMap, DigitVec & digits)
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
        } else {
          digitsMap[uniqSensorID][uniqPixelID] = digit;
        }
      } else {
        digitsMap[uniqSensorID][uniqPixelID] = digit;
      }
    }
  }
  
  //
  // Method for calculating electronic effects on signal digits
  //
  void SiPixDigitizer::ProduceNoiseEffects(DigitsMap & digitsMap)
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
        double charge = digit->charge;
           
	    // For big charge assume Gaussian distr.
	    if (charge > (1000.*e)) {
	      double sigma  = sqrt(charge);
	      digit->charge = myRng->Gaus(charge,sigma);     
	    } else {
	      digit->charge = myRng->Poisson(charge);    
	    }
        
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
  void SiPixDigitizer::ProduceNoiseDigits(DigitsMap & digitsMap)
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
    }// 3xFor: layers, ladders, sensors
  }
  
  //
  // Method returning collected charge in ADC units
  //
  int SiPixDigitizer::getInADCUnits(double charge) {
    static double ADCunit = _ADCRange/short(pow(2, _ADCBits));
    return int(charge/ADCunit);
  }

  //
  // Method providing mathematical round
  //
  int SiPixDigitizer::round(double num) { return (int)(num+0.5);}
   
  //
  // Method printing processor parameters
  //
  void SiPixDigitizer::printProcessorParams() const
  {
     
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "SiPixDigitizer development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
    
  }

} // Namespace






