// SiPixDigitizer
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "SiPixDigitizer.h"

// Include TBTools  
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
                             m_SimTrackerHitCollectionName, std::string ("SimTrackerHits") );
      
    registerOutputCollection(LCIO::TRACKERDATA, "DigitCollectionName",
                             "Collection name for Digits",
                              m_digitCollectionName, std::string("zsdata"));
    
    std::vector<int> initFilterIDs;
    registerProcessorParameter ("FilterIDs",
                                "Apply digitization only to SimTrackerHits for sensors having DAQ IDs in this list",
                                m_filterIDs, initFilterIDs);
      
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
   
    registerProcessorParameter( "FrontEndType",
                                "Choose front-end electronics: ADC(0) or Comparator(1) or None(2) ",
                                m_frontEndType,
                                int(1));
    
    registerProcessorParameter( "ADCRange",
                                "Set analog-to-digital converter range 0 - ? (in e)",
                                m_ADCRange,
                                int(50000));
    
    registerProcessorParameter( "ADCBits",
                                "Set how many bits the ADC uses",
                                m_ADCBits,
                                int(8));
    
    registerProcessorParameter( "ComparatorThrehold",
                                "Set threshold for comparator (in e)",
                                m_ComparatorThr,
                                float(1000));
    
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
    m_ComparatorThr        *= e;
    m_bulkDoping           *= 1/um/um/um; 
    m_topVoltage           *= V; 
    m_backVoltage          *= V; 
    m_maxSegmentLength     *= mm;
    m_eGroupSize           *= e; 
    m_eStepTime            *= ns; 
    m_startIntegration     *= ns;
    m_stopIntegration      *= ns;
    m_uSideBorderLength    *=um;
    m_vSideBorderLength    *=um;
    m_elNoise              *= e; 
    
    // Threshold of zero is never valid 
    if (m_zsThreshold==0) m_zsThreshold++;
    // Comparator mode needs threshold one
    if (m_frontEndType == 1)  m_zsThreshold = 1;
    
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
    
    // Print event number
    if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                  << (evt->getEventNumber())
                                                                  << std::endl << std::endl;
    
    //
    // Open collections
    try {
      
      // Maybe this event got recorded because of a fake trigger 
      bool isFakeTrigger = false;
      try {
        evt->getCollection( "FakeTrigger" );
        isFakeTrigger = true;
      } catch (lcio::DataNotAvailableException& e) {}  
      
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
      streamlog_out(MESSAGE2) << " Producing all SimTrackerHits ..." << std::endl;
      
      for (int i=0; i<nSimHits; ++i) {
        
        SimTrackerHit * simTrkHit = dynamic_cast<SimTrackerHit*>(simHitCol->getElementAt(i));
        
        // Set current - layer ID, ladder ID and sensor ID
        m_sensorID = cellIDDec(simTrkHit)["sensorID"];
        m_ipl = m_detector.GetPlaneNumber(m_sensorID);
         
        streamlog_out(MESSAGE1) << " Found SimTrackerHit with sensorID: " << m_sensorID  
                                << std::setprecision(10)
                                << " at time[s]:" << simTrkHit->getTime()
                                << std::setprecision(0) << std::endl;

        streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                                << std::setprecision(10)
                                << "Integration start/stop[s]: " << m_startIntegration << "/" << m_stopIntegration
                                << std::setprecision(0) << std::endl;
                              
 
        // Cut on simHit creation time --> simulate integration time of a sensor (if option switched on))
	    if ((simTrkHit != 0) && (m_integrationWindow)) {
	      if (simTrkHit->getTime() < m_startIntegration || simTrkHit->getTime() > m_stopIntegration) {
            streamlog_out(MESSAGE1) << " Skipped simHit out of integration time" << std::endl; 
	        continue;
	      }		
	    }
        
        // Cut on simHit creation time --> events with fake trigger never have simhits at t=0 
	    if ((simTrkHit != 0) && (isFakeTrigger)) {
	      if (simTrkHit->getTime() == 0 ) {
            streamlog_out(MESSAGE1) << " Skipped simHit at t=0 because of a fake trigger" << std::endl; 
	        continue;
	      }		
	    }
        
        if ( std::find(m_filterIDs.begin(), m_filterIDs.end(), m_sensorID) == m_filterIDs.end() ) {
          streamlog_out(MESSAGE2) << " Ignore SimTrackerHit with sensorID: "  << m_sensorID << std::endl;
          continue;
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
        ProduceSignalDigits(signalPoints, digits);
        
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
      
      if (m_electronicEffects) {
        
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
      LCCollectionVec * digitCol = new LCCollectionVec(LCIO::TRACKERDATA);       
        
      WriteDigitsToLCIO( digitCol  , digitsMap); 
	  
      // Store stuff in LCIO file 
      evt->addCollection(digitCol, m_digitCollectionName);
        
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
    
    m_nEvt ++ ;
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
    Hep3Vector momentum(simTrkHit->getMomentum()[0],simTrkHit->getMomentum()[1],simTrkHit->getMomentum()[2]);
    
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
    
    Hep3Vector entryPoint = hitLocal.position - hitLocal.direction*trackLength/2.;
    Hep3Vector exitPoint  = hitLocal.position + hitLocal.direction*trackLength/2.;
      
    // Check entry and exit point are within sensor boundaries  
    if (  m_detector.GetDet(m_ipl).isPointOutOfSensor( entryPoint.getX(), entryPoint.getY() , entryPoint.getZ() ) ) {
      streamlog_out(MESSAGE2) << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(3)
                              << "SiPixDigitizer::ProduceIonisationPoints - ionPoint: " << entryPoint/mm << " out of sensor!!!"
                              << std::setprecision(0) << std::endl
                              << "SensorID is " << m_sensorID
                              << std::endl;
    
            
      return; 
    }
    if (  m_detector.GetDet(m_ipl).isPointOutOfSensor( exitPoint.getX(), exitPoint.getY() , exitPoint.getZ() ) ) {
      streamlog_out(MESSAGE2) << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(3)
                              << "SiPixDigitizer::ProduceIonisationPoints - ionPoint: " << exitPoint/mm << " out of sensor!!!"
                              << std::setprecision(0) << std::endl
                              << "SensorID is " << m_sensorID
                              << std::endl;
    
            
      return; 
    }
      
    // For charged particles: energy loss is smeared along the step in 
    // small segments. 
     
    int numberOfSegments = int(trackLength/m_maxSegmentLength) + 1;
     
    // Energy deposit along G4 step 
    double Edep = (simTrkHit->getdEdx()*GeV);
    // Calculate mean energy loss in each segment
    double dEMean = Edep/((double)numberOfSegments);
     
    ionisationPoints.resize(numberOfSegments);
          
    // Set ionisation points
    for (int i=0; i<numberOfSegments; i++) {
       
      IonisationPoint * iPoint = new IonisationPoint;
       
      iPoint->position = entryPoint + hitLocal.direction*trackLength/numberOfSegments*(i+0.5);
      iPoint->eLoss    = dEMean; 
      ionisationPoints[i] = iPoint;
      
      // Print
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
    
    
    
    for (int i=0; i<numberOfIonPoints; i++) {
      
      // Get current ionisation point
      IonisationPoint * iPoint = ionisationPoints[i];
      
      //  Charge cloud created at distance to top plane 
      double w =  m_detector.GetDet(m_ipl).GetSensitiveThickness()/2. - iPoint->position.getZ();
      w = sqrt( w * w );
      
      // Potential valley is at distance to top plane
      double w0 = m_detector.GetDet(m_ipl).GetSensitiveThickness()/2. - 0.02;
      double dw = 0.003;  
      
      // Drift time into potential valley -  
      double td = 1.4*1.4*um*um/( 2* Utherm * e_mobility);
        
      if ( w > w0 + dw ) {
        td += Perm_Si * ::log( (w - w0)/dw ) / ( e_mobility * e * m_bulkDoping);
      } else if ( w < w0 - dw )  {
        td += Perm_Si * ::log( (w0 - w)/dw ) / ( e_mobility * e * m_bulkDoping);
      }
      
      // Diffusive spread in lateral plane
      double sigmaDiffus = sqrt( 2 * Utherm * e_mobility * td ); 
      
      double sigmaU = sigmaDiffus;
      double sigmaV = sigmaDiffus; 
      
      //  After Lorentz shift
      double onPlaneU = iPoint->position.getX() + m_tanLorentzAngle * (w - w0);
      double onPlaneV = iPoint->position.getY();
        
      // Save info in signal point
      SignalPoint * sPoint = new SignalPoint;
      
      sPoint->position.setX(onPlaneU);
      sPoint->position.setY(onPlaneV);
      sPoint->position.setZ(w0);
      
      sPoint->sigma.setX(sigmaU);
      sPoint->sigma.setY(sigmaV);
      sPoint->sigma.setZ(0);
      
      // Charge in electrons
      sPoint->charge = iPoint->eLoss/Eeh * e;
      
      // Save signal point
      signalPoints[i] = sPoint;
       
      // Print
      streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(3)
                              << "   Pos [mm]:( " << sPoint->position.getX()/mm << ", " << sPoint->position.getY()/mm << ", " << sPoint->position.getZ()/mm << " )"
                              << " , q [e-]: "    << sPoint->charge
                              << std::setprecision(0)
                              << std::endl;
       
      streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(3)
                              << "   Sigma [um]:( " << sPoint->sigma.getX()/um << ", " << sPoint->sigma.getY()/um << ", " << sPoint->sigma.getZ()/um << " )"
                              << std::setprecision(0)
                              << std::endl;
      
    }
  }
  
  //
  // Method producing digits from signal points
  //
  void SiPixDigitizer::ProduceSignalDigits(const SignalPointVec & signalPoints, DigitVec & digits) 
  {
     
    // Calculate number of signal points
    int numberOfSigPoints = signalPoints.size(); 
    
    // Run over all signal points and create digits
    digits.clear();
    
    // Random walk step length 
    double sigmaDiffus = sqrt( 2*Utherm*e_mobility*m_eStepTime ) ;
    
    streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(3)
                              << "   Random walk step length [um]: " << sigmaDiffus/um 
                              << std::setprecision(0)
                              << std::endl;
     
    // Internal pixel borders   
    double halfwidthU   = m_detector.GetDet(m_ipl).GetPitchU()/2. - m_uSideBorderLength;          
    double halfwidthV   = m_detector.GetDet(m_ipl).GetPitchV()/2. - m_vSideBorderLength;         
    
    streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                            << std::setprecision(3)
                            << "   u halfwidth [mm]: " <<  halfwidthU << ", v halfwidth [mm]: " << halfwidthV
                            << std::setprecision(0)
                            << std::endl;
    
    for (int i=0; i<numberOfSigPoints; ++i) {
      
      // Get current signal point
      SignalPoint * sPoint = signalPoints[i];
      
      // Calculate centre of gaussian charge cloud
      double centreU = sPoint->position.getX();
      double centreV = sPoint->position.getY();
      
      // Calculate width of charge cloud
      double sigmaU  = sPoint->sigma.getX();
      double sigmaV  = sPoint->sigma.getY();
      
      // Get number of electrons in cloud
      double clusterCharge = sPoint->charge;
      
      // Now, the signal point is split into groups of electons spread around the 
      // signal point center. Each group is tracked into an internal gate.
      int numberOfGroups = int(clusterCharge/m_eGroupSize) + 1;
      double groupCharge = clusterCharge/((double)numberOfGroups);
      
      for (int iGroup=0; iGroup<numberOfGroups; ++iGroup)  {         
        
        // Initial group position   
        double groupPosU = centreU + gRandom->Gaus(0, sigmaU);
        double groupPosV = centreV + gRandom->Gaus(0, sigmaV);
            
        int iV = m_detector.GetDet(m_ipl).GetRowFromCoord( groupPosU, groupPosV );
        int iU = m_detector.GetDet(m_ipl).GetColumnFromCoord( groupPosU, groupPosV ); 
        
        double pixelPosV = m_detector.GetDet(m_ipl).GetPixelCenterCoordV(iV, iU); 
        double pixelPosU = m_detector.GetDet(m_ipl).GetPixelCenterCoordU(iV, iU); 
           
        // control variables 
        double collectionTime = 0;
      
        
        for (int iStep = 0; iStep < 400; ++iStep) {
           
          // Calculate border of internal gate region
          double deltaV = groupPosV - pixelPosV; 
          double deltaU = groupPosU - pixelPosU;  
          
          // Check if charge is inside IG region
          if (   deltaU > - halfwidthU &&  deltaU < halfwidthU &&
                 deltaV > - halfwidthV &&  deltaV < halfwidthV        )
          {
            streamlog_out(MESSAGE1) << "   random walk finished: (iU: " << iU << ", iV: " << iV << ", istep: " << iStep << ")"  << std::endl;
            break;     
          }
            
  	      // Update position of group, possibly leaving current pixel 
          collectionTime += m_eStepTime;
          groupPosV += gRandom->Gaus(0, sigmaDiffus);
          groupPosU += gRandom->Gaus(0, sigmaDiffus);

          
          // Update charge cloud posisiton 
          iV = m_detector.GetDet(m_ipl).GetRowFromCoord( groupPosU, groupPosV );
          iU = m_detector.GetDet(m_ipl).GetColumnFromCoord( groupPosU, groupPosV ); 
        
          pixelPosV = m_detector.GetDet(m_ipl).GetPixelCenterCoordV(iV, iU); 
          pixelPosU = m_detector.GetDet(m_ipl).GetPixelCenterCoordU(iV, iU); 
           
                 
        } // end of group tracking
         
        // Now, pixel cell [iU,iV] collects group of electrons         
        // Update digit info if exists, otherwise create new one
        DigitVec::iterator iterDVec;
        bool digitFound = false;
               
        // Update digit info
        for (iterDVec=digits.begin(); iterDVec!=digits.end(); iterDVec++) {
              
          Digit * digit = (*iterDVec);
                 
          if ( (digit->cellIDU == iU) && (digit->cellIDV == iV) ) {
                   
            digit->charge += groupCharge;
            digitFound = true;
            break;
          }
        }
                 
        // Such a digit not found - create new digit
        if (digitFound == false) {
               
          Digit * digit = new Digit;
                
          digit->cellIDV = iV;
          digit->cellIDU = iU;
          digit->charge   = groupCharge;                
          digits.push_back(digit);
        } 
      } // For groups  
    } // For signalPoints
  }


  //
  // Method producing sparsified pixels output from all digits (noise + signal)
  //
  void SiPixDigitizer::WriteDigitsToLCIO(LCCollectionVec * digitCollection , const DigitsMap & digitsMap) 
  {
  
    // Print
    streamlog_out(MESSAGE2) << " Writing digits " << std::endl;
    
    // Prepare an encoder for zero suppressed pixels
    CellIDEncoder<TrackerDataImpl> sparseDataEncoder( "sensorID:6,sparsePixelType:5" , digitCollection ); 
     
    // Prepare navigation in digits map 
    DigitsMap::const_iterator iterDigitsMap;
    DigitMap::const_iterator iterSensorMap;
    
    // Go through all sensors
    for (iterDigitsMap=digitsMap.begin(); iterDigitsMap!=digitsMap.end(); iterDigitsMap++) {
        
      // Get sensorID
      int sensorID = iterDigitsMap->first;
      
      // Get handle to DigitsMap for sensor
      DigitMap digitsSensorMap = iterDigitsMap->second;
       
      // Prepare a TrackerData to store digit
      TrackerDataImpl* digitVec = new TrackerDataImpl;
      
      // Set cellID
      sparseDataEncoder["sensorID"] = sensorID; 
      sparseDataEncoder["sparsePixelType"] = 0;
      sparseDataEncoder.setCellID( digitVec  );
      
      // Loop over all digits on this sensor 
      for (iterSensorMap=digitsSensorMap.begin(); iterSensorMap!=digitsSensorMap.end(); iterSensorMap++) {
       
        Digit * currentDigit = iterSensorMap->second;
       
        double signal  = currentDigit->charge;
        
        // Note: zero suppressed pixels have signal == 0, just skip 
        if (signal == 0) continue;  
      
        int iV = currentDigit->cellIDV;
        int iU = currentDigit->cellIDU;
      
        // Store digits  
        digitVec->chargeValues().push_back( iU );
        digitVec->chargeValues().push_back( iV );
        digitVec->chargeValues().push_back( signal );   

        streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                                << std::setprecision(3)
                                << " (sensorID: " <<  sensorID  << ", iU: " << iU << ", iV: " << iV
                                << " , q[ADU]: " << signal    << " )"
                                << std::setprecision(0)
                                << std::endl;	 
       
        
      
      }
      
      // Add zs pixels to collection 
      digitCollection->push_back( digitVec  );
        
    } // end for: loop over sensor
  }

  //
  // Method updating PXD map (containing all digits) using given vector of digits
  //
  void SiPixDigitizer::UpdateDigitsMap(DigitsMap & digitsMap, DigitVec & digits)
  {
    // Go through all digits in given vector 
    for (auto iterDigitVec=digits.begin(); iterDigitVec!=digits.end(); iterDigitVec++) {
          
      Digit * digit = *iterDigitVec;
      
      int uniqSensorID = m_sensorID;
      int uniqPixelID  = m_detector.GetDet(m_ipl).encodePixelID(digit->cellIDV, digit->cellIDU);    
      
      streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                                << std::setprecision(3)
                                << " (sensorID: " <<  uniqSensorID  << ", iU: " << digit->cellIDU << ", iV: " << digit->cellIDV
                                << " , q[e-]: " << digit->charge    << " )"
                                << std::setprecision(0)
                                << std::endl;	
      
      // Find if sensor already has some signal
      if (digitsMap.find(uniqSensorID)!=digitsMap.end()) {
         
        // Find if pixel already has some signal --> update
        if(digitsMap[uniqSensorID].find(uniqPixelID)!=digitsMap[uniqSensorID].end()) {
          // Digit exists already
          digitsMap[uniqSensorID][uniqPixelID]->charge += digit->charge;
          // Release digit memery
          delete digit;
        } else {
          // store digit
          digitsMap[uniqSensorID][uniqPixelID] = digit;
        }
      } else {
        // store digits
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
      
      // Get sensor noise level
      double elNoise = m_elNoise;  
      
      // Denote PXD sensor map as
      DigitMap digitsSensorMap = iterDigitsMap->second;
      
      // Go through all signal digits and add Poisson smearing and electronics effects
      for (iterSensorMap=digitsSensorMap.begin(); iterSensorMap!=digitsSensorMap.end(); iterSensorMap++) {
      
        Digit * digit = iterSensorMap->second;
        
        // Add Poisson smearing 
        double charge = digit->charge;
           
	    // For big charge assume Gaussian distr.
	    if (charge > (1000.*e)) {
	      double sigma  = sqrt(charge);
	      digit->charge = gRandom->Gaus(charge,sigma);     
	    } else {
	      digit->charge = gRandom->Poisson(charge);    
	    }
          
	    // Add Gaussian readout noise
	    double noise   = gRandom->Gaus(0,elNoise);     
	    digit->charge += noise;
          
        // Very basic simulation of front end electronics
        if ( m_frontEndType == 0) {
          // Simulate ADC 
          digit->charge = (double) getInADCUnits(digit->charge);
        } else if ( m_frontEndType == 1) { 
          // Simulate Comparator       
          if ( digit->charge < m_ComparatorThr ) digit->charge = 0;
          else digit->charge = 1;
        }
        
        // Print
        streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                                << std::setprecision(3)
                                << " (iU: " << digit->cellIDU << ", iV: " << digit->cellIDV
                                << " , q[ADU]: " << digit->charge    << " )"
                                << std::setprecision(0)
                                << std::endl;	   
         
        // Perform zero suppression
        if ( digit->charge < m_zsThreshold ) digit->charge = 0; 
	     
      } // end loop: sensor map   
    } // end loop: digits map
  }
  
  //
  // Method producing noise digits and updating PXD map containing all digits
  //
  void SiPixDigitizer::ProduceNoiseDigits(DigitsMap & digitsMap)
  {
    // Initialize
    Digit * digit = 0; 
    
    // Go through all sensors & generate noise digits
    for (short int ipl=0; ipl < m_detector.GetNSensors(); ipl++) {
      
      m_ipl = ipl;     
      m_sensorID = m_detector.GetDet(m_ipl).GetDAQID();
      
      if ( std::find(m_filterIDs.begin(), m_filterIDs.end(), m_sensorID) == m_filterIDs.end() ) {
        streamlog_out(MESSAGE2) << " Do not create noise hits on sensorID: "  << m_sensorID << std::endl;
        continue;
      }
      
      // Average number of noise pixels
      double meanNoisePixels = m_noiseFraction * m_detector.GetDet(m_ipl).GetNColumns() * m_detector.GetDet(m_ipl).GetNRows();
            
      // Total number of noise pixels has Poison distribution
      int fractionPixels = gRandom->Poisson(meanNoisePixels);  

      // Generate noise digits
      for (int iNoisePixel=0; iNoisePixel<fractionPixels; iNoisePixel++) {
                    
        int iU  = int(gRandom->Uniform( m_detector.GetDet(m_ipl).GetNColumns() ));
        int iV  = int(gRandom->Uniform( m_detector.GetDet(m_ipl).GetNRows() ));
                  
        // Describe pixel by unique ID
        int uniqPixelID  = m_detector.GetDet(m_ipl).encodePixelID(iV, iU);     
                
        // Find if pixel doesn't already have some signal+noise or just noise
        if( !(digitsMap[m_sensorID].find(uniqPixelID)!=digitsMap[m_sensorID].end()) ) {
                     
          // Create a noise charge  
          double charge = m_zsThreshold; 
                     
          // Create new noise digit      
          digit = new Digit;        		       	
          digit->cellIDU  = iU;
          digit->cellIDV  = iV;
          digit->charge   = charge;
                            
          // Record it
          digitsMap[m_sensorID][uniqPixelID] = digit;
                  
        }
      } // For noise pixels   
    }
  }
  
  //
  // Method returning collected charge in ADC units
  //
  int SiPixDigitizer::getInADCUnits(double charge) {
    static double ADCunit = m_ADCRange/int(pow(2, m_ADCBits));
    static int ADCmax  = (int) pow(2, m_ADCBits) - 1.0;
    int val = int(charge/ADCunit);
    
    // lower end of dynamic range
    if (val <= 0) return 0;
    // upper end of dynamic range
    if (val > ADCmax) return ADCmax;
    return val;
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






