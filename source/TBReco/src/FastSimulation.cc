// FastSimulation
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// user includes
#include "FastSimulation.h"

// C++ includes
#include <iostream>
#include <iomanip>

// Include LCIO classes
#include <lcio.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/MCParticleImpl.h> 
#include <IMPL/SimTrackerHitImpl.h> 
#include <IMPL/LCFlagImpl.h>

// ROOT includes
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>


// TBTools includes 
#include "MaterialEffect.h"
#include "HelixTrackModel.h"
#include "StraightLineTrackModel.h"
#include "TBDetector.h"

#include <Eigen/Core>
using Eigen::Vector3d;
using Eigen::Matrix;

// Used namespaces
using namespace std; 
using namespace lcio ;
using namespace marlin ;
using namespace std::string_literals;

namespace depfet {

  //
  // Instantiate this object
  //
  FastSimulation aFastSimulation ;


  //
  // Constructor
  //
  FastSimulation::FastSimulation() : Processor("FastSimulation")
  {
   
    // Processor description
    _description = "Fast simulation processor for propagation of beam particels trough a test beam telescope";
   

    //
    // Input collections  
    registerInputCollection (LCIO::MCPARTICLE, "MCParticleCollectionName",
                            "Collection name for MCParticles",
                            m_MCParticleCollectionName, string ("MCParticles") );

    //
    // Output collections  
    registerOutputCollection(LCIO::SIMTRACKERHIT,"SimTrackerHitCollectionName",
                             "Collection name for SimTrackerHits",
                             m_SimTrackerHitCollectionName, string ("SimTrackerHits"));
    
    registerProcessorParameter ("ScatterModel", "Choose model for multiple scattering: Highland(0), SumOfSingleScatterings(1)",
                                m_scatterModel,  static_cast < int > (0));
    
    registerProcessorParameter ("DoEnergyLossStraggling", "Flag (true/false) for simulating energy loss straggling",
                                m_doEnergyLossStraggling,  static_cast < bool > (false));
    
    registerProcessorParameter ("DoFractionalBetheHeitlerEnergyLoss", "Flag (true/false) for simulating fractional Bethe Heitler energy loss",
                                m_doFractionalBetheHeitlerEnergyLoss,  static_cast < bool > (false));
    
                                 
  }


  void FastSimulation::init () {
  
    // Initialize variables
    m_nRun = 0 ;
    m_nEvt = 0 ;
  
    // Print set parameters
    printProcessorParams();
  
    // CPU time start
    m_timeCPU = clock()/1000;
  }

  //
  // Method called for each run
  //
  void FastSimulation::processRunHeader(LCRunHeader * run)
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
  void FastSimulation::processEvent(LCEvent * evt)
  {

    // Print event number
    if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                  << (evt->getEventNumber())
                                                                  << std::endl << std::endl;
    
    //////////////////////////////////////////////////////////////////////  
    // Process next event
    ++m_nEvt;
    
    //
    // Open MCParticle collection
    LCCollectionVec* mcVec = nullptr;
    try {
      mcVec = dynamic_cast < LCCollectionVec * >  ( evt->getCollection(m_MCParticleCollectionName) );
    } catch (DataNotAvailableException& e) {
      throw SkipEventException(this);
    }
      
    // Init the track model for simulation
    GenericTrackModel* TrackModel; 
      
    Vector3d field;
    field<< TBDetector::GetInstance().GetBx(), TBDetector::GetInstance().GetBy(), TBDetector::GetInstance().GetBz();
      
    if ( field.norm() == 0 ) {
      TrackModel = new StraightLineTrackModel();    
    } else {
      TrackModel = new HelixTrackModel(field);  
    }   

    int nSensors = TBDetector::GetInstance().GetNSensors();
    
    // Create SimTrackerHit collection  
    LCCollectionVec * simHitVec = new LCCollectionVec(LCIO::SIMTRACKERHIT) ;
    
    // Create cellID encoder
    CellIDEncoder<SimTrackerHitImpl> cellIDEnc( "sensorID:6,isEntry:2,isExit:2" , simHitVec ) ;
    
    // Needed to store momentum to lcio file 
    LCFlagImpl aFlag(0);
    aFlag.setBit( LCIO::THBIT_MOMENTUM );
    simHitVec->setFlag( aFlag.getFlag() );
     
    // Loop on all MCParticles
    for (unsigned int iMC = 0; iMC < mcVec->size(); iMC++) 
    { 
         
      // Read back MCParticle 
      MCParticle* mcp = dynamic_cast<MCParticle* > ( mcVec->getElementAt(iMC) )  ;   
      
      streamlog_out(MESSAGE1) << "Processing generator particle " << mcp->getGeneratorStatus() 
                              << " at time[s] " << std::setprecision(10) << mcp->getTime() << std::setprecision(0) << endl; 
      
      
      // Create a complete trajectory state for the MCParticle needed to 
      // propagate it through the telescope 
      
      ReferenceFrame frame;
      
      // Origin of frame coincides with particle vertex


      Vector3d origin;
      origin<< mcp->getVertex()[0], mcp->getVertex()[1], mcp->getVertex()[2];

      frame.SetPosition(origin); 
      
      // Track state defined at reference frame
      double p1 = mcp->getMomentum()[0];
      double p2 = mcp->getMomentum()[1];
      double p3 = mcp->getMomentum()[2];  
     
      TrackState state; // du/dw,  dv/dw,  u,   v, q/p
      state<< p1/p3 , p2/p3, 0, 0,  mcp->getCharge()/std::sqrt(p1*p1 + p2*p2 + p3*p3);

        
      
      // Enter main loop to propagate particle through telescope
      int ipl = -1;
      
      // Average mometum before and after energy loss. 
      double weighted_mean_mom = 0;

      do  { 
        
        streamlog_out(MESSAGE1) << " Particle is at plane " << ipl << " with state " << state << endl; 
                                
        // Scatter at a material or sensor plane 
        
        if (ipl >= 0 ) { 
                 
          const Det& current_det = TBDetector::Get(ipl);
          
          // Local track state on sensor == defined as state before scattering and energy loss 
          // --------------------------- 
          double dudw = state[0];
          double dvdw = state[1];
          double u = state[2];
          double v = state[3];
          double mom = mcp->getCharge()/state[4];
          double l0 = current_det.GetThickness(u,v)*std::sqrt(1 + dudw*dudw + dvdw*dvdw);  
          
          // Simulate energy loss in material layer
          double eDep = materialeffect::GetMostProbableEnergyLossInSilicon(state, l0, mcp->getMass(), mcp->getCharge()); 
          if ( m_doEnergyLossStraggling ) {
            double lambda = gRandom->Landau();
            eDep = materialeffect::SimulateEnergyLossInSilicon(state, l0, mcp->getMass(), mcp->getCharge(), lambda); 
          }          
          
          // Create a new LCIO SimTracker hit
          SimTrackerHitImpl * simHit = new SimTrackerHitImpl;
          simHit->setdEdx(eDep); 
          simHit->setTime(mcp->getTime());
          simHit->setPathLength(l0); 
          // Set local hit position
          double hitPos[3] =  {u , v , 0} ;
          simHit->setPosition(hitPos);  
          // Set local hit momentum 
          float hitMom[3] = { float(dudw*mom/std::sqrt(dudw*dudw + dvdw*dvdw +1)), float(dvdw*mom/std::sqrt(dudw*dudw + dvdw*dvdw +1)) , float(1.0*mom/std::sqrt(dudw*dudw + dvdw*dvdw +1))};
          simHit->setMomentum(hitMom);
          
          // Set CellID
          cellIDEnc["sensorID"] = current_det.GetSensorID();  
          cellIDEnc["isEntry"] = 0;
          cellIDEnc["isExit"] = 0;
          cellIDEnc.setCellID(simHit);
          
          // Set particle that has interacted in the detector
          simHit->setMCParticle( mcp );
           
          // Put the hit into LCIO collection
          simHitVec->push_back( simHit );
           
          // Simulate energy loss by bremsstrahlung (Bethe Heitler theory)
		  // A weighted mean approach is used similar to the momentum calculation in
		  // the X/X0 imaging and calibration scripts.
		  // The momentum during the multiple scattering process is calculated as
		  // a weighted mean with the weight parameter epsilon, which is a steering parameter
		  // of this processor. 
		  // First term of the weighted mean
          weighted_mean_mom = (1.0-materialeffect::Epsilon_Weightfactor)*mcp->getCharge()/state[4];

   		  streamlog_out ( MESSAGE1 ) << endl;
   		  streamlog_out ( MESSAGE1 ) << "Materialeffects of plane " << ipl << endl;
   		  streamlog_out ( MESSAGE1 ) << "Momentum before transition: " << mcp->getCharge()/state[4] << " GeV" << endl;
   		  streamlog_out ( MESSAGE1 ) << "First term of weighted mean: " << weighted_mean_mom << " GeV" << endl;
		  
          if (m_doFractionalBetheHeitlerEnergyLoss) {
            double t = l0/current_det.GetRadLength(u,v);
            double rndm = gRandom->Rndm(1);
            materialeffect::SimulateBetherHeitlerEnergyLoss(state, t, mcp->getMass(), mcp->getCharge(), rndm);    
          }
          // Second term of the weighted mean
          weighted_mean_mom += materialeffect::Epsilon_Weightfactor*mcp->getCharge()/state[4];

   		  streamlog_out ( MESSAGE1 ) << "Momentum after transition: " << mcp->getCharge()/state[4] << " GeV" << endl;
   		  streamlog_out ( MESSAGE1 ) << "Second term of weighted mean: " << materialeffect::Epsilon_Weightfactor*mcp->getCharge()/state[4] << " GeV" << endl;
   		  streamlog_out ( MESSAGE1 ) << "Overall weighted mean: " << weighted_mean_mom << " GeV" << endl;
           
          if(  m_scatterModel==0  ) { 
            // Highland model scattering
            double theta2 = materialeffect::GetScatterTheta2(weighted_mean_mom, l0, current_det.GetRadLength(u,v), mcp->getMass(), mcp->getCharge() );      
            double kink_u = gRandom->Gaus(0, TMath::Sqrt( theta2 ));
            double kink_v = gRandom->Gaus(0, TMath::Sqrt( theta2 ));      
            // Scatter track ('in' state -> 'out' state)
            materialeffect::ScatterTrack(state, kink_u, kink_v); 
          } else if ( m_scatterModel==1 ) {
            double kink_u = materialeffect::GetScatterKink_SC(l0, current_det.GetRadLength(u,v), current_det.GetAtomicNumber(u,v), current_det.GetAtomicMass(u,v), mcp->getMass(), mcp->getCharge(), weighted_mean_mom   );
            double kink_v = materialeffect::GetScatterKink_SC(l0, current_det.GetRadLength(u,v), current_det.GetAtomicNumber(u,v), current_det.GetAtomicMass(u,v), mcp->getMass(), mcp->getCharge(), weighted_mean_mom   );
            // Scatter track ('in' state -> 'out' state)
            materialeffect::ScatterTrack(state, kink_u, kink_v);     
          } 
          
         
        }
         
        // Propagate particle to next detector surface
        // ----------------------------------------
        
        int jpl = ipl+1;      
            
        // Exit condition   
        if (jpl == nSensors ) {break;}
        
        // Next sensor along beam line 
        const ReferenceFrame& next_frame = TBDetector::Get(jpl).GetNominal();
               
        // Check that track really intersects surface
        if (! TrackModel->CheckHitsSurface(state, frame, next_frame) ) {
          break;    
        } 
        
        // Now, compute the fligth length to next surface
        double length = TrackModel->GetSignedStepLength(state, frame, next_frame);
        
        // Extrapolate half step 
        TrackModel->Extrapolate(state, frame, 0.5*length); 
        
        // Scatter in air between detectors
        // ------------------------------- 
        
        // Simulate energy loss by bremsstrahlung (Bethe Heitler theory)
		// A weighted mean approach is used similar to the momentum calculation in
		// the X/X0 imaging and calibration scripts.
		// The momentum during the multiple scattering process is calculated as
		// a weighted mean with the weight parameter epsilon, which is a steering parameter
		// of this processor. 
		// First term of the weighted mean
        weighted_mean_mom = (1-materialeffect::Epsilon_Weightfactor)*mcp->getCharge()/state[4];

   		streamlog_out ( MESSAGE1 ) << endl;
   		streamlog_out ( MESSAGE1 ) << "Materialeffects of air plane " << endl;
   		streamlog_out ( MESSAGE1 ) << "Momentum before transition: " << mcp->getCharge()/state[4] << " GeV" << endl;
   		streamlog_out ( MESSAGE1 ) << "First term of weighted mean: " << weighted_mean_mom << " GeV" << endl;

        if (m_doFractionalBetheHeitlerEnergyLoss) {
          double t = length/materialeffect::X0_air;
          double rndm = gRandom->Rndm(1);
          materialeffect::SimulateBetherHeitlerEnergyLoss(state, t, mcp->getMass(), mcp->getCharge(), rndm); 
        } 
		// Second term of the weighted mean
        weighted_mean_mom += materialeffect::Epsilon_Weightfactor*mcp->getCharge()/state[4];

   		streamlog_out ( MESSAGE1 ) << "Momentum after transition: " << mcp->getCharge()/state[4] << " GeV" << endl;
   		streamlog_out ( MESSAGE1 ) << "Second term of weighted mean: " << materialeffect::Epsilon_Weightfactor*mcp->getCharge()/state[4] << " GeV" << endl;
   		streamlog_out ( MESSAGE1 ) << "Overall weighted mean: " << weighted_mean_mom << " GeV" << endl;
        
        if( m_scatterModel==0 ) { 
          // Highland model scattering     
          double theta2 = materialeffect::GetScatterTheta2(weighted_mean_mom, length, materialeffect::X0_air, mcp->getMass(), mcp->getCharge()) ;   
          double kink_u = gRandom->Gaus(0, TMath::Sqrt( theta2 ));
          double kink_v = gRandom->Gaus(0, TMath::Sqrt( theta2 )); 
          // Scatter track ('in' state -> 'out' state)
          materialeffect::ScatterTrack(state, kink_u, kink_v);     
        } else if ( m_scatterModel==1 ) {
          double kink_u = materialeffect::GetScatterKink_SC(length, materialeffect::X0_air, materialeffect::AtomicNumber_air, materialeffect::AtomicMass_air, mcp->getMass(), mcp->getCharge(), weighted_mean_mom   );
          double kink_v = materialeffect::GetScatterKink_SC(length, materialeffect::X0_air, materialeffect::AtomicNumber_air, materialeffect::AtomicMass_air, mcp->getMass(), mcp->getCharge(), weighted_mean_mom   );
          // Scatter track ('in' state -> 'out' state)
          materialeffect::ScatterTrack(state, kink_u, kink_v);     
        } 
          
        // Get state on next detector
        bool error = false; 
        TrackState next_state = TrackModel->Extrapolate(state, frame, next_frame, error);
        
        if (error) { // check for loopers
          break;      
        } 
        
        // Update track snapshot  
        ipl = jpl;  
        frame = next_frame;
        state = next_state;  
           
      } while ( ipl < nSensors ); 
        
    }
    
    delete TrackModel; 
        
    evt->addCollection(simHitVec, m_SimTrackerHitCollectionName); 

  }

  //
  // Method called after each event to check the data processed
  //
  void FastSimulation::check( LCEvent * ) {}

  //
  // Method called after all data processing
  //
  void FastSimulation::end()
  {
     
    streamlog_out ( MESSAGE3 ) << endl;
    streamlog_out ( MESSAGE3 ) << "Successfully finished" << endl;
    
    // CPU time end
    m_timeCPU = clock()/1000 - m_timeCPU;
   
    // Print message
    streamlog_out(MESSAGE3) << std::endl
                             << " "
                             << "Time per event: "
                             << std::setiosflags(std::ios::fixed | std::ios::internal )
                             << std::setprecision(3)
                             << m_timeCPU/m_nEvt
                             << " ms"
                             << std::endl
                             << std::setprecision(3)
                             << std::endl
                             << " "
                             << "Processor succesfully finished!"
                             << std::endl;

 
  }


  //
  // Method printing processor parameters
  //
  void FastSimulation::printProcessorParams() const
  {

    streamlog_out(MESSAGE3)  << std::endl
                              << " "
                              << "FastSimulation development Version, be carefull!!"
                              << " "
                              << std::endl  << std::endl;   


  }



} // Namespace


