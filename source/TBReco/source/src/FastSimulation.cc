// FastSimulation
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// user includes
#include "FastSimulation.h"

// C++ includes
#include <iostream>

// Include LCIO classes
#include <lcio.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/MCParticleImpl.h> 
#include <IMPL/SimTrackerHitImpl.h> 

// ROOT includes
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>

// TBTools includes 
#include "MaterialEffect.h"
#include "HelixTrackModel.h"
#include "StraightLineTrackModel.h"

// CLHEP includes
#include "CLHEP/Matrix/Vector.h"
//#include "CLHEP/Vector/Rotation.h"
//#include "CLHEP/Vector/ThreeVector.h"


// Used namespaces
using namespace std; 
using namespace lcio ;
using namespace marlin ;
using namespace CLHEP;

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
    
    registerProcessorParameter ("ScatterModel", "Choose model for multiple scattering: Highland(0)",
                                m_scatterModel,  static_cast < int > (0));
   
    registerProcessorParameter ("ElossModel", "Choose model for energy loss: G4UniversalFluctuation(0)",
                                m_eLossModel,  static_cast < int > (0));
   
    
                                 
  }


  void FastSimulation::init () {
  
    // Initialize variables
    m_nRun = 0 ;
    m_nEvt = 0 ;
  
    // Print set parameters
    printProcessorParams();
  
    // CPU time start
    m_timeCPU = clock()/1000;

    // Read detector constants from gear file
    m_detector.ReadGearConfiguration();  
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
      
    HepVector field(3,0);
    field[0] = m_detector.GetBx();
    field[1] = m_detector.GetBy();
    field[2] = m_detector.GetBz(); 
      
    if ( field.norm() == 0 ) {
      TrackModel = new StraightLineTrackModel();    
    } else {
      TrackModel = new HelixTrackModel(field);  
    }   
    
    // Create SimTrackerHit collection  
    LCCollectionVec * simHitVec = new LCCollectionVec(LCIO::SIMTRACKERHIT) ;
     
    // Create cellID encoder
    CellIDEncoder<SimTrackerHitImpl> cellIDEnc( "layer:7,ladder:5,sensor:4,isEntry:2,isExit:2" , simHitVec ) ;
         
    // Loop on all MCParticles
    for (unsigned int iMC = 0; iMC < mcVec->size(); iMC++) 
    { 
         
      // Read back MCParticle 
      MCParticle* mcp = dynamic_cast<MCParticle* > ( mcVec->getElementAt(iMC) )  ;   
      
      streamlog_out(MESSAGE1) << "Processing generator particle " << mcp->getGeneratorStatus() 
                              << " at time " << mcp->getTime() << endl; 
      
      
      // Create a complete trajectory state for the MCParticle needed to 
      // propagate it through the telescope 
      
      ReferenceFrame frame;
      
      // Origin of frame coincides with particle vertex
      HepVector origin(3);
      origin[0] = mcp->getVertex()[0];
      origin[1] = mcp->getVertex()[1];
      origin[2] = mcp->getVertex()[2];
      frame.SetPosition(origin); 
      
      // Track state defined at reference frame
      HepMatrix state(5,1); 
        
      double p1 = mcp->getMomentum()[0];
      double p2 = mcp->getMomentum()[1];
      double p3 = mcp->getMomentum()[2];  
      double momentum = std::sqrt(p1*p1 + p2*p2 + p3*p3);
      
      state[0][0] = p1/p3;                              // du/dw
      state[1][0] = p2/p3;                              // dv/dw
      state[2][0] = 0;                                  // u 
      state[3][0] = 0;                                  // v  
      state[4][0] = mcp->getCharge()/momentum;          // q/p
        
      
      // Enter main loop to propagate particle through telescope
      int ipl = -1;
      
      do  { 
        
        streamlog_out(MESSAGE1) << " Particle is at plane " << ipl << " with state " << state << endl; 
                                
        // Scatter at a material or sensor plane 
        
        if (ipl >= 0 ) { 
                 
          Det& current_det = m_detector.GetDet(ipl);
          
          // Local track state on sensor
          // --------------------------- 
          double dudw = state[0][0];
          double dvdw = state[1][0];
          double u = state[2][0]; 
          double v = state[3][0]; 
          double mom = mcp->getCharge()/state[4][0]; 
          double l0 = current_det.GetThickness(u,v)*std::sqrt(1 + dudw*dudw + dvdw*dvdw);  
          double eDep = l0*82*3.64*1e-6;
             
          // Create a new LCIO SimTracker hit
          SimTrackerHitImpl * simHit = new SimTrackerHitImpl;
          simHit->setdEdx(eDep); 
          simHit->setTime(mcp->getTime());
          simHit->setPathLength(l0); 
          // Set local hit position
          double hitPos[3] =  {u , v , 0} ;
          simHit->setPosition(hitPos);  
          // Set local hit momentum 
          float hitMom[3] = { dudw*mom/std::sqrt(dudw*dudw + dvdw*dvdw +1), dvdw*mom/std::sqrt(dudw*dudw + dvdw*dvdw +1) , 1.0*mom/std::sqrt(dudw*dudw + dvdw*dvdw +1) };
          simHit->setMomentum(hitMom);
          
          // Set CellID
          cellIDEnc["layer"]  = ipl;
          cellIDEnc["ladder"] = ipl;
          cellIDEnc["sensor"] = ipl;
          cellIDEnc["isEntry"] = 0;
          cellIDEnc["isExit"] = 0;
          cellIDEnc.setCellID(simHit);
          
          // Set particle that has interacted in the detector
          simHit->setMCParticle( mcp );
    
          // Put the hit into LCIO collection
          simHitVec->push_back( simHit );
           
          if(  m_scatterModel==0  ) { 
            // Highland model scattering
            double theta2 = materialeffect::GetScatterTheta2(l0, current_det.GetRadLength(u,v), mcp->getMass(), mcp->getCharge(), mom );      
            double kink_u = gRandom->Gaus(0, TMath::Sqrt( theta2 ));
            double kink_v = gRandom->Gaus(0, TMath::Sqrt( theta2 ));      
            // Scatter track ('in' state -> 'out' state)
            materialeffect::ScatterTrack(state, kink_u, kink_v); 
          }
        }
         
        // Propagate particle to next detector surface
        // ----------------------------------------
        
        int jpl = ipl+1;      
            
        // Exit condition   
        if (jpl == m_detector.GetNSensors() ) {break;}
        
        // Next sensor along beam line 
        ReferenceFrame next_frame = m_detector.GetDet(jpl).GetNominal();
               
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
        
        if( m_scatterModel==0 ) { 
          // Highland model scattering     
          double theta2 = materialeffect::GetScatterTheta2(length, materialeffect::X0_air, mcp->getMass(), mcp->getCharge(), momentum ) ;   
          double kink_u = gRandom->Gaus(0, TMath::Sqrt( theta2 ));
          double kink_v = gRandom->Gaus(0, TMath::Sqrt( theta2 )); 
          // Scatter track ('in' state -> 'out' state)
          materialeffect::ScatterTrack(state, kink_u, kink_v);     
        } 
          
        // Get state on next detector
        bool error = false; 
        HepMatrix next_state = TrackModel->Extrapolate(state, frame, next_frame, error);  
        
        if (error) { // check for loopers
          break;      
        } 
        
        // Update track snapshot  
        ipl = jpl;  
        frame = next_frame;
        state = next_state;  
           
      } while ( ipl < m_detector.GetNSensors() ); 
        
    }
    
    delete TrackModel; 
        
    evt->addCollection(simHitVec, m_SimTrackerHitCollectionName); 

  }

  //
  // Method called after each event to check the data processed
  //
  void FastSimulation::check( LCEvent * evt ) {}

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
                              << "BeamEnergyCorrector Development Version, be carefull!!"
                              << " "
                              << std::endl  << std::endl;   


  }



} // Namespace


