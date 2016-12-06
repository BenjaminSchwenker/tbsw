// FastSimulation
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// user includes
#include "FastSimulation.h"

#include "DEPFET.h" 
#include <CLHEP/Random/RandGamma.h>

// C++ includes
#include <iostream>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCCollectionVec.h>
#include "IMPL/MCParticleImpl.h" 



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
    _description = "Particle gun processor for simulation of a directed particle beam";
   
    //
    // Output collections  
    registerOutputCollection(LCIO::TRACK,"MCParticleCollectionName",
                             "Collection name for MCParticles",
                             m_MCParticleCollectionName, string ("MCParticles"));
    
    registerProcessorParameter ("ParticleMomentum", "Particle momentum [GeV]",
                                m_GunMomentum,  static_cast < double > (4.0));
   
    registerProcessorParameter ("ParticleMass", "Particle mass [GeV]",
                                m_GunMass,  static_cast < double > (0.139));
   
    registerProcessorParameter ("ParticleCharge", "Particle charge [e]",
                                m_GunCharge,  static_cast < double > (+1));
  
    registerProcessorParameter ("GunPositionX", "X position of particle gun [mm]",
                                m_GunXPosition,  static_cast < double > (0));
  
    registerProcessorParameter ("GunPositionY", "Y position of particle gun [mm]",
                                m_GunYPosition,  static_cast < double > (0));

    registerProcessorParameter ("GunPositionZ", "Z position of particle gun [mm]",
                                m_GunZPosition,  static_cast < double > (-5000));
  
    registerProcessorParameter ("GunRotationX", "X rotation of particle gun [rad]",
                                m_GunRotX,  static_cast < double > (0)); 
  
    registerProcessorParameter ("GunRotationY", "Y rotation of particle gun [rad]",
                                m_GunRotY,  static_cast < double > (0)); 

    registerProcessorParameter ("GunRotationZ", "Z rotation of particel gun [rad]",
                                m_GunRotZ,  static_cast < double > (0)); 
  
    registerProcessorParameter ("GunSpotSizeX", "Smearing of X vertex position at beam collimator [mm]",
                                m_GunSpotSizeX,  static_cast < double > (1)); 

    registerProcessorParameter ("GunSpotSizeY", "Smearing of Y vertex position at beam collimator [mm]",
                                m_GunSpotSizeY,  static_cast < double > (1)); 
  
    registerProcessorParameter ("GunDivergenceX", "RMS track slope in XZ plane [rad]",
                                m_GunDivergenceX,  static_cast < double > (0.0001)); 
  
    registerProcessorParameter ("GunDivergenceY", "RMS track slope in YZ plane [rad]",
                                m_GunDivergenceY,  static_cast < double > (0.0001)); 
  
    registerProcessorParameter ("GunCorrelationX", "Beam correlation coefficient X",
                                m_GunCorrelationX,  static_cast < double > (0.0)); 
  
    registerProcessorParameter ("GunCorralationY", "Beam correlation coefficient Y",
                                m_GunCorrelationY,  static_cast < double > (0.0)); 
   
    registerProcessorParameter ("GunIntensity", "Number of particles per second",
                                m_GunBeamIntensity,  static_cast < double > (10000)); 
  
    registerProcessorParameter ("GunTimeWindow", "A simulated event contains one particle at t=0 and extends for given time window in seconds",
                                m_GunTimeWindow,  static_cast < double > (0.0001)); 
   
    registerProcessorParameter ("BetheHeitlerSmearing", "Thickness of material before telescope for Bremsstrahlung [X/X0]",
                                m_GunBetheHeitlerT0,  static_cast < double > (0.0)); 
   
                                 
  }


  void FastSimulation::init () {
  
    // Initialize variables
    _nRun = 0 ;
    _nEvt = 0 ;
  
    // Print set parameters
    printProcessorParams();
  
    // CPU time start
    _timeCPU = clock()/1000;
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

    _nRun++ ;
  
  }

  //
  // Method called for each event
  //
  void FastSimulation::processEvent(LCEvent * evt)
  {
    
    //////////////////////////////////////////////////////////////////////  
    // Process next event
    ++_nEvt;
       
    // Create output MCParticle collection
    LCCollectionVec * mcVec = new LCCollectionVec( LCIO::MCPARTICLE );
    
    //
	//  Create a MCParticle and fill it from stdhep info
	//
    MCParticleImpl* mcp = new MCParticleImpl();
			
    //
	//  PDGID
	//
	mcp->setPDG(IDHEP);
	//
	//  Momentum vector
	//
	float p0[3] = {PHEP1,PHEP2,PHEP3};
	mcp->setMomentum(p0);
	//
	//  Mass
	//
	mcp->setMass(PHEP5);
	//
	//  Vertex
	// 
	double v0[3] = {VHEP1, VHEP2, VHEP3};
	mcp->setVertex(v0);
	//
	//  Generator status
	//
	mcp->setGeneratorStatus(ISTHEP);
	//
	//  Simulator status 0 until simulator acts on it
	//
	mcp->setSimulatorStatus(0);
	//
	//  Creation time (note the units)
	// 
	mcp->setTime(VHEP4/c_light);
	//
    //  Add the particle to the collection vector
	//
	mcVec->push_back(mcp);
    //
    // Add MCParticle collection
    // 
    lcEvt->addCollection( (LCCollection*) lcMCVec , "MCParticle" ) ;
     


    /*
    PHEP1 = mcp->getMomentum()[0];
    PHEP2 = mcp->getMomentum()[1];
    PHEP3 = mcp->getMomentum()[2];

    PHEP5 = mcp->getMass();

      std::cout << " HepLCIOInterface::GeneratePrimaryVertex: adding particle with pdg: " 
// 		<< mcp->getPDG() << " and energy: " << mcp->getEnergy()  
// 		<< " status : ["  << mcp->getGeneratorStatus() <<"]" << std::endl ;

     //      //use the vertex information from the MCParticles
//	  const double* mcpVertex = mcp->getVertex();
//      G4ThreeVector vertex_position(mcpVertex[0], mcpVertex[1], mcpVertex[2]);


    //	  // remove vertex information from the MCParticle
//	  MCParticleImpl* mcpImpl = dynamic_cast<MCParticleImpl*>(mcp);
//
//	  double mcpv[3] = {0.0, 0.0, 0.0};
//	  mcpImpl->setVertex(mcpv);
//	  mcpImpl->setTime(0.0);

   if( theParticle )
    mcp_impl->setCharge(theParticle->GetPDGCharge());
  else
    mcp_impl->setCharge(-1000);

    */

    
    /*
    // Particles are created in the plane of the particle
    // beam collimator opening 
    
    // First, we must construct of reference frame uvw plane for the
    // collimator opening 
    // The beam spot center is at u=v=0. The direction of the beam is 
    // the ew direction. 
    
    ReferenceFrame CollimatorFrame;
    
    // Position of center of collimator 
    HepVector CollimatorPosition(3);
    CollimatorPosition[0] = positionX; 
    CollimatorPosition[1] = positionY;
    CollimatorPosition[2] = positionZ;
    CollimatorFrame.SetPosition(CollimatorPosition); 
    
    // Rotation matrix of collimator plane  
    HepMatrix CollimatorRotation;
    double alpha = rotationX; 
    double beta = rotationY; 
    double gamma = rotationZ; 
    FillRotMatrixKarimaki(CollimatorRotation, alpha,beta,gamma);
    CollimatorFrame.SetRotation(CollimatorRotation);    
    
    // Next, we create a local track state in the collimator frame.
    HepMatrix GunState(5,1); 
    
    // Sample from a correlated particle position and direction 
    // variables 
    // 
    // Mean position and directions are always zero in collimator 
    // frame. 
    // 
    // Variables are decomposed into two uncorrelated pairs, namely 
    // (dx/dz,x) and (dy/dz,y). 
    //
    // We use a chelesky decomposition method to sample from the 
    // correlated pairs. 
    
    double X1 = myRng->Gaus(0, 1); 
    double X2 = myRng->Gaus(0, 1); 
    double Y1 = myRng->Gaus(0, 1); 
    double Y2 = myRng->Gaus(0, 1);

    double covX = correlationX*divergenceX*spotsizeX;
    double covY = correlationY*divergenceY*spotsizeY;  
    double DX = std::sqrt( std::pow(divergenceX,2)*std::pow(spotsizeX,2) - covX*covX );
    double DY = std::sqrt( std::pow(divergenceY,2)*std::pow(spotsizeY,2) - covY*covY );
    
    GunState[0][0] = divergenceX*X1;                     // dx/dz
    GunState[1][0] = divergenceY*Y1;                     // dy/dz
    GunState[2][0] = (covX*X1 + DX*X2) / divergenceX;    // x 
    GunState[3][0] = (covY*Y1 + DY*Y2) / divergenceY;    // y  
    GunState[4][0] = charge/mom;    // q/p
    
    // Finally, we put everything together :)
    TrackSnapShot MyTrack; 
    MyTrack.Surf = CollimatorFrame;
    MyTrack.State = GunState;
    */



    
        
    evt->addCollection(outputCollection, m_MCParticleCollectionName); 

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
    _timeCPU = clock()/1000 - _timeCPU;
   
    // Print message
    streamlog_out(MESSAGE3) << std::endl
                             << " "
                             << "Time per event: "
                             << std::setiosflags(std::ios::fixed | std::ios::internal )
                             << std::setprecision(3)
                             << _timeCPU/_nEvt
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


