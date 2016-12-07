// ParticleGunGenerator
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// user includes
#include "ParticleGunGenerator.h"

// C++ includes
#include <iostream>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCCollectionVec.h>
#include "IMPL/MCParticleImpl.h" 

// ROOT includes
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>



// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;


namespace depfet {

  //
  // Instantiate this object
  //
  ParticleGunGenerator aParticleGunGenerator ;


  //
  // Constructor
  //
  ParticleGunGenerator::ParticleGunGenerator() : Processor("ParticleGunGenerator")
  {
   
    // Processor description
    _description = "Particle gun processor for simulation of a directed particle beam";
   
    //
    // Output collections  
    registerOutputCollection(LCIO::MCPARTICLE,"MCParticleCollectionName",
                             "Collection name for MCParticles",
                             m_MCParticleCollectionName, string ("MCParticles"));
    
    registerProcessorParameter ("ParticleMass", "Particle mass [GeV]",
                                m_GunMass,  static_cast < double > (0.139));
    
    registerProcessorParameter ("ParticleCharge", "Particle charge [e]",
                                m_GunCharge,  static_cast < double > (+1));

    registerProcessorParameter ("PDG", "PDG code",
                                m_GunPDG,  static_cast < int > (22));
  
    registerProcessorParameter ("ParticleMomentum", "Particle momentum [GeV]",
                                m_GunMomentum,  static_cast < double > (4.0));

    registerProcessorParameter ("GunPositionX", "X position of particle gun [mm]",
                                m_GunPositionX,  static_cast < double > (0));
  
    registerProcessorParameter ("GunPositionY", "Y position of particle gun [mm]",
                                m_GunPositionY,  static_cast < double > (0));

    registerProcessorParameter ("GunPositionZ", "Z position of particle gun [mm]",
                                m_GunPositionZ,  static_cast < double > (-5000));
  
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
   
    
   
                                 
  }


  void ParticleGunGenerator::init () {
  
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
  void ParticleGunGenerator::processRunHeader(LCRunHeader * run)
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
  void ParticleGunGenerator::processEvent(LCEvent * evt)
  {
    
    //////////////////////////////////////////////////////////////////////  
    // Process next event
    ++_nEvt;
       
    // Create output MCParticle collection
    LCCollectionVec * mcVec = new LCCollectionVec( LCIO::MCPARTICLE );
    
    double time = 0; 
    int nParticle = 0; 

    while ( time < m_GunTimeWindow ) { 

      
      // Sample vertex and momentum from a directed particle beam.
      // 
      // Mean position and directions are always zero in collimator 
      // frame. 
      // 
      // Variables are decomposed into two uncorrelated pairs, namely 
      // (dx/dz,x) and (dy/dz,y). 
      //
      // We use a cholesky decomposition method to sample from the 
      // correlated pairs. 
      
      double X1 = gRandom->Gaus(0, 1); 
      double X2 = gRandom->Gaus(0, 1); 
      double Y1 = gRandom->Gaus(0, 1); 
      double Y2 = gRandom->Gaus(0, 1);
      
      double covX = m_GunCorrelationX*m_GunDivergenceX*m_GunSpotSizeX;
      double covY = m_GunCorrelationY*m_GunDivergenceY*m_GunSpotSizeY;  
      double DX = std::sqrt( std::pow(m_GunDivergenceX,2)*std::pow(m_GunSpotSizeX,2) - covX*covX );
      double DY = std::sqrt( std::pow(m_GunDivergenceY,2)*std::pow(m_GunSpotSizeY,2) - covY*covY );
      
      double P1 = m_GunDivergenceX*X1;                     // dx/dz
      double P2 = m_GunDivergenceY*Y1;                     // dy/dz
      double P3 = 1; 
      
      double Norm = m_GunMomentum/std::sqrt( P1*P1 + P2*P2 + 1 );
      P1 *= Norm; 
      P2 *= Norm; 
      P3 *= Norm; 
      
      double V1 = m_GunPositionX + (covX*X1 + DX*X2) / m_GunDivergenceX;    // x 
      double V2 = m_GunPositionY + (covY*Y1 + DY*Y2) / m_GunDivergenceY;    // y  
      double V3 = m_GunPositionZ; 
      
      //
	  //  Create a MCParticle and fill it from stdhep info
	  //
      MCParticleImpl* mcp = new MCParticleImpl();	
      //
	  //  PDGID
	  //
	  mcp->setPDG(m_GunPDG);
	  //
	  //  Momentum vector
	  //
	  float p0[3] = {P1,P2,P3};
	  mcp->setMomentum(p0);
	  //
	  //  Mass
	  //
	  mcp->setMass(m_GunMass);
      //
	  //  Charge
	  //
	  mcp->setCharge(m_GunCharge);
	  //
	  //  Vertex
	  // 
	  double v0[3] = {V1, V2, V3};
	  mcp->setVertex(v0);
	  //
	  //  Generator status
	  //
	  mcp->setGeneratorStatus(nParticle);
	  //
	  //  Simulator status 0 until simulator acts on it
	  //
	  mcp->setSimulatorStatus(nParticle);
	  //
	  //  Creation time (note the units)
	  // 
	  mcp->setTime(time);
	  //
      //  Add the particle to the collection vector
	  //
	  mcVec->push_back(mcp);
      
      //
      //  Advance time until the next beam particle
      // 
      time += gRandom->Exp(1.0/m_GunBeamIntensity);  
      nParticle += 1;  
    
    }
    
    //
    // Add MCParticle collection
    //     
    evt->addCollection(mcVec, m_MCParticleCollectionName); 

  }

  //
  // Method called after each event to check the data processed
  //
  void ParticleGunGenerator::check( LCEvent * evt ) {}

  //
  // Method called after all data processing
  //
  void ParticleGunGenerator::end()
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
  void ParticleGunGenerator::printProcessorParams() const
  {

    streamlog_out(MESSAGE3)  << std::endl
                              << " "
                              << "ParticleGunGenerator development Version, be carefull!!"
                              << " "
                              << std::endl  << std::endl;   


  }



} // Namespace


