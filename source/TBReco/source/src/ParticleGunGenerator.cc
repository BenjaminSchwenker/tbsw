// ParticleGunGenerator
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// user includes
#include "ParticleGunGenerator.h"

// C++ includes
#include <iostream>
#include <cmath>
#include <iomanip>


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
using namespace CLHEP;

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
                                m_ParticleMass,  static_cast < double > (0.139));
    
    registerProcessorParameter ("ParticleCharge", "Particle charge [e]",
                                m_ParticleCharge,  static_cast < double > (+1));

    registerProcessorParameter ("PDG", "PDG code",
                                m_ParticlePDG,  static_cast < int > (22));
  
    registerProcessorParameter ("BeamMomentum", "Beam momentum [GeV]",
                                m_BeamMomentum,  static_cast < double > (4.0));

    registerProcessorParameter ("BeamVertexX", "Beam vertex position in X direction [mm]",
                                m_BeamVertexX,  static_cast < double > (0));
  
    registerProcessorParameter ("BeamVertexY", "Beam vertex position in Y direction [mm]",
                                m_BeamVertexY,  static_cast < double > (0));

    registerProcessorParameter ("BeamVertexZ", "Beam vertex position in Z direction [mm]",
                                m_BeamVertexZ,  static_cast < double > (-5000));
    
    registerProcessorParameter ("BeamVertexXSigma", "Beam vertex sigma in X direction [mm]",
                                m_BeamVertexXSigma,  static_cast < double > (1)); 
    
    registerProcessorParameter ("BeamVertexYSigma", "Beam vertex sigma in Y direction [mm]",
                                m_BeamVertexYSigma,  static_cast < double > (1)); 
     
    registerProcessorParameter ("BeamSlopeXSigma", "Beam slope sigma in XZ plane [rad]",
                                m_BeamSlopeXSigma,  static_cast < double > (0.0001)); 
    
    registerProcessorParameter ("BeamSlopeYSigma", "Beam slope sigma in YZ plane [rad]",
                                m_BeamSlopeYSigma,  static_cast < double > (0.0001)); 
    
    registerProcessorParameter ("BeamMomentumSigma", "Beam momentum sigma [GeV]",
                                m_BeamMomentumSigma,  static_cast < double > (0.01)); 
    
    registerProcessorParameter ("CorrelationVertexXvsSlopeX", "Beam correlation factor between vertex position and slope",
                                m_BeamCorrelationVertexXvsSlopeX,  static_cast < double > (0.0)); 
    
    registerProcessorParameter ("CorrelationVertexYvsSlopeY", "Beam correlation factor between vertex position and slope",
                                m_BeamCorrelationVertexYvsSlopeY,  static_cast < double > (0.0)); 
    
    registerProcessorParameter ("CorrelationVertexXvsMomentum", "Beam correlation factor between vertex position and momentum",
                                m_BeamCorrelationVertexXvsMomentum,  static_cast < double > (0.0)); 
    
    registerProcessorParameter ("CorrelationVertexYvsMomentum", "Beam correlation factor between vertex position and momentum",
                                m_BeamCorrelationVertexYvsMomentum,  static_cast < double > (0.0)); 
    
    registerProcessorParameter ("BeamIntensity", "Number of beam particles per second",
                                m_BeamIntensity,  static_cast < double > (10000)); 
    
    registerProcessorParameter ("BeamTimeWindow", "A simulated event contains one particle at t=0 and extends for given time window in seconds",
                                m_BeamTimeWindow,  static_cast < double > (0.0001)); 
                                 
  }
  
  
  void ParticleGunGenerator::init () {
  
    // Initialize variables
    _nRun = 0 ;
    _nEvt = 0 ;
    _nParticles = 0;
  
    // Print set parameters
    printProcessorParams();
  
    // CPU time start
    _timeCPU = clock()/1000;

    if ( m_BeamSlopeXSigma <= 0 ) m_BeamSlopeXSigma = 0.00001;     // in rad 
    if ( m_BeamSlopeYSigma <= 0 ) m_BeamSlopeYSigma = 0.00001;     // in rad 
    if ( m_BeamVertexXSigma <= 0 ) m_BeamVertexXSigma = 0.001;     // in mm
    if ( m_BeamVertexYSigma <= 0 ) m_BeamVertexYSigma = 0.001;     // in mm
    if ( m_BeamMomentumSigma <= 0 ) m_BeamMomentumSigma = 0.00001; // in GeV
    
    // Construct beam covariance matrix S 
    //
    // The particle state is defined relativ to the Z= m_GunPositionZ
    // plane. The 5 dimensional state vector of the beam particle is  
    // (dx/dz,dy/dz,x,y,p) in global XYZ coordinates. 
    HepSymMatrix S(5,0);
    S[0][0] = std::pow(m_BeamSlopeXSigma,2);
    S[1][1] = std::pow(m_BeamSlopeYSigma,2);
    S[2][2] = std::pow(m_BeamVertexXSigma,2);
    S[3][3] = std::pow(m_BeamVertexYSigma,2);
    S[4][4] = std::pow(m_BeamMomentumSigma,2);    
    S[0][2] = m_BeamCorrelationVertexXvsSlopeX*m_BeamSlopeXSigma*m_BeamVertexXSigma;
    S[2][4] = m_BeamCorrelationVertexXvsMomentum*m_BeamVertexXSigma*m_BeamMomentumSigma;
    S[1][3] = m_BeamCorrelationVertexYvsSlopeY*m_BeamVertexYSigma*m_BeamSlopeYSigma;
    S[3][4] = m_BeamCorrelationVertexYvsMomentum*m_BeamVertexYSigma*m_BeamMomentumSigma;     
    
    // Construct average track state of beam particle 
    m_Mean = HepVector(5);
    m_Mean(1) = 0;
    m_Mean(2) = 0;
    m_Mean(3) = m_BeamVertexX;
    m_Mean(4) = m_BeamVertexY;
    m_Mean(5) = m_BeamMomentum;
    
    // Decompose S and store results to 
    // later generate random deviates    
    m_U  = HepMatrix(5,1);
    m_Sigmas = HepVector(5);
    
    HepSymMatrix tempS ( S ); 
    
    m_U = diagonalize ( &tempS );               // S = U Sdiag U.T()
    HepSymMatrix D = S.similarityT(m_U);        // D = U.T() S U = Sdiag
    for (int i = 1; i <= S.num_row(); i++) {
      double s2 = D(i,i);
      if ( s2 > 0 ) {
	    m_Sigmas(i) = std::sqrt ( s2 );
      } else {
        std::cerr << "In ParticelGunGenerator: " <<
                    "      Beam covariance matrix is not positive definite.  Eigenvalues are:\n";
        for (int ixx = 1; ixx <= S.num_row(); ixx++) {
          std::cerr << "      " << D(ixx,ixx) << std::endl;
        }
        std::cerr << "---Exiting to System\n";
        exit(1);
      }
    } 
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
    
    // Print event number
    if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                  << (evt->getEventNumber())
                                                                  << std::endl << std::endl;
    
    //////////////////////////////////////////////////////////////////////  
    // Process next event
    ++_nEvt;
       
    // Create output MCParticle collection
    LCCollectionVec * mcVec = new LCCollectionVec( LCIO::MCPARTICLE );
    
    double time = 0; 
    int nParticle = 0; 

    while ( time < m_BeamTimeWindow ) { 

      // Sample 5 dimensional state vector for new beam particle 
      HepVector state = deviates();
    
      // The sampling may return a negative momentum 
      if ( state(5) <= 0.001 ) state(5) = 0.001;

      double momentum = state(5);    // momentum in GeV        
      double tx = state(1);          // direction tangent dx/dz
      double ty = state(2);          // direction tangent dy/dz 
      
      double Norm = momentum/std::sqrt( tx*tx + ty*ty + 1 );
      double P1 = tx*Norm; 
      double P2 = ty*Norm; 
      double P3 = Norm; 
      
      double V1 = state(3);          // x vertex position  
      double V2 = state(4);          // y vertex position   
      double V3 = m_BeamVertexZ;     // z vertex position  
          
      //
	  //  Create a MCParticle and fill it from stdhep info
	  //
      MCParticleImpl* mcp = new MCParticleImpl();	
      //
	  //  PDGID
	  //
	  mcp->setPDG(m_ParticlePDG);
	  //
	  //  Momentum vector
	  //
	  double p0[3] = {P1,P2,P3};
	  mcp->setMomentum(p0);
	  //
	  //  Mass
	  //
	  mcp->setMass(m_ParticleMass);
      //
	  //  Charge
	  //
	  mcp->setCharge(m_ParticleCharge);
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
      time += gRandom->Exp(1.0/m_BeamIntensity);  
      nParticle += 1;  
    
    }

    //
    // Add MCParticle collection
    //     
    evt->addCollection(mcVec, m_MCParticleCollectionName); 
    
    // 
    // Increment counter for particles 
    // 
    _nParticles += nParticle;
    
    streamlog_out(MESSAGE2) << "Number of particels in this event is: " << nParticle << std::endl << std::endl;

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

    streamlog_out ( MESSAGE3 ) << "Number of simulated particles is " << _nParticles  << endl;
    
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
  
  HepVector ParticleGunGenerator::deviates( ) const
  {
    // Returns vector of gaussian randoms based on sigmas, rotated by U,
    // with means of 0. 
    HepVector v(5);  // The vector to be returned
    for ( int i = 1; i <= 5; i++ ) {
      v(i) = gRandom->Gaus(0, m_Sigmas(i));  
    }
    return m_Mean + m_U*v;
  } 
  
} // Namespace


