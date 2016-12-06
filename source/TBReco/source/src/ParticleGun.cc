// ParticleGun
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// user includes
#include "ParticleGun.h"


#include <CLHEP/Random/RandGamma.h>

// C++ includes
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <Exceptions.h>
#include <IMPL/LCFlagImpl.h>

#include <marlin/Global.h>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace CLHEP;

namespace depfet {

  //
  // Instantiate this object
  //
  ParticleGun aParticleGun ;


  //
  // Constructor
  //
  ParticleGun::ParticleGun() : Processor("ParticleGun")
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
                                m_GunTimeWindow,  static_cast < int > (0.0001)); 
   
    registerProcessorParameter ("BetheHeitlerSmearing", "Thickness of material before telescope for Bremsstrahlung [X/X0]",
                                m_GunBetheHeitlerT0,  static_cast < double > (0.0)); 
   
                                 
  }


void ParticleGun::init () {
  
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
void ParticleGun::processRunHeader(LCRunHeader * run)
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
void ParticleGun::processEvent(LCEvent * evt)
{
    
  //////////////////////////////////////////////////////////////////////  
  // Process next event
  ++_nEvt;
      
  // Create output MCParticle collection
  LCCollectionVec * outputCollection = new LCCollectionVec(LCIO::TRACK);
    
  // Set flag for storing track hits in track collection
  LCFlagImpl flag(outputCollection->getFlag());
  flag.setBit( LCIO::TRBIT_HITS );
  outputCollection->setFlag(flag.getFlag());
  
        
  evt->addCollection(outputCollection, m_MCParticleCollectionName); 

}

//
// Method called after each event to check the data processed
//
void ParticleGun::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void ParticleGun::end()
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
void ParticleGun::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "BeamEnergyCorrector Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}



} // Namespace


