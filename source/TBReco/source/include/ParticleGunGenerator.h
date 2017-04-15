// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    ParticleGunGenerator - Marlin Processor                                                  //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef ParticleGunGenerator_H
#define ParticleGunGenerator_H 1


// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <marlin/ProcessorMgr.h>
#include <marlin/Exceptions.h>

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

// lcio includes <.h>

// system includes <>
#include <string>
#include <vector>

namespace depfet
{
  
  /** The ParticleGunGenerator Processor
   * The processor provides a particle gun for simulation of 
   * a directed particle beam.
   * 
   * Author: Benjamin Schwenker, Göttingen University 
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de> 
   */
  
  class ParticleGunGenerator : public marlin::Processor
  {
   public:
     
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new ParticleGunGenerator ; }
    
    //!Constructor - set processor description and register processor parameters
    ParticleGunGenerator();
    
    //!Method called at the beginning of data processing - used for initialization
    virtual void init();
    
    //!Method called for each run - used for run header processing
    virtual void processRunHeader(LCRunHeader * run);
    
    //!Method called for each event - used for event data processing
    virtual void processEvent(LCEvent * evt);
    
    //!Method called after each event - used for data checking
    virtual void check(LCEvent * evt);
    
    //!Method called after all data processing
    virtual void end();
    
    //!Method printing processor parameters
    void printProcessorParams() const;
      
   protected:
    
    //!Method returns vector of gaussian randoms based on sigmas, rotated by U,
    // with means of 0. 
    CLHEP::HepVector deviates() const;
     
    //! Output MCParticle collection name
    std::string m_MCParticleCollectionName;
     
    //! Particle beam paramaters
    double m_BeamVertexX;
    double m_BeamVertexY;
    double m_BeamVertexZ;
    double m_BeamMomentum;
    double m_BeamVertexXSigma;
    double m_BeamVertexYSigma;
    double m_BeamSlopeXSigma;
    double m_BeamSlopeYSigma;
    double m_BeamMomentumSigma; 
    double m_BeamCorrelationVertexXvsSlopeX;
    double m_BeamCorrelationVertexYvsSlopeY;
    double m_BeamCorrelationVertexXvsMomentum;
    double m_BeamCorrelationVertexYvsMomentum;
    double m_BeamIntensity;   
    double m_BeamTimeWindow;
    double m_ParticleMass;
    double m_ParticleCharge;
    int m_ParticlePDG;
    
   private: 
    
    // Vector of sigmas and rotation matrix U for computing 
    // random deviates
    CLHEP::HepMatrix m_U;
    CLHEP::HepVector m_Sigmas;	
    CLHEP::HepVector m_Mean;
    
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
    int    _nParticles; //!< Number of particles  
            
  };
 
}            
#endif
