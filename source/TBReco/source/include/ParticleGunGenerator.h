// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    ParticleGun2 - Marlin Processor                                                  //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef ParticleGun2_H
#define ParticleGun2_H 1


// TBTools includes 
#include "TBDetector.h"

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <marlin/ProcessorMgr.h>
#include <marlin/Exceptions.h>

// lcio includes <.h>

// system includes <>
#include <string>
#include <vector>

namespace depfet
{
  
  /** The ParticleGun2 Processor
   * The processor provides a particle gun for simulation of 
   * a directed particle beam.
   * 
   * Author: Benjamin Schwenker, Göttingen University 
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de> 
   */
  
  class ParticleGun2 : public marlin::Processor
  {
   public:
     
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new ParticleGun2 ; }
    
    //!Constructor - set processor description and register processor parameters
    ParticleGun2();
    
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
     
    //! Output track collection name
    std::string m_MCParticleCollectionName;
     
    //! Particle gun
    double m_GunXPosition;
    double m_GunYPosition;
    double m_GunZPosition; 
    double m_GunRotX;
    double m_GunRotY;
    double m_GunRotZ;
    double m_GunSpotSizeX;
    double m_GunSpotSizeY;
    double m_GunDivergenceX;
    double m_GunDivergenceY;
    double m_GunCorrelationX; 
    double m_GunCorrelationY;    
    double m_GunBeamIntensity;   
    double m_GunTimeWindow; 
    double m_GunBetheHeitlerT0;
    double m_GunMomentum;
    double m_GunMass;
    double m_GunCharge;
    
   private: 

    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
             
    
  };
 
}            
#endif
