// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    FastSimulation - Marlin Processor                                                  //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef FastSimulation_H
#define FastSimulation_H 1

// Include TBTools header files
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
  
  /** The FastSimulation Processor
   * The processor tracks all particles found in a collection of type 
   * lcio::MCParticle trough a test beam telescope geometry and adds 
   * a new lcio::SimTrackerHits collection to the event. 
   * 
   * The processor uses the very same track extrapolation code as the 
   * track fitter. It is not a Geant4 based simulation.
   * 
   * Author: Benjamin Schwenker, Göttingen University 
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de> 
   */
  
  class FastSimulation : public marlin::Processor
  {
   public:
     
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new FastSimulation ; }
    
    //!Constructor - set processor description and register processor parameters
    FastSimulation();
    
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
     
    //! Input MCParticle collection name
    std::string m_MCParticleCollectionName;
     
    //! Output SimTrackerHit collection name
    std::string m_SimTrackerHitCollectionName;
    
    //! Alignment DB file name 
    std::string _alignmentDBFileName;
    
    // Choose model for multiple scattering ( Highland:0 )
    int m_scatterModel;  
 
    // Flag for simulating energy loss straggling from ionization 
    bool m_doEnergyLossStraggling;  
    
    // Flag for simulating fractional bethe heitler energy loss 
    bool m_doFractionalBetheHeitlerEnergyLoss; 
   
   private: 
    
    // Handle to detector data sheets 
    TBDetector m_detector;  
    
    double m_timeCPU; //!< CPU time
    int    m_nRun ;   //!< Run number
    int    m_nEvt ;   //!< Event number
             
    
  };
 
}            
#endif
