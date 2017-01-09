// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    TriggerGenerator - Marlin Processor                                                  //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef TriggerGenerator_H
#define TriggerGenerator_H 1

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
  
  class TrgScinti {
   	
   public:
    
    // Constructors 
    TrgScinti(int sensID, float minU, float maxU, float minV, float maxV) : 
      m_sensID(sensID), m_minU(minU), m_maxU(maxU), m_minV(minV), m_maxV(maxV) {}
        
    int GetDAQID(){ return m_sensID; };
  
    bool isPointInSensor( double u , double v); 

   private:
    
    // DAQ Id
    int m_sensID;
    // min/max U
    float m_minU; 
    float m_maxU;
    // min/max V
    float m_minV; 
    float m_maxV;
  };
  
  /** The TriggerGenerator Processor
   * The processor simulates a trigger system for a test beam 
   * experiment. The trigger signal is raised when the first 
   * MCParticle in an event generates a coincidence from one 
   * or more trigger scintillators. Up to four scintillators 
   * can be used. For simplicity, each scinitillators is attached
   * to a sensor plane and covers a rectangular region from 
   * (umin,vmin) to (umax,vmax) in local sensor coordinates. 
   * The scintillator gets actived, if a SimTrackerHit from 
   * the first MCParticel (at t=0) intercepts this region. 
   * 
   * Author: Benjamin Schwenker, Göttingen University 
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de> 
   */
  
  class TriggerGenerator : public marlin::Processor
  {
   public:
     
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new TriggerGenerator ; }
    
    //!Constructor - set processor description and register processor parameters
    TriggerGenerator();
    
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
     
    //! Output SimTrackerHit collection name
    std::string m_SimTrackerHitCollectionName;
    
    //! Every fakeTriggerPeriod event has a fake trigger
    int m_fakeTriggerPeriod;
 
    //! Scinti parameters: DAQ ID, Umin, Vmin, Umax, Vmax (leave empty to deactivate)
    std::vector<float >  _scintiNo1;
    std::vector<float >  _scintiNo2;
    std::vector<float >  _scintiNo3;
    std::vector<float >  _scintiNo4;
    
   private: 

    // Active scintis
    std::vector< TrgScinti > m_scintiVec;
    
    // Handle to detector data sheets 
    TBDetector m_detector;  
    
    double m_timeCPU; //!< CPU time
    int    m_nRun ;   //!< Run number
    int    m_nEvt ;   //!< Event number
             
  };
 
}            
#endif
