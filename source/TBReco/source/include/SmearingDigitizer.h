// //////////////////////////////////////////////////////////////////////////  //
//                                                                             //
//    SmearingDigitizer - Marlin Processor                                        //
//                                                                             //
// //////////////////////////////////////////////////////////////////////////  //

#ifndef SmearingDigitizer_H
#define SmearingDigitizer_H 1

// Include TBTools header files
#include "TBDetector.h"

// Include basic C
#include <vector>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <marlin/Exceptions.h>


namespace depfet {

  
  /** The SmearingDigitizer Processor
   * Smearing digitizer creates TrackerHits by applying noise to SimTrackerHits
   * 
   * Author: Benjamin Schwenker, Göttingen University 
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de> 
   */
   
  class SmearingDigitizer : public marlin::Processor {
   
   public:
    
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new SmearingDigitizer ; }
     
    //!Constructor - set processor description and register processor parameters
    SmearingDigitizer();
    
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

   protected:

    //!Method printing processor parameters
    void printProcessorParams() const;

    // VARIABLES
    
    //! SimTrackerHit collection name
    std::string m_SimTrackerHitCollectionName;

    //! Output hit collection name
    std::string  m_hitCollectionName;
        
    // Readout/DAQ parameters - set by user
    std::vector<int >  m_filterIDs;              //!< Digitize only sensors in this list
    double m_sigmaU;                             //!< Sigma for covariance matrix [mm]
    double m_sigmaV;                             //!< Sigma for covariance matrix [mm]                                               
    double m_tanLorentzAngle;                    //!< Tangent of Lorentz angle
    bool   m_integrationWindow;                  //!< Use integration window?
    double m_startIntegration;                   //!< Start time of integration of the sensors in ns (everything before this value will not be digitized)
    double m_stopIntegration;                    //!< Stop time of integration of the sensors in ns(everything after this value will not be digitized)
    
   private:
    
    // Handle to detector data sheets 
    TBDetector m_detector;  
    
    double m_timeCPU; //!< CPU time
    int    m_nRun ;   //!< Run number
    int    m_nEvt ;   //!< Event number 
    
  }; // Class

} // Namespace

#endif // SmearingDigitizer_H



