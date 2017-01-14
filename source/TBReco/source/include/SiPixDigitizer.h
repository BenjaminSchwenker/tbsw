// //////////////////////////////////////////////////////////////////////////  //
//                                                                             //
//    SiPixDigitizer - Marlin Processor                                        //
//                                                                             //
// //////////////////////////////////////////////////////////////////////////  //

#ifndef SiPixDigitizer_H
#define SiPixDigitizer_H 1

// Include TBTools header files
#include "TBDetector.h"

// Include basic C
#include <vector>

// Include CLHEP classes
#include <CLHEP/Vector/ThreeVector.h>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerDataImpl.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>


#include <TRandom.h>
#include <TRandom3.h>
#include <TMath.h>

namespace depfet {

  // Structures
  struct Digit
  {
    int cellIDV;
    int cellIDU;
    double charge;
  };

  struct IonisationPoint
  {
    CLHEP::Hep3Vector position;
    double eLoss;
  };

  struct SignalPoint
  {
    CLHEP::Hep3Vector position;
    CLHEP::Hep3Vector sigma;
    double charge;
  };

  struct SpacePoint
  {
    CLHEP::Hep3Vector position;
    CLHEP::Hep3Vector direction;  
  };


  // TypeDefs
  typedef std::vector<Digit*>                   DigitVec;
  typedef std::vector<IonisationPoint*>         IonisationPointVec;
  typedef std::vector<SignalPoint*>             SignalPointVec;
  typedef std::map<int, Digit*>                 DigitMap;
  typedef std::map<int, std::map<int, Digit*> > DigitsMap; // unique sensorID, unique pixelID, Digit

  /** The SiPixDigitizer Processor
   * A rather generic digitizer for silicon pixel detectors. 
   * 
   * Author: Benjamin Schwenker, Göttingen University 
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de> 
   */
   
  class SiPixDigitizer : public marlin::Processor {
   
   public:
    
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new SiPixDigitizer ; }
     
    //!Constructor - set processor description and register processor parameters
    SiPixDigitizer();
    
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

    // MAIN METHODS
    
    //!Method calculating ionisation points along the path from given SimTrackerHit
    void ProduceIonisationPoints(const SimTrackerHit * simTrkHit, IonisationPointVec & ionisationPoints);
    
    //!Method calculating signal points from ionisation points
    void ProduceSignalPoints(const IonisationPointVec & ionisationPoints, SignalPointVec & signalPoints);
    
    //!Method calculating digits from signal points - diffusion, Poisson and electronics effects added
    void ProduceSignalDigits(const SignalPointVec & signalPoints, DigitVec & digits);
    
    //!Method for calculating electronic effects on signal digits
    void ProduceNoiseEffects(DigitsMap & digitsMap);
    
    //!Method producing noise digits and updating digitsmap containing all digits
    void ProduceNoiseDigits(DigitsMap & digitsMap);
    
    //!Method producing digits output from all digits (noise + signal)
    void WriteDigitsToLCIO(LCCollectionVec * digitCol , const DigitsMap & digitsMap);
    
    
    // OTHER METHODS
    
    //!Method transforming hit to local coordinate system using SiPxlGeom interface
    void TransformToLocal(const SimTrackerHit * simTrkHit, SpacePoint & hitLocal);
    
    //!Method updating PXD map (containing all digits) using given vector of digits
    void UpdateDigitsMap(DigitsMap & digitsMap, DigitVec & digits);
     
    //!Method returning collected charge in ADC units
    int getInADCUnits(double charge);
    
    //!Method providing mathematical round
    int round(double num);

    //!Method printing processor parameters
    void printProcessorParams() const;

    // VARIABLES
    
    //! SimTrackerHit collection name
    std::string m_SimTrackerHitCollectionName;
    
    //! Digit collection name
    std::string m_digitCollectionName;	
         
    // Readout/DAQ parameters - set by user
    std::vector<int >  m_filterIDs;              //!< Digitize only sensors in this list
    float m_noiseFraction;                       //!< Fraction of noise hits
    float m_zsThreshold;                         //!< ZS threshold for zero suppression  
    int   m_frontEndType;                        //!< Front-end electronics type
    int   m_ADCRange;                            //!< ADC range is from 0 - ? (in electrons)
    int   m_ADCBits;                             //!< ADC has 0 - (2^m -1) digital values
    float m_ComparatorThr;                       //!< Comparator threshold (in e)
                                                   
    // Digitization parameters - set by user     
    double m_bulkDoping;                         //!< Net bulk doping concentration in sensors, in um^-3
    double m_topVoltage; 			             //!< Top plane voltage wrt. source, in V
    double m_backVoltage;                        //!< Back plane voltage wrt. source, in V
    double m_uSideBorderLength;                  //!< Border region with small drift fields, in um   
    double m_vSideBorderLength;                  //!< Border region with small drift fields, in um  
    double m_eGroupSize; 			             //!< Split Signalpoints in smaller groups of N electrons (in e)
    double m_eStepTime; 			             //!< Time step for tracking electron groups in readout plane (in ns)
    bool   m_electronicEffects;                  //!< Define if noise added or not
    double m_elNoise;                            //!< El. noise of individual pixels in ENC
    double m_maxSegmentLength;                   //!< Max length of path segment - space precision
    double m_tanLorentzAngle;                    //!< Tangent of Lorentz angle
    bool   m_integrationWindow;                  //!< Use integration window?
    double m_startIntegration;                   //!< Start time of integration of the sensors in ns (everything before this value will not be digitized)
    double m_stopIntegration;                    //!< Stop time of integration of the sensors in ns(everything after this value will not be digitized)
    
    // Magnetic field - obtained from Gear xml file and transformed into local system
    CLHEP::Hep3Vector m_magField;    //!< Magnetic field in T in detector reference system
    
    // Handle to detector data sheets 
    TBDetector m_detector;  
    
    // Current sensor number
    int m_ipl ;   
    int m_sensorID; 
    
   private:

    double m_timeCPU; //!< CPU time
    int    m_nRun ;   //!< Run number
    int    m_nEvt ;   //!< Event number 
   
    
  }; // Class

} // Namespace

#endif // SiPixDigitizer_H



