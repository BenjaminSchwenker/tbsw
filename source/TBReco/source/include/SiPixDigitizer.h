// //////////////////////////////////////////////////////////////////////////  //
//                                                                             //
//    SiPxlDigi - Marlin Processor - Digitizing data from DEPFET pixel sensors //
//                                                                             //
// //////////////////////////////////////////////////////////////////////////  //

/**
 *
 *
 * SiPxlDigi (a Marlin processor) represents a detailed
 * DEPFET detector digitizer. Additional information and validation of the simulation 
 * code is provided by the DEPFET Colaboration (http://www.depfet.org). 
 * 
 * THIS VERSION IS FOR THINNED DEPFETS ONLY!!!!!
 * 
 * First, a primary particle passing a DEPFET Pixel sensor produces a collection of  
 * SimTrackerHits. Each SimTrackerHits represents a Geant4 step. The number of 
 * SimTrackerHits is controlled via the Geant4 range cut. The range cut is set in the Mokka 
 * steering file. The default name for the input collection of SimTrackerHits is 
 * "PXDCollection". 
 * Second, all the digitization procedure is done within the sensor local coordinate system.
 * Each SimTrackerHit, more precisely its position and momentum, is transformed into local 
 * coordinates. If the track length of a SimTrackerHit exceeds a user defined maximum value, 
 * the SimTrackerHit is subdivided into smaller track segments called Ionisation points and 
 * the Geant4 energy loss in the SimTrackerHit is uniformly split between Ionisationpoints. 
 * Ionisationpoints are then drifted to the readout plane, shifted in the presence of magnetic
 * field (Lorentz shift) and smeared by a Gaussian distribution to account for diffusion during 
 * the drift time. Finally, the total charge of an Ionisation point is split into carrier groups 
 * of ~50 electrons. For each carrier group, we sample a random walk in the readout plane until 
 * a region close to the the internal gate of a pixel cell is reached. At this point, we assume 
 * that all charge of the carrier group is finally drifted into the internal gate of the 
 * corresponding pixel cell. 
 * Readout noise is modeled as a Gaussian signal component with mean zero. The rms value of readout 
 * noise is controlled via a marlin steering file. Readout noise is added to all pixels that collected
 * signal charge. Moreover, the sensor is populated with random noise hits. The number of noise hits is 
 * controlled by the signal threshold for zero supression. The collected charge can be stored either 
 * analog (in electrons) or digital in ADC units. Dynamic range and number of bits of ADC are configurable
 * by user. BelleII simulations should use 5bit ADC while testbeam simulations with CURO readout should use 
 * 14bit. 
 * The output of the Digitizer is a collection of zero suppressed pixels. Optionally, full analog frames
 * can be created for testing/debugging. Optionally, the digitizer can create mc truth hits for sensor
 * resolution studies.   
 *
 * If you have any questions or comments, please write an email to the authors of this code:
 * <a class="email" href="mailto:benjamin.schwenker@phys.uni-goettingen.de"> Benjamin Schwenker</a>,
 * Georg August University, Goettingen.
 * <a class="email" href="mailto:drasal@ipnp.troja.mff.cuni.cz"> Zbynek Drasal</a>,
 * Charles University, Prague.
 *
 *\namespace CLHEP        Namespace of Class Library in High Energy Physics
 *\namespace lcio         Namespace of Linear Collider InputOutput
 *\namespace marlin       Namespace of Modular Analysis and Reconstruction tool for LINear collider
 *\namespace std          Standard library namespace
 */

#ifndef SIPXLDIGI_H
#define SIPXLDIGI_H 1

// Define ROOT output if needed
#define ROOT_OUTPUT

// Include basic C
#include <vector>

// Include CLHEP classes
#include <CLHEP/Vector/ThreeVector.h>

// Include Digi header files
#include "SiPxlGeom.h"


// Include LCIO classes
#include <lcio.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>

// Include ROOT classes
#ifdef ROOT_OUTPUT
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#endif

#include <TRandom.h>
#include <TRandom3.h>
#include <TMath.h>

namespace sipxl {

// Define constants
#define MAXLAYERS 20
#define MAXPADS 2000

// Structures
struct Digit
{
   int cellIDZ;
   int cellIDRPhi;
   double cellPosZ;
   double cellPosRPhi;
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

//! Marlin processor intended for detailed digitization of DEPFET pixel detectors
//!
//! @author of original algorithm A. Raspereza, MPI Munich (ILCsoft VTXDigitizer)
//! @author of the code and geometry Z. Drasal, Charles University, Prague
//! @author of the new algorithm + validation B. Schwenker, University of Göttingen

class SiPxlDigi : public marlin::Processor {

 public:

//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new SiPxlDigi ; }

//!Constructor - set processor description and register processor parameters
   SiPxlDigi();

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
   void ProduceDigits(const SignalPointVec & signalPoints, DigitVec & digits);

//!Method for calculating electronic effects on signal digits
   void ProduceNoiseEffects(DigitsMap & digitsMap);

//!Method producing noise digits and updating PXD map containing all digits
   void ProduceNoiseDigits(DigitsMap & digitsMap);

//!Method producing sparsified pixels output from all digits (noise + signal)
   void ProduceSparsePixels(LCCollectionVec * pixels , const DigitsMap & digitsMap);

//!Method producing full analog output matrix output from all digits (noise + signal)
   void ProduceFullMatrix(LCCollectionVec * matrix , const DigitsMap & digitsMap);

//!Method producing local truth hits from geant4 steps 
   void ProduceTruthHits(LCCollectionVec * truthhit, LCCollection * simTrkHits);
   


// OTHER METHODS

//!Method transforming hit to local coordinate system using SiPxlGeom interface
   void TransformToLocal(const SimTrackerHit * simTrkHit, SpacePoint & hitLocal);

//!Method transforming hit to global coordinates using SiPxlGeom interface
   void TransformToGlobal(const SpacePoint &  hitLocal, SpacePoint & hitGlobal);

//!Method updating PXD map (containing all digits) using given vector of digits
   void UpdateDigitsMap(DigitsMap & digitsMap, DigitVec & digits);

//!Method returning collected charge in ADC units
   int getInADCUnits(double charge);

//!Method providing mathematical round
   int round(double num);


// PRINT METHODS

//!Method printing processor parameters
   void printProcessorParams() const;


// VARIABLES

// Collection names
   std::string _inColName;       		//!< LCIO input collection name (SimTrackerHits)
   std::string _sparseDataCollectionName;	//!< LCIO output collection name (Sparsified Pixels) 

// Readout/DAQ parameters - set by user
   float _Fraction;                       //!< Fraction of noise hits
   float _ZSthr;                          //!< ZS threshold for zero suppression  
   bool  _ADC;                            //!< Simulate ADC?
   int   _ADCRange;                       //!< ADC range is from 0 - ? (in electrons)
   int   _ADCBits;                        //!< ADC has 0 - (2^? -1) digital values
   
   bool   _createFullMatrix;              //!< Define if full analog output is produced
   bool   _createTruthHits;               //!< Define if truth hits is produced


// Digitization parameters - set by user
   bool   _bricked;                       //!< Define if bricked structure of pixels used
   bool   _doublePixel;                   //!< Define if double pixel structure is used
   double _bulkDoping;                    //!< Net bulk doping concentration in sensors, in um^-3
   double _Utop; 			  //!< Top plane voltage wrt. source, in V
   double _Uback;                         //!< Back plane voltage wrt. source, in V
   FloatVec _sourceBorderLength;          //!< Border region with small drift fields, in um   
   FloatVec _drainBorderLength;           //!< Border region with small drift fields, in um  
   FloatVec _clearBorderLength;           //!< Border region with small drift fields, in um  
   double _eGroupSize; 			  //!< Split Signalpoints in smaller groups of N electrons (in e)
   double _eStepTime; 			  //!< Time step for tracking electron groups in readout plane (in ns)
   bool   _electronicEffects;             //!< Define if noise added or not
   FloatVec _elNoise;                     //!< El. noise of individual pixels in ENC
   bool   _PoissonSmearing;               //!< Define if Poisson smearing set or unset (only if electronics effects defined)
   double _maxSegmentLength;              //!< Max length of path segment - space precision
   double _tanLorentzAngle;               //!< Tangent of Lorentz angle
   bool   _integrationWindow;             //!< Use integration window?
   double _startIntegration;              //!< Start time of integration of the sensors in ns (everything before this value will not be digitized)
   double _stopIntegration;               //!< Stop time of integration of the sensors in ns(everything after this value will not be digitized)
    
   
   
// Digitization parameters - obtained from Gear xml file
   short int _currentLayerID;      //!< Actual layer ID
   short int _currentLadderID;     //!< Actual ladder ID
   short int _currentSensorID;     //!< Actual sensor ID
   float     _sensorThick;         //!< Actual sensor - Si wafer thickness in system of units
   float     _sensorWidth;         //!< Actual sensor width in system of units
   float     _sensorLength;        //!< Actual sesnor length in system of units

// Magnetic field - obtained from Gear xml file and transformed into local system
   CLHEP::Hep3Vector _magField;    //!< Magnetic field in T in detector reference system

// Geometry parameters
   SiPxlGeom * _geometry;          //!< All geometry information from Gear xml file

// Random generator  
   TRandom3 * myRng;


#ifdef ROOT_OUTPUT

   TFile * _rootFile;
   TTree * _rootTreeHits;
   TTree * _rootTreeSimHits;
      
   int _rootLayerID;
   int _rootLadderID;
   int _rootSensorID;

   // branches in Hit tree
   int _rootPrimaryPDG;         // PDG code of primary particle           
   double _rootX;               // truth hit x [px]
   double _rootY;               // truth hit y [px] 
   double _rootZ;                
   double _rootTOT;             // eLoss  [keV]  
  
  
   // branches in _rootTreeSimHits
   int    _rootPDG;              // PDG code of geant4 step   
   double _rootTrackLength;      // length of geant4 step [um]
   double _rootMomentum;         // momentum of geant4 step  [in keV] 
   double _rootdEdx;             // eLoss along geant4 step [in keV]
   double _rootXPos; 
   double _rootYPos;
   double _rootZPos;


   TH1D * _rootChargeCollectionTime; // charge collection time [in ns]
   TH1D * _rootRandomWalkSteps;      // number of random walk steps 
   TH1D * _rootKilledRandomWalks;    // number of killed random walks  

#endif

 private:

   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number 
   

}; // Class

} // Namespace

#endif // SIPXLDIGI_H



