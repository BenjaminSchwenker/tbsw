// ////////////////////////////////////////////////////////////////////////// //
//                                                                            //
//    BeamEnergyCorrector - Marlin Processor                                  //
// ////////////////////////////////////////////////////////////////////////// //

#ifndef BeamEnergyCorrector_H
#define BeamEnergyCorrector_H 1


// DEPFETTrackTools includes   
#include "TBDetector.h"

// Include ROOT classes
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include basic C
#include <string>


namespace depfet {

//! BeamEnergyCorrector processor 
/*! The task of this processor is to to apply a beam energy correction to 
 *  a collection of tracks. The processor creates a new collection of 
 *  corrected tracks that can be used instead of the orignal tracks.  
 * 
 *  A beam energy correction may be needed in test beams without a magnetic
 *  field where the track momentum and particle ID are given as external 
 *  parameters. Typically, parameter values are provided by the beam line
 *  instrumentation. This processor allows a fine tuning of the beam momentum
 *  by specifying a linear momentum field M(u,v)
 * 
 *  M(u,v) = M0 + Mu * u + Mv * v 
 *
 *  where u, v are local coordinates on a reference plane in the telescope. 
 *  The reference plane may be any plane specified in the gear file.  
 *   
 *  Author: Benjamin Schwenker, Göttingen University 
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */


class BeamEnergyCorrector : public marlin::Processor {
  
 public:
  
//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new BeamEnergyCorrector ; }
   
//!Constructor - set processor description and register processor parameters
   BeamEnergyCorrector();
   
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
   
//! Processor Parameters 
   
//! Input track collection name
   std::string _inputTrackCollectionName;

//! Output track collection name
   std::string _outputTrackCollectionName;
      
//! AlignmentDB file name 
   std::string _alignmentDBFileName;
   
//! Reference plane for applying correction
   int _iref;    

//! Momentum at centre of reference plane
   double  _M0;    

//! Momentum gradient along u direction
   double _Mu;    

//! Momentum gradient along v direction
   double _Mv;   
      
 private:
   
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   
   // Handle to detector data 
   TBDetector  _detector;        
   
   
   
}; // Class

} // Namespace

#endif 

