// ////////////////////////////////////////////////////////////////////////// //
//                                                                            //
//    KalmanAligner - Marlin Processor                                        //
// ////////////////////////////////////////////////////////////////////////// //

#ifndef KalmanAligner_H
#define KalmanAligner_H 1


// DEPFETTrackTools includes   
#include "TBDetector.h"
#include "AlignEvent.h"

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
#include <map>

namespace depfet {

//! KalmanAligner processor 
/*! The task of this processor is carry out a track based alignment of 
 *  a tracking telescope. The KalmanAligner processor estimates the deviations
 *  between the current detector positions and the real detector positions 
 *  inside the global telescope coordinate system. The current detector 
 *  positions are obtained from merging the nominal detector positions 
 *  from the gear file and all known corrections from an alignmentDB file. 
 *  After the alignment is finished, the alignmentDB file is overwritten 
 *  with new corrections.  
 *  
 *  The Kalman Alignment Algorithm (KAA) method is iterative and based on the 
 *  Kalman filter equations for minimization of the track residuals chisqu as 
 *  a function of the detector alignment constants and the track parameters. 
 *  
 *  In particular, the KAA offers the possibility to use any prior information
 *  about the detector alignment from mechanical measurements and/or previous 
 *  alignment procedures. It easily allows to fix certain alignment parameters
 *  to define a global reference frame. The method is also suitable for alignment 
 *  of individual DUT sensors relative to the tracking telescope. 
 *  
 *  The processor has three major inputs: 
 *  
 *  1) A root file containing the current nominal alignment DB  
 *  2) A text based config file for setting initial errors and outlier rejection 
 *  3) A lcio file containing a collection of tracks for detector alignment.  
 *   
 *  The processor fills the data structures for the alignment algorithm, runs the 
 *  alignment code and updates the alignment DB lcio file. A root file is written 
 *  containing data driven validation plots for the new alignment constants. 
 *   
 *  The Kalman Alignment Algorithm goes back to R. Fruehwirth and E. Widl. A 
 *  detailed explaination is given in the Phd thesis "Global Alignment of the 
 *  CMS Tracker" by E. Widl. 
 *   
 *  Author: Benjamin Schwenker, Göttingen University 
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */


class KalmanAligner : public marlin::Processor {
  
 public:
  
//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new KalmanAligner ; }
   
//!Constructor - set processor description and register processor parameters
   KalmanAligner();
   
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
   
//! Track collection name
   std::string _trackCollectionName;
      
//! AlignmentDB file name 
   std::string _alignmentDBFileName;
   
//! AlignConfig file name 
   std::string _alignConfigFileName;
      
//! Update alignmentDB   
   bool _updateAlignment;

//! New alignment  
/*! Don't use current alignment data base, but start from scratch   
 */
   bool _newAlignment;
      
      
 private:
   
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   
   // Count reco tracks used for alignment 
   int _nKAATracks;
   
   // Handle to detector data 
   TBDetector  _detector;        
   
   //! Alignment validation data containers 
   // The file contains trees and histos to 
   // monitor alignment quality and progress.
    
   TFile * alignment_data;
      
   //! Alignment data containers 
   // Temporary buffers for alignment data 
   // needed by KAA.
   
   AlignEvent * myEvent;  
   TTree * AlignTree;  
   
}; // Class

} // Namespace

#endif 

