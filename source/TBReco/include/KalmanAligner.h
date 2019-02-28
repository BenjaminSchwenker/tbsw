// ////////////////////////////////////////////////////////////////////////// //
//                                                                            //
//    KalmanAligner - Marlin Processor                                        //
// ////////////////////////////////////////////////////////////////////////// //

#ifndef KalmanAligner_H
#define KalmanAligner_H 1


// TBTools includes   
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
 *  The processor has two inputs: 
 *  
 *  1) A root file containing the current nominal alignment DB   
 *  2) A lcio file containing a collection of tracks for detector alignment.  
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
      
//! Update alignmentDB   
   bool _updateAlignment;

//! New alignment  
/*! Don't use current alignment data base, but start from scratch   
 */
   bool _newAlignment;

//! Initial errors on alignment x shift [mm] for sensors ordered along beam line
   std::vector<float >  _errorsShiftX;
   
//! Initial errors on alignment y shift [mm] for sensors ordered along beam line
   std::vector<float >  _errorsShiftY;   
      
//! Initial errors on alignment z shift [mm] for sensors ordered along beam line
   std::vector<float >  _errorsShiftZ;

//! Initial errors on alignment alpha [rad] for sensors ordered along beam line
   std::vector<float >  _errorsAlpha;

//! Initial errors on alignment beta [rad] for sensors ordered along beam line
   std::vector<float >  _errorsBeta;

//! Initial errors on alignment gamma [rad] for sensors ordered along beam line
   std::vector<float >  _errorsGamma;

//! Number of tracks before the annealign is turned OFF
   int  _annealingTracks;

//! Scale factor for annealing schedule
   double  _annealingFactor;

//! LogLever during alignment
   int  _logLevel;

//! Maximum number of tracks passed to alignemnt
   int  _maxTracks;

//! Use beam model to constrain track fitting during alignment
   bool  _useBC;

//! P-Value cut for tracks used during alignment
   double  _pValueCut;

//! Reject alignment corrections exceeding DeviationCut*Sigma
   double  _deviationCut;

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

