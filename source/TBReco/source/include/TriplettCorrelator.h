// /////////////////////////////////////////////////////////////////////////////////////// //
//                                                                                         //
//    TriplettCorrelator - Marlin Processor                                                //
// /////////////////////////////////////////////////////////////////////////////////////// //


#ifndef TriplettCorrelator_H
#define TriplettCorrelator_H 1

// Include DEPFETTrackTools header files
#include "TBDetector.h"

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <marlin/Exceptions.h>

// Include ROOT classes
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

// Include basic C
#include <string>
#include <cmath>
#include <map>
#include <vector>

namespace depfet {

//! TriplettCorrelator 
  /*! This processor provides a pre-alignment in XY for a tracking telescope in a 
   *  test beam setup. 
   *  
   *  The input to this processor is a sample of trackletts from an aligned part 
   *  of the tracking telescope. For example this can be the first arm of the
   *  reference telescope. The other input is a set of collections with hits
   *  from the other sensors. 
   * 
   *  The goal of the pre- alignment process is to correlate the xy of hits with 
   *  the extrapolations from the trackletts. Only the xy position of detectors 
   *  which do not have hits in the trackletts are corrected.
   *  
   *  The use of trackletts (tripletts) for pre-alignment offers two main advantages 
   *  over the standard correlator: 
   *
   *  1) We can cut on the chi2 value of trackletts to suppress combinatorial background. 
   *     This makes it easier to see real correlation bands. 
   *
   *  2) Trackletts offer a much more accurate track extrapolation from one end of 
   *     the telescope to the other.
   *   
   * 
   *  The two arguments apply especially for testbeams with multi GeV electrons when there 
   *  is a large distance between the arms of the reference telescope. 
   *
   *  Author: B.Schwenker, Universität Göttingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */

class TriplettCorrelator : public marlin::Processor {

 public:

//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new TriplettCorrelator ; }

//!Constructor - set processor description and register processor parameters
   TriplettCorrelator();
   
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
   
//! Histogram booking
/*! This method is used to book all the correlation histograms. It
 *  is called by processEvent when processing the first event.
 */
   void bookHistos();
    
protected:
   
//! Input hit collection names
/*! A vector containing all the collection names to be used in hit 
 *  correlations.
 */
   std::vector< std::string >  _inputHitCollectionNameVec; 
   
//! Track collection name
/*! These tracks provide a reference for the alignment. Typical these 
 *  are tripplet tracks from a telescope arm  
 */
   std::string _trackCollectionName;

//! Alignment DB file name
/*! Stores all alignment information (positions + orientations) of 
 *  detector setup. 
 */ 
   std::string _alignmentDBFileName;
   
//! Update alignment  
/*! Update the alignment data base using correlation band offset 
 *  values.  
 */
   bool _updateAlignment;

    
//! Output root file name  
/*!
 */
    std::string _rootFileName;       
           
 private:
    
//! ROOT_Output
    TFile * _rootFile;

//! Active flag for sensors 
   std::vector<bool> _isActive;
  
//! Correlation histograms
/*! 
 */
    std::map< unsigned int , TH2D* >  _hitUCorrelationMatrix;
    std::map< unsigned int , TH2D* >  _hitVCorrelationMatrix;
    
    std::map< unsigned int , TH1D* >  _hitUResidualHisto;
    std::map< unsigned int , TH1D* >  _hitVResidualHisto;
        
//! Run number
    int _iRun;
    
//! Event number
    int _iEvt;
                 
//! Handle to detector data 
    TBDetector _detector;    
     
}; // Class

} // Namespace

#endif 

