// /////////////////////////////////////////////////////////////////////////////////////// //
//                                                                                         //
//    Correlator - Marlin Processor                                                  //
// /////////////////////////////////////////////////////////////////////////////////////// //


#ifndef Correlator_H
#define Correlator_H 1

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

//! Correlator 
  /*! This processor provides a XY alignment for all sensors in a  
   *  test beam setup. In contrast to more precise track based alignment 
   *  methods, the correlator algorithm does not need tracks. This is 
   *  particularly useful, if the mis- alignemnt prohibits efficient track 
   *  finding. 
   *  
   *  The correlator algorithm utilizes two special attributes of beam test 
   *  experiments
   *  
   *  a) The position resolution even of a single hit is much better than 
   *     the mechanical alignment of the detectors.  
   * 
   *  b) The beam is highly collimated. The beam opening angle is small 
   *     (opening angle ~1° or smaller) and multiple scattering is frequently 
   *     limited due to high beam energies and/or small material budget.  
   *  
   *  A global coordinate system can be introduced by indentifying the z axis 
   *  with the beam axis. The misalinment of sensors in XY is the largest
   *  source of misalignment. Misalgnment in sensor tilts from the nominal values 
   *  will be neglected at this point. 
   *  
   *  The correlator algorithm creates seed track candidates from a single hit on 
   *  a user defined reference sensor. The seed track is extrapolated from sensor 
   *  to sensor using straight lines parallel to the z axis. The algorithms 
   *  now fills histograms with u and v residual seperately for all sensors. 
   * 
   *  After processing many events, the residuals histograms should show a peak 
   *  over a flat combinatorial background. Alignment corrections are computed by 
   *  centering the position of the residual peak to zero. 
   *
   *   
   *  Author: B.Schwenker, Universität Göttingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */

class Correlator : public marlin::Processor {

 public:

//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new Correlator ; }

//!Constructor - set processor description and register processor parameters
   Correlator();
   
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
   
//! Input hit collection names.
/*! A vector containing all the collection names to be used in hit 
 *  correlations.
 */
   std::vector< std::string >  _inputHitCollectionNameVec; 
   
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

//! New alignment  
/*! Don't use current alignment data base, but start from scratch   
 */
   bool _newAlignment;
    
//! Reference Plane
/*! Plane number of reference pixel module. Plane numbers are counted
 *  starting with zero, i.e. first module in beam.
 */
   int _refPlane;

//! Parameters for beam particles 
   double _momentum;
   double _mass;
   double _charge;
    
//! Output root file name  
/*!
 */
    std::string _rootFileName;       
           
 private:
    
//! ROOT_Output
    TFile * _rootFile;
  
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

