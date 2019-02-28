// /////////////////////////////////////////// //
//                                             //
//    BeamSpotAligner - Marlin Processor      //
// /////////////////////////////////////////// //


#ifndef BeamSpotAligner_H
#define BeamSpotAligner_H 1

// Include DEPFETTrackTools header files
#include "TBDetector.h"


// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <marlin/Exceptions.h>

// Include ROOT classes
#include <TFile.h>
#include <TH2D.h>

// Include basic C
#include <string>
#include <map>


namespace depfet {

//! BeamSpotAligner 
/*! This processor estimates the center of the beam spot in local sensor 
 *  coordinates from a hit collection. Aftwards, the sensor is moved in 
 *  global coordinates to make the global z axis hit the center of the 
 *  beam spot. The idea of this alignment step is to let the z axis 
 *  coincide with the beam axis. 
 *    
 *  One of the problem with telescope alignment is the weak shearing mode in the 
 *  X and Y direction. The information that all tracks originate from a 
 *  highly collimated particle beam allows to set a constraint on shearing
 *  modes. On possibility is to estimate the position of the beam axis 
 *  relative to the center of two sensors in the telescope.  
 *  
 *  The method is not applicable in all cases. It is needed that the beam 
 *  axis is really crossing at least one sensor. It is also needed to 
 *  have an approximately Gaussian hit density where the position of the 
 *  beam axis can be estimated as a maximum of the hit map. This can be 
 *  obscured when the trigger is computed from the coincidence of scintilators 
 *  at both ends of the telescope. 
 * 
 *  The proposed alignment schema is to align the telescope from an air run 
 *  where nothing is in between the telesocpe arms. The beam collimators should 
 *  be tuned to give a narrow spot. Only one big scintillator behind the last 
 *  plane of the telescope should be used for triggering.    
 *   
 *  Author: B.Schwenker, Universität Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class BeamSpotAligner : public marlin::Processor {

 public:

//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new BeamSpotAligner ; }

//!Constructor - set processor description and register processor parameters
   BeamSpotAligner();
   
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

//! Histogram booking
/*! This method is used to book all the correlation histograms. It
 *  is called by processEvent when processing the first event.
 */
   void bookHistos();
   
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

//! New alignment  
/*! Don't use current alignment data base, but start from scratch   
 */
   bool _newAlignment;
        
//! Output root file name  
/*!
 */
    std::string _rootFileName;       
           
 private:
    
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
  
    // Handle to detector data 
    TBDetector  _detector;    
    
    // Handle to root file
    TFile * _rootFile;
    std::map< std::string, TH2D *> _histoMap2D;  
    
    
}; // Class

} // Namespace

namespace beamspotfitter { 

 /** Bivariate gaussian fit model in sensor u,v coordinates 
   *  Summary of model parameters: 
   *
   * par[0] :	amplitude factor 
   * par[1] : 	mean in u in mm
   * par[2] :  	standard deviation in u in mm
   * par[3] : 	mean in v in mm
   * par[4] :   standard deviation in v in mm
   * par[5] : 	correlation coefficient
   *
   * Summary of variables 
   *
   * x[0]   :   u coordinate in mm
   * x[1]   :	v coordinate in mm
   */ 
  double gauss2D(double *x, double *par);   

   /** Bivariate gaussian fit model in sensor u,v coordinates 
   *  Summary of model parameters: 
   *
   * par[0] :	amplitude factor 
   * par[1] : 	mean in u in mm
   * par[2] :  	standard deviation in u in mm
   * par[3] : 	mean in v in mm
   * par[4] :   standard deviation in v in mm
   *
   * Summary of variables 
   *
   * x[0]   :   u coordinate in mm
   * x[1]   :	v coordinate in mm
   */ 
  double gauss2Dsimple(double *x, double *par);  
}


#endif 

