// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    X0ImageProducer - Marlin Processor                                                  //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef X0ImageProducer_H
#define X0ImageProducer_H 1

// DEPFETTrackTools includes
#include "TBDetector.h"


// Include basic C
#include <vector>
#include <string>
#include <map>

// Include LCIO classes
#include "lcio.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include ROOT classes
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>

namespace depfet {

//! X0ImageProducer Processor 
/*! 
 *  
 *  The processor allows a position resolved measurement of X/X0 using 
 *  reconstructed scatter kinks of telescope tracks.
 *
 *  The processor fills a root tree with measured scattering kinks at a 
 *  planar scattering target. The two projected scattering kinks are 
 *  measured along with their local intersection coordinates from a pair 
 *  of Kalman filters. A forward (backward) Kalman filter processes hits 
 *  in an upstream (downstream) arm of a high resolution reference 
 *  telescope. The root tree contains all information needed to compute 
 *  a 2D image of the radiation length X/X0 at the target plane. 
 *  
 *  The algorithms used are described in "Radiation length imaging with 
 *  high-resolution telescopes" by U. Stolzenberg, B. Schwenker et. al. 
 *  (Printed in the proceedings of VCI 2016). 
 * 
 *  Author: U.Stolzenberg, Universität Göttingen
 *  <mailto:ulf.stolzenberg@stud.uni-goettingen.de>
 *  
 *  Author: B.Schwenker, Universität Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
   
  
class X0ImageProducer : public marlin::Processor {
   
 public:
   
//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new X0ImageProducer ; }
    
//!Constructor - set processor description and register processor parameters
   X0ImageProducer();
   
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
   void bookHistos();
   
//!Method printing processor parameters
   void printProcessorParams() const;
   
// Processor Parameters
   
//! Input Track collection name
   std::string _downStreamTrackColName; 
   std::string _upStreamTrackColName; 
  
      
//! Alignment DB file name 
   std::string _alignmentDBFileName;
   
//! ROOT output file name  
   std::string _rootFileName;  
   
//! DUT plane number 
   int _idut; 

//! Max distance
   double _maxDist; 
     
     
// ROOT_OUTPUT 
   TFile * _rootFile;
   
   TTree * _rootMscTree;
   TTree * _rootEventTree;
    
   // Event tree variables 
   int _rootEventNumber;
   int _rootRunNumber;   
   int _rootnDownTracks;    
   int _rootnUpTracks;    
   int _rootNMatched;

   // Msc tree variables    
   double _rootTrackProbUp;
   double _rootTrackProbDown;
   double _rootTrackProbCombo;
 
   double _root_u; 
   double _root_v; 
   double _root_u_var; 
   double _root_v_var; 
   double _root_u_in; 
   double _root_v_in; 
   double _root_u_out; 
   double _root_v_out; 
   double _root_dudw;
   double _root_dvdw;
   double _root_angle1;
   double _root_angle2;
   double _root_angle1_var;
   double _root_angle2_var;
   double _root_momentum;

   double _root_vertex_mean_chi2;
   double _root_vertex_mean_prob;

   // Vertexing tree variables

   double _root_vertex_u;
   double _root_vertex_v;
   double _root_vertex_w;
   double _root_vertex_u_var;
   double _root_vertex_v_var;
   double _root_vertex_w_var;
   double _root_vertex_chi2;
   double _root_vertex_prob;
   double _root_vertex_u_res;
   double _root_vertex_v_res;
   int _root_vertex_multiplicity;
   int _root_vertex_id;
  
 private:
   
   // Handle to detector data 
   TBDetector  _detector;    
   
   // Few counter to show the final summary   

   
   
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   
}; // Class

} // Namespace

#endif 


