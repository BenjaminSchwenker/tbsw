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

//! track matching switch 
   bool _vertexfitswitch; 

//! toy study switch
   bool _m_toy;

//! toy study switch
   double _m_reco_error;

//! toy Bethe Heitler switch
   bool _m_ToyBetheHeitler;

//! Max vertex chi2/ndof
   double _maxVertexChi2;

//! Max distance
   double _maxDist;  
     
     
// ROOT_OUTPUT 
   TFile * _rootFile;
   
   TTree * _rootMscTree;
   TTree * _rootEventTree;
    
   // Event tree variables 
   // Each entry corresponds to a new event
   int _rootEventNumber;							// Event number of the track
   int _rootRunNumber;   							// Run number of the track
   int _rootnDownTracks;    						// Number of downstream tracks in this event
   int _rootnUpTracks;    							// Number of upstream tracks in this event
   int _rootNMatched;								// Number of matched upstream and downstream tracks in this event
													// If there are multiple downstream tracks that have been matched to a upstream track, it still only counts as 1 matched track 

   // Upstream track tree variables
   // Each entry corresponds to a upstream-downstream track combination that passed the chi2 criterium of the vertex fit  
   // This means that the root tree may contain multiple entries with the same upstream but different downstream tracks (see _rootVertexMultiplicity)

   double _rootTrackProbUp;							// p value of upstream track
   double _rootTrackProbDown;						// p value of downstream track
   double _rootTrackProbCombo;						// p value of combination of down and upstream track

   double _root_u; 									// u intersection coordinate on central plane, weighted mean of up- and downstream track intersection (mm)
   double _root_v; 									// v intersection coordinate on central plane, weighted mean of up- and downstream track intersection (mm)
   double _root_u_var; 								// u intersection variance on central plane, determined via weighted mean of up- and downstream track u intersection variance (mm^2)
   double _root_v_var; 								// v intersection variance on central plane, determined via weighted mean of up- and downstream track v intersection variance (mm^2)
   double _root_u_in; 								// u intersection coordinate on central plane, estimated from the upstream track (in-state) (mm)
   double _root_v_in;  								// v intersection coordinate on central plane, estimated from the upstream track (in-state) (mm)
   double _root_u_out; 								// u intersection coordinate on central plane, estimated from the downstream track (out-state) (mm)
   double _root_v_out; 								// v intersection coordinate on central plane, estimated from the downstream track (out-state) (mm)
   double _root_angle1;								// Projected scattering angle on target in u-w plane (rad)
   double _root_angle2;								// Projected scattering angle on target in v-w plane (rad)
   double _root_angle1_var;							// Variance of the projected scattering angle on target in u-w plane (rad^2)
   double _root_angle2_var;							// Variance of the projected scattering angle on target in u-w plane (rad^2)
   double _root_momentum;							// Nominal momentum of beam particle

   double _root_chi2;								// Chi2 value calculated from the up-downstream track residuals
   double _root_prob;								// p value value calculated from the up-downstream track residuals			

   int _rootVertexMultiplicity;						// Number of downstream tracks, which passed the vertexfit for this upstream track
   double _root_vertex_u;							// Vertex u position in the local coordinate system of the target plane (mm)
   double _root_vertex_v;							// Vertex v position in the local coordinate system of the target plane (mm)
   double _root_vertex_w;							// Vertex w position in the local coordinate system of the target plane (mm)
   double _root_vertex_u_var;						// Vertex u position variance in the local coordinate system of the target plane (mm^2)
   double _root_vertex_v_var;						// Vertex v position variance in the local coordinate system of the target plane (mm^2)
   double _root_vertex_w_var;						// Vertex w position variance in the local coordinate system of the target plane (mm^2)

   double _root_vertex_x;							// Vertex x position in the global coordinate system (mm)
   double _root_vertex_y;							// Vertex y position in the global coordinate system (mm)
   double _root_vertex_z;							// Vertex z position in the global coordinate system (mm)
   double _root_vertex_x_var;						// Vertex x position variance in the global coordinate system (mm^2)
   double _root_vertex_y_var;						// Vertex y position variance in the global coordinate system (mm^2)
   double _root_vertex_z_var;						// Vertex z position variance in the global coordinate system (mm^2)

   double _root_vertex_chi2ndf;						// chi2 value per degree of freedom of the vertex fit
   double _root_vertex_prob;						// p value of the vertex fit
   double _root_vertex_u_res;						// u residual of the vertex fit (mm)
   double _root_vertex_v_res;						// v residual of the vertex fit (mm)
   int _root_vertex_multiplicity;					// Number of downstream tracks, which passed the vertex fit for this upstream track
   int _root_vertex_id;								// Vertex ID, corresponds to the upstream track number 

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


