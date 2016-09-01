// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    X0ImageSimulator - Marlin Processor - Validation of material budget estimation            //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef X0ImageSimulator_H
#define X0ImageSimulator_H

// Include DEPFETTrackTools header files
#include "TBDetector.h"
#include "TBTrack.h"

// Include ROOT classes
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>

// marlin includes ".h"
#include "marlin/DataSourceProcessor.h"

// lcio includes <.h>

// system includes <>
#include <vector>

namespace depfet
{
  
  /** The X0ImageSimulator Module
   * The processor provides a simple test environment for 
   * the validation material estimation. A gear file
   * is used to specify the geometrical configuration of the 
   * test beam setup. The processor propagates a number of 
   * particles through the setup and creates track hits on 
   * the detectors. Afterwards, tracks are built and 
   * a track fit is carried out. The user is supplied a root
   * output file containing histograms for track chisq,
   * p-values, pulls and so on.
   *  
   * Author: Benjamin Schwenker, Göttingen University 
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de> 
   */
  
  class X0ImageSimulator:public marlin::DataSourceProcessor 
  {
   public:
     
    //! Default constructor
    X0ImageSimulator();
    virtual X0ImageSimulator * newProcessor ();
    virtual void readDataSource (int Ntrig);
    virtual void init ();
    virtual void end ();
    
    // Compute tracking pulls, difference to truth dividided by tracking errors
    void fill_tracking_pulls(TBTrack & reco, TBTrack &  truth);
    

    //! Histogram booking
    void bookHistos();
    
   protected:
    
    // Processor parameters 
    
    //! ROOT file name Output 
    std::string _rootFileName;  

    //! Alignment DB file name 
    std::string _alignmentDBFileName;
    
    //! Number of simulated tracks
    int _maxTracks; 
 
    //! position of the DUT in the beamline
    int _idut; 

    //! Choose model of multiple scattering for the telescope ( Highland:0 or Frühwirth:1 )
    int _mscmodel;

    // Choose digitizer model for pixel sensors ( Mimosa26:0 or Gaussian:1 or Box:2)
    int _digi_type;
    
    //! Hit detection efficiency 
    double _hitEfficiency; 
        
    //! Particle gun
    double _GunXPosition;
    double _GunYPosition;
    double _GunZPosition; 
    double _GunRotX;
    double _GunRotY;
    double _GunRotZ;
    double _GunSpotSizeX;
    double _GunSpotSizeY;
    double _GunDivergenceX;
    double _GunDivergenceY;
    double _GunCorrelationX; 
    double _GunCorrelationY;      

    double _BetheHeitlerT0;
    bool _useTrueEnergy; 
     
    double _momentum;
    double _momentum_error;
    double _mass;
    double _charge;

   private: 
             
    // Handle to detector data sheets 
    TBDetector _detector;   
    
    TFile * _rootFile;
     
    std::map< std::string, TH1D* > _histoMap;
   
    TTree * _rootMscTree;

    int _rootEventNumber;   
    int _rootRunNumber;      
    int _rootDaqID; 
    int _rootPlaneID;
    int _rootTrackHits;
    double _rootTrackChi2;
    double _rootTrackProb;
    double _rootTrackProbUp;
    double _rootTrackProbDown;
    double _rootTrackProbCombo;
    double _root_x; 
    double _root_y; 
    double _root_dxdz;
    double _root_dydz;
    double _root_u; 
    double _root_v; 
    double _root_u_in;    
    double _root_v_in;     
    double _root_u_out;    
    double _root_v_out;    
    double _root_dudw;
    double _root_dvdw;
    double _root_angle1;
    double _root_angle2;
    double _root_angle1_err;
    double _root_angle2_err;
    double _root_momentum;
    double _root_truth_momentum;
    double _root_angle1_truth;    
    double _root_angle2_truth;    
     
  };
  
  X0ImageSimulator gX0ImageSimulator;
 
}            
#endif
