// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    AlignmentSimulator - Marlin Processor                                                  //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef AlignmentSimulator_H
#define AlignmentSimulator_H


// Include DEPFETTrackTools header files
#include "TBDetector.h"
#include "TBTrack.h"
#include "AlignEvent.h"

// Include ROOT classes
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom.h>
#include <TH1.h>
#include <TH2.h>

// marlin includes ".h"
#include "marlin/DataSourceProcessor.h"

// lcio includes <.h>

// system includes <>
#include <vector>

namespace depfet
{
  
  /** The AlignmentSimulator Module
   * The processor provides a simple test environment for the 
   * validation of alignment algorithms. Here, the focus in on 
   * track based alignment of high resolution pixel tracking 
   * telescopes. A typical example is the EUDET/AIDA telescope. 
   *  
   * A gear file is used to specify the geometrical configuration of 
   * the detector. It is precisely the same data model as used for 
   * the reconstruction of test beam data from the EUDET telescope. 
   * The model features a full rigid body alignment with 3 shifts 
   * and 3 rotations per sub- detector.
   * 
   * A Marlin steering file is used to control the simulation.
   * There exist parameters for number of tracks and spectra
   * of particle vertex positions, slopes and energies. 
   * Material effects like multiple scattering may be turned 
   * on/off. The user can define the method for smearing truth
   * hits with detector resultions (BOX/GAUSSIAN).   
   *  
   * The processor reads an initial gear file with a nominal detector
   * geometry. Here, nominal means the initial guess for the true 
   * detector positions with misalignment. Then, a misaligned (true) 
   * detector geometry is created by randomly displacing detectors from
   * their nominal positions. A sample of particle tracks is simulated in 
   * misaligned(!) detector and fitted in the nominal(!) geometry.   
   * Then, this track sample is used for running alignment code. Root 
   * file is created showing the quality of the alignment constants. 
   *  
   * THERE IS NO EXCUSE FOR BAD PULLS :)
   * 
   * Author: Benjamin Schwenker, Göttingen University 
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de> 
   */
     
  class AlignmentSimulator:public marlin::DataSourceProcessor 
  {
   public:
     
    //! Default constructor
    AlignmentSimulator();
    virtual AlignmentSimulator * newProcessor ();
    virtual void readDataSource (int Ntrig);
    virtual void init ();
    virtual void end ();
    
    // Compute align pulls, difference to truth dividided by align error  
    void fill_align_pulls(CLHEP::HepVector& reco_state, CLHEP::HepSymMatrix& reco_cov, 
                          CLHEP::HepVector& truth_state );
    
    // Compute tracking parameter resolutions
    void fill_tracking_histos(TBTrack & reco);
     
    //! Histogram booking
    void bookHistos();
    
   protected:
    
    // Processor parameters 
    
    //! ROOT file name 
    std::string _rootFileName;  
    
    //! AlignConfig file name 
    std::string _alignConfigFileName;
    
    //! Generate n misaligned detectors 
    int _nMisalign;
    
    //! Generate m tracks per detector 
    int _nTracks; 
    
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
    
    
    //! Particle ID 
    double _momentum;
    double _mass;
    double _charge;
    
    double _beamIntensity; 
    int _trgmode; 
    
    //! Skip certain planes from misalignment 
    std::vector<int >  _DontMisAlignPlanes; 
    
    
    
    //! Shift of sub-detector origin (XYZ) 
    double _PositionAlignmentErrorX;    
    double _PositionAlignmentErrorY;
    double _PositionAlignmentErrorZ;
             
    //! Rotation arounf sub-detector axes (UVW) 
    double _RotationAlignmentErrorU; 
    double _RotationAlignmentErrorV;
    double _RotationAlignmentErrorW;  
       
   private: 
             
    // Handle to detector data sheets 
    TBDetector _detector;      
    
    // Flags reference sensor planes 
    std::vector<bool> _isMisAligned; 
    std::vector<bool> _isAligned; 
    
    TFile * _rootFile;
    
    // Buffer tracks in alignment 
    TTree * myAlignTree; 
    AlignEvent * myEvent; 

    // Buffer hits in pre-alignment
    TTree* myHitTree; 
    double myHitU;           
    double myHitV;
    int myPlane;     
    
    // Store for final histograms 
    std::map< std::string, TH1D* > _histoMap;
    std::map< std::string, TH2D* > _histoMap2D;
    
  };
  
  AlignmentSimulator gAlignmentSimulator;
 
}            
#endif
