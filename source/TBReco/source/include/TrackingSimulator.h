// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    TrackingSimulator - Marlin Processor                                                  //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef TrackingSimulator_H
#define TrackingSimulator_H

// Include DEPFETTrackTools header files
#include "TBDetector.h"
#include "TBTrack.h"

// Include ROOT classes
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

// marlin includes ".h"
#include "marlin/DataSourceProcessor.h"

// lcio includes <.h>

// system includes <>
#include <vector>

namespace depfet
{
  
  /** The TrackingSimulator Module
   * The processor provides a simple test environment for 
   * the validation of track fitting algorithms. A gear file
   * is used to specify the geometrical configuration of the 
   * test beam setup. The processor propagates a number of 
   * particles through the setup and creates track hits on 
   * the detectors. Afterwards, tracks are built and 
   * a track fit is carried out. The user is supplied a root
   * output file containing histograms for track chisq,
   * p-values, pulls and so on.
   *  
   * A Marlin steering file is used to control the simulation.
   * There exist parameters for number of tracks and spectra
   * of particle vertex positions, slopes and energies. 
   * Material effects like multiple scattering may be turned 
   * on/off. The user can define the method for smearing truth
   * hits with detector resultions (BOX/GAUSSIAN).  
   * 
   * The gear file is used to built a simplified geometry. All 
   * tracking detectors are treated as thin planes with given 
   * Gaussian resolutions and radiation length. All Passive 
   * materials are traeted as a thin scatterers with Gaussian
   * scattering.
   * 
   * The simplified simulation means that the assumptions for 
   * Kalman tracking are exactly true: 
   *  	
   * - known/variable contamination with misses and outliers   
   * - scatterers are thin and exactly known
   * - measurment errors are gaussion and known 
   * - same materials and scaterring formula in sim and fit
   * - no misalignment of the detector
   * 
   * In other words: there is no(!) excuse for bad pulls :) 
   * 
   * Author: Benjamin Schwenker, Göttingen University 
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de> 
   */
  
  class TrackingSimulator:public marlin::DataSourceProcessor 
  {
   public:
     
    //! Default constructor
    TrackingSimulator();
    virtual TrackingSimulator * newProcessor ();
    virtual void readDataSource (int Ntrig);
    virtual void init ();
    virtual void end ();
    
    // Compute tracking pulls, difference to truth dividided by tracking errors
    void fill_tracking_pulls(TBTrack & reco, TBTrack &  truth);
    
    // Make a histogram of track chi2, chi2/ndof and chi2-Probability. 
    void fill_track_chi2(TBTrack & reco);
    
    //! Histogram booking
    void bookHistos();
    
   protected:
    
    // Processor parameters 
    
    //! ROOT file name Output 
    std::string _rootFileName;  
    
    //! Number of simulated tracks
    int _maxTracks;  
    
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
    
    double _momentum;
    double _mass;
    double _charge;

    double _BetheHeitlerT0;
    bool _useTrueEnergy; 

    double _beamIntensity; 
    int _trgmode; 
    int _numit; 
      
    int _minHits;

   private: 
             
    // Handle to detector data sheets 
    TBDetector _detector;      
    
    TFile * _rootFile;
     
    std::map< std::string, TH1D* > _histoMap;
    std::map< std::string, TH2D* > _histoMap2D;
    
  };
  
  TrackingSimulator gTrackingSimulator;
 
}            
#endif
