// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    M26DigitizerSimulator - Marlin Processor                                              //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef M26DigitizerSimulator_H
#define M26DigitizerSimulator_H

// Include DEPFETTrackTools header files
#include "TBDetector.h"

// Include ROOT classes
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

// marlin includes ".h"
#include "marlin/DataSourceProcessor.h"

// lcio includes <.h>

// system includes <>
#include <vector>

namespace depfet
{
  
  /** The M26DigitizerSimulator Module
   * The processor provides a simple test environment for 
   * non Gaussian digitizers. The most important example 
   * is the Mimosa26Digitizer
   * 
   * Author: Benjamin Schwenker, Göttingen University 
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de> 
   */
  
  class M26DigitizerSimulator:public marlin::DataSourceProcessor 
  {
   public:
     
    //! Default constructor
    M26DigitizerSimulator();
    virtual M26DigitizerSimulator * newProcessor ();
    virtual void readDataSource (int Ntrig);
    virtual void init ();
    virtual void end ();
    
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
  
   private: 
             
    // Handle to detector data sheets 
    TBDetector _detector;      
    
    TFile * _rootFile; 
    TTree * _rootHitTree;
    
    int _rootDetectorID;
    int _rootHitQuality;       // GoodCluster == 0
    double _rootHitU;         // in mm        
    double _rootHitV;         // in mm  
    double _rootHitCharge; 
    double _rootHitSeedCharge;
    int _rootHitSize;  
    int _rootHitSizeCol;     
    int _rootHitSizeRow;    
    int _rootHitCol;
    int _rootHitRow; 
    int _rootHitHasTrack;             // Matched to track != 0         
    double _rootHitFitU;              // Track impact point [mm] - in local frame       
    double _rootHitFitV;  
    double _rootHitTruthU; 
    double _rootHitTruthV;            
    double _rootHitFitdUdW;           // Track direction tangent [rad] - in local frame       
    double _rootHitFitdVdW;      
    double _rootHitFitErrorU;          // rms error u
    double _rootHitFitErrorV;          // rms error v 
    double _rootHitPullResidualU; 
    double _rootHitPullResidualV;             
    int _rootHitFitCol;               // DUT pixel - readout column         
    int _rootHitFitRow;               // DUT pixel - readout row   
    double _rootHitFitPixU;        // Hit pixel, center in u[mm]        
    double _rootHitFitPixV;        // Hit pixel, center in v[mm] 
    double _rootHitTrackChi2 ; 
    int _rootHitTrackNDF; 
    double _rootHitLocalChi2; 
    
  };
  
  M26DigitizerSimulator gM26DigitizerSimulator;
 
}            
#endif
