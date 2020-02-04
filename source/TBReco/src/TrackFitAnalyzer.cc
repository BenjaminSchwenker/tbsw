// TrackFitAnalyzer implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// Local includes 
#include "TrackFitAnalyzer.h"

// TBTools includes
#include "TBDetector.h"
#include "TBTrack.h"
#include "TrackInputProvider.h"
#include "GenericTrackFitter.h"
#include "PixelCluster.h"
#include "ThreeDModel.h"

// C++ includes
#include <iomanip>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <Exceptions.h>

#include <TMath.h>

// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;
using namespace std::string_literals;
namespace depfet {

//
// Instantiate this object
//
TrackFitAnalyzer aTrackFitAnalyzer ;

//
// Constructor
//
TrackFitAnalyzer::TrackFitAnalyzer() : Processor("TrackFitAnalyzer")
{
   
// Processor description
  _description = "TrackFitAnalyzer: Create tree for monitoring the track fit quality";
   

  //
  // Input collections  
  registerInputCollection(LCIO::TRACK,"InputTrackCollectionName",
                          "Track input collection",
                          _inputTrackCollectionName,std::string("tracks"));
   
  // 
  // Processor parameters
  registerProcessorParameter( "RootFileName",
                              "Output root file name",
                              _rootFileName, 
                              std::string("Histos.root"));

  registerProcessorParameter ("ReferencePlane",
                              "Plane number of teference plane. Use only tracks having a hit on the reference plane to measure DUT efficiency. Put -1 to deactivate.",
                              _iref,  static_cast < int > (-1));
                                 
}

//
// Method called at the beginning of data processing
//
void TrackFitAnalyzer::init() {
  
  // Initialize variables
  _nRun = 0 ;
  _nEvt = 0 ;
   
  if (_iref < 0 && _iref >= TBDetector::GetInstance().GetNSensors()) {
    streamlog_out ( MESSAGE3 )  << "Steering parameter reference plane has invalid value and will be ignored."
                                << endl << endl;  
    _iref = -1;
  }
  
  // Print set parameters
  printProcessorParams();
  
  // CPU time start
  _timeCPU = clock()/1000;
  
  // Book all needed histograms 
  bookHistos();
}

//
// Method called for each run
//
void TrackFitAnalyzer::processRunHeader(LCRunHeader * run)
{

// Print run number
   streamlog_out(MESSAGE3) << "Processing run: "
                           << (run->getRunNumber())
                           << std::endl << std::endl;

   _nRun++ ;

}


//
// Method called for each event
//
void TrackFitAnalyzer::processEvent(LCEvent * evt)
{
    
  //////////////////////////////////////////////////////////////////////  
  // Process next event
  ++_nEvt;
   
  if ( _nEvt % 1000 == 0 ) {
    streamlog_out( MESSAGE3 ) << "Processing event "
                              << evt->getEventNumber() << " in run "
                              << evt->getRunNumber() << endl; 
                               
  }

  TrackInputProvider TrackIO; 
  
  GenericTrackFitter TrackFitter(TBDetector::GetInstance());
  TrackFitter.SetNumIterations(1); 
   
  LCCollection* inputCollection;
  try {
      inputCollection = evt->getCollection(_inputTrackCollectionName);
  } catch (DataNotAvailableException& e) {
      throw SkipEventException(this);
  }
  
  // Main loop over all tracks
  int nTracks = inputCollection->getNumberOfElements(); 
  
  
  for (int itrk = 0; itrk < nTracks; itrk++) {
    
    // Retrieve track from LCIO 
    Track * inputtrack = dynamic_cast<Track*> (inputCollection->getElementAt(itrk));
    
    // Convert LCIO -> TB track  
    TBTrack track = TrackIO.MakeTBTrack( inputtrack, TBDetector::GetInstance() );  
    
    // ReFit track 
    bool trkerr = TrackFitter.Fit(track);
    if ( trkerr ) {
      continue;
    }  
    
    // Fill tree  
    // ===============
    
    _rootEventNumber = evt->getEventNumber();           
    _rootRunNumber = evt->getRunNumber();                  
    _rootTrackChi2 = track.GetChiSqu();          
    _rootTrackNDF = track.GetNDF();                
    _rootTrackNHits = track.GetNumHits();            

    _rootTrackWithRefHit = -1;
    // Check track has a hit on reference (timing) plane
    if ( _iref >= 0 ) {
      if ( track.GetTE(_iref).HasHit() ) {
        streamlog_out ( MESSAGE2 ) << "Track has hit on reference plane." << endl;
        _rootTrackWithRefHit = 0;
      }  
    } 

    //
    // Sensor level histograms 
  
    // Get number of sensors
    int nSens = TBDetector::GetInstance().GetNSensors(); 


    for (int ipl= 0; ipl< nSens; ++ipl) {  
        
      // Get sensor data 
      //------------------------
      const TBTrackElement& TE = track.GetTE(ipl);  
      const Det& dut = track.GetTE(ipl).GetDet();
           
      // Get local track parameters 
      double trk_tu = TE.GetState().GetPars()[0];  // rad
      double trk_tv = TE.GetState().GetPars()[1];  // rad
      double trk_u = TE.GetState().GetPars()[2];   // mm
      double trk_v = TE.GetState().GetPars()[3];   // mm
      double trk_qp = TE.GetState().GetPars()[4];   // 1/GeV
         
      double trk_charge = track.GetCharge();
      double trk_mom = std::abs(trk_charge/trk_qp); 

      // Get readout channels  
      int fitcellu = dut.GetUCellFromCoord( trk_u, trk_v );     
      int fitcellv = dut.GetVCellFromCoord( trk_u, trk_v );    
      
      _rootSensorID = dut.GetSensorID();             
      _rootPlaneNumber = ipl;
      _rootFitMomentum = trk_mom;        
      _rootFitU = trk_u;                   
      _rootFitV = trk_v;                     
      _rootFitdUdW = trk_tu;        
      _rootFitdVdW = trk_tv;            
      _rootFitErrorU = TMath::Sqrt( TE.GetState().GetCov()(2,2) );           
      _rootFitErrorV = TMath::Sqrt( TE.GetState().GetCov()(3,3) );                  
      _rootFitCellU = fitcellu;                     
      _rootFitCellV = fitcellv;                
      _rootFitCellUCenter = dut.GetPixelCenterCoordU( fitcellv, fitcellu );              
      _rootFitCellVCenter = dut.GetPixelCenterCoordV( fitcellv, fitcellu );          
      _rootTrackPixelType = dut.GetPixelType(fitcellv, fitcellu);         

      // Set defaults for hit varibles
      _rootHasHit = -1;            
      _rootHitQuality = -1;             
      _rootHitLocalChi2 = -1;
      _rootPullResidualU = -1;    
      _rootPullResidualV = -1;      
      _rootHitU = -1;                  
      _rootHitV = -1;               
      _rootHitClusterCharge = -1;    
      _rootHitSeedCharge = -1;     
      _rootHitSeedPixelType = -1;     
      _rootHitSize = -1;              
      _rootHitSizeU = -1;               
      _rootHitSizeV = -1;              
      _rootHitCellU = -1;             
      _rootHitCellV = -1;             
      _rootHitSeedCellU = -1;            
      _rootHitSeedCellV = -1;           
	
       
      // Skip sensor w/o measurment
      if ( !TE.HasHit() ) continue;  
       
      // Get pixel residuals 
      double hit_u = TE.GetHit().GetCoord()[0]; // mm 
      double hit_v = TE.GetHit().GetCoord()[1]; // mm
               
      double pull_u = (hit_u - trk_u) / TMath::Sqrt( TE.GetState().GetCov()(2,2) + TE.GetHit().GetCov()(0,0) ) ; 
      double pull_v = (hit_v - trk_v) / TMath::Sqrt( TE.GetState().GetCov()(3,3) + TE.GetHit().GetCov()(1,1) ) ;  
      
      PixelCluster Cluster = TE.GetHit().GetCluster(); 

      _rootHasHit = 0;            
      _rootHitQuality = TE.GetHit().GetQuality();           
      _rootHitLocalChi2 = TrackFitter.GetPredictedChi2(TE.GetState().GetPars(), TE.GetState().GetCov(), TE.GetHit());             
      _rootPullResidualU = pull_u;    
      _rootPullResidualV = pull_v;      
      _rootHitU = hit_u;                  
      _rootHitV = hit_v;               
      _rootHitClusterCharge = Cluster.getCharge();    
      _rootHitSeedCharge = Cluster.getSeedCharge();         
      _rootHitSize = Cluster.getSize();             
      _rootHitSizeU = Cluster.getUSize();              
      _rootHitSizeV = Cluster.getVSize();     
      _rootHitCellU = dut.GetUCellFromCoord( hit_u, hit_v );               
      _rootHitCellV = dut.GetVCellFromCoord( hit_u, hit_v );                  
      _rootHitSeedCellU = Cluster.getUSeed();         
      _rootHitSeedCellV = Cluster.getVSeed();        
      _rootHitSeedPixelType = dut.GetPixelType(_rootHitSeedCellV, _rootHitSeedCellU);           
      
      
      // Fill tree with set variables 
      _rootFile->cd("");
      _rootHitTree->Fill();  
                  
    } // end sensor loop   
    
                       
  } // end track loop 
}

//
// Method called after each event to check the data processed
//
void TrackFitAnalyzer::check( LCEvent * )
{
}

//
// Method called after all data processing
//
void TrackFitAnalyzer::end()
{

 
  
  _rootFile->Write();
  _rootFile->Close();
   
  delete _rootFile;
   
  streamlog_out ( MESSAGE3 ) << endl;
  streamlog_out ( MESSAGE3 ) << "Successfully finished" << endl;
  
  // CPU time end
  _timeCPU = clock()/1000 - _timeCPU;
   
  // Print message
  streamlog_out(MESSAGE3) << std::endl
                           << " "
                           << "Time per event: "
                           << std::setiosflags(std::ios::fixed | std::ios::internal )
                           << std::setprecision(3)
                           << _timeCPU/_nEvt
                           << " ms"
                           << std::endl
                           << std::setprecision(3)
                           << std::endl
                           << " "
                           << "Processor succesfully finished!"
                           << std::endl;
  
}


//
// Method printing processor parameters
//
void TrackFitAnalyzer::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "TrackFitAnalyzer Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}


void TrackFitAnalyzer::bookHistos()
{   
  
  _rootFile = new TFile( _rootFileName.c_str(),"recreate");
  _rootFile->cd("");
  
  // 
  // Hit Tree  
  _rootHitTree = new TTree("Hit","Hit info");
  _rootHitTree->Branch("iRun"            ,&_rootRunNumber        ,"iRun/I");
  _rootHitTree->Branch("iEvt"            ,&_rootEventNumber      ,"iEvt/I");
  _rootHitTree->Branch("trackChi2"       ,&_rootTrackChi2      ,"trackChi2/D");
  _rootHitTree->Branch("trackNdof"       ,&_rootTrackNDF       ,"trackNdof/I");
  _rootHitTree->Branch("trackNHits"      ,&_rootTrackNHits     ,"trackNHits/I");   
  _rootHitTree->Branch("hasRefHit"       ,&_rootTrackWithRefHit,"hasRefHit/I");    
  _rootHitTree->Branch("sensorID"        ,&_rootSensorID       ,"sensorID/I");
  _rootHitTree->Branch("ipl"             ,&_rootPlaneNumber    ,"ipl/I");
  _rootHitTree->Branch("momentum"        ,&_rootFitMomentum      ,"momentum/D");    
  _rootHitTree->Branch("u_fit"           ,&_rootFitU             ,"u_fit/D");
  _rootHitTree->Branch("v_fit"           ,&_rootFitV             ,"v_fit/D"); 
  _rootHitTree->Branch("dudw_fit"        ,&_rootFitdUdW          ,"dudw_fit/D");
  _rootHitTree->Branch("dvdw_fit"        ,&_rootFitdVdW          ,"dvdw_fit/D");    
  _rootHitTree->Branch("u_fiterr"        ,&_rootFitErrorU        ,"u_fiterr/D");
  _rootHitTree->Branch("v_fiterr"        ,&_rootFitErrorV        ,"v_fiterr/D");   
  _rootHitTree->Branch("cellU_fit"       ,&_rootFitCellU           ,"cellU_fit/I");
  _rootHitTree->Branch("cellV_fit"       ,&_rootFitCellV           ,"cellV_fit/I");
  _rootHitTree->Branch("cellUCenter_fit" ,&_rootFitCellUCenter  ,"cellUCenter_fit/D");
  _rootHitTree->Branch("cellVCenter_fit" ,&_rootFitCellVCenter  ,"cellVCenter_fit/D");    
  _rootHitTree->Branch("pixeltype_track" ,&_rootTrackPixelType     ,"pixeltype_track/I");          
  _rootHitTree->Branch("hasHit"          ,&_rootHasHit        ,"hasHit/I");  
  _rootHitTree->Branch("clusterQuality"  ,&_rootHitQuality   ,"clusterQuality/I");
  _rootHitTree->Branch("localChi2"       ,&_rootHitLocalChi2       ,"localChi2/D"); 
  _rootHitTree->Branch("pull_resu"       ,&_rootPullResidualU    ,"pull_resu/D");
  _rootHitTree->Branch("pull_resv"       ,&_rootPullResidualV    ,"pull_resv/D");   
  _rootHitTree->Branch("u_hit"           ,&_rootHitU             ,"u_hit/D");
  _rootHitTree->Branch("v_hit"           ,&_rootHitV             ,"v_hit/D");     
  _rootHitTree->Branch("clusterCharge"   ,&_rootHitClusterCharge ,"clusterCharge/D");
  _rootHitTree->Branch("seedCharge"      ,&_rootHitSeedCharge       ,"seedCharge/D");
  _rootHitTree->Branch("sizeU"           ,&_rootHitSizeU        ,"sizeU/I");
  _rootHitTree->Branch("sizeV"           ,&_rootHitSizeV        ,"sizeV/I");
  _rootHitTree->Branch("size"            ,&_rootHitSize         ,"size/I");
  _rootHitTree->Branch("cellU_hit"       ,&_rootHitCellU              ,"cellU_hit/I");
  _rootHitTree->Branch("cellV_hit"       ,&_rootHitCellV              ,"cellV_hit/I");
  _rootHitTree->Branch("cellU_seed"      ,&_rootHitSeedCellU          ,"cellU_seed/I");
  _rootHitTree->Branch("cellV_seed"      ,&_rootHitSeedCellV          ,"cellV_seed/I");
  _rootHitTree->Branch("pixeltype"       ,&_rootHitSeedPixelType     ,"pixeltype/I");                                 
   
}

} // Namespace
