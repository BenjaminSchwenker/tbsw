// DEPFETResoCalc Processor  
// 
// See DEPFETResoCalc.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "DEPFETResoCalc.h"

// DEPFETTrackTools includes
#include "PixelCluster.h"
#include "TBHit.h"
#include "Det.h"
#include "MatrixDecoder.h" 

// Include basic C
#include <iostream>
#include <limits>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <string>
#include <algorithm>

// Include LCIO classes
#include <UTIL/CellIDDecoder.h>
#include <IMPL/LCFlagImpl.h>

// Include CLHEP classes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/SymMatrix.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Random/RandGauss.h>
#include <CLHEP/Random/RandPoisson.h>
#include <CLHEP/Random/RandFlat.h>

// Used namespaces
using namespace std; 
using namespace CLHEP; 
using namespace lcio ;
using namespace marlin ;

namespace depfet {

//
// Instantiate this object
//
DEPFETResoCalc aDEPFETResoCalc ;

//
// Constructor
//
DEPFETResoCalc::DEPFETResoCalc() : Processor("DEPFETResoCalc")
{

// Processor description
  _description = "DEPFETResoCalc: Calculate position resolution for pixel detector digitizers" ;
   
//
// Input collections 
   
  registerInputCollection( LCIO::TRACKERHIT,
                           "RecoHitCollection" ,
                           "Name of reco hit collection"  ,
                           _recoHitColName ,
                           std::string("rechit") ) ;
   
  registerInputCollection( LCIO::TRACKERHIT,
                           "TruthHitCollection" ,
                           "Name of truth hit collection"  ,
                           _truthHitColName ,
                           std::string("truthhit") ) ;
         
  // Processor parameters:
  
  registerProcessorParameter ("DUTPlane",
                              "Plane number of DUT along the beam line",
                              _idut,  
                              static_cast < int > (3));
  
  registerProcessorParameter ("RunNumber",
                              "Monte Carlo run number",
                              _runnumber,  
                              static_cast < int > (99));

  registerProcessorParameter ("TelescopeErrorU",
                              "Telescope pointing resolution [mm]",
                              _telerrorU,  
                              static_cast < double > (0.0));

  registerProcessorParameter ("TelescopeErrorV",
                              "Telescope pointing resolution [mm]",
                              _telerrorV,  
                              static_cast < double > (0.0));
     
  registerProcessorParameter ("HitDistMax",
                              "Maximum hit dist for matching [mm]",
                              _hitdistmax,  
                              static_cast < double > (0.05));
                               
  registerProcessorParameter( "RootFileName",
                              "Output root file name",
                              _rootFileName, 
                              std::string("reso_histos.root"));



   
}

//
// Method called at the beginning of data processing
//
void DEPFETResoCalc::init() {
   
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   _timeCPU = clock()/1000 ;
   
   // Initialize all the counters
   _noOfTruthHits = 0;       
   _noOfRecoHits = 0;     
   _noOfMatchedHits = 0;     
   
   // Print set parameters
   printProcessorParams();
   
   // Read detector constants from gear file
   _detector.ReadGearConfiguration();    
   
   // Load DUT module    
   Det & dut = _detector.GetDet(_idut); 
          
   // Print out geometry information  
   streamlog_out ( MESSAGE3 )  << "D.U.T. plane  ID = " << dut.GetDAQID()
                               << "  at position = " << _idut 
                               << endl << endl;
       
         
   bookHistos();
   
}

//
// Method called for each run
//
void DEPFETResoCalc::processRunHeader(LCRunHeader * run)
{

// Print run number
   streamlog_out(MESSAGE3) << "Processing run: "
                           << (_runnumber)
                           << std::endl << std::endl;
   
   _nRun++ ;

}

//
// Method called for each event
//
void DEPFETResoCalc::processEvent(LCEvent * evt)
{
  
  // Print event number
  if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                << (evt->getEventNumber())
                                                                << std::endl << std::endl;
  
  streamlog_out(MESSAGE2) << "Events processed: " << (evt->getEventNumber())
                                                   << std::endl << std::endl;
  
  _nEvt ++ ;
     
  // Load DUT module    
  Det & dut = _detector.GetDet(_idut);   
    
  // Decoding of DUT matrix
  MatrixDecoder matrixDecoder(dut.GetNColumns(), dut.GetNRows()); 
  
      
  //
  // Get truth hit collection 
  //
  
  LCCollection* truthhitcol = NULL;
  bool isTruthOk = true;   
  try {
    truthhitcol = evt->getCollection( _truthHitColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out(MESSAGE2) << "Not able to get collection "
                            << _truthHitColName
                            << " from event " << evt->getEventNumber()
                            << " in lcio run " << evt->getRunNumber()  << endl << endl;   
    
    isTruthOk = false; 
  }  
  
  // 
  // Get reco hit collection 
  // 
  
  LCCollection* recohitcol = NULL;
  bool isRecoOk=true;
  try {
    recohitcol = evt->getCollection( _recoHitColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out(MESSAGE2) << "Not able to get collection "
                            << _recoHitColName
                            << " from event " << evt->getEventNumber()
                            << " in lcio run " << evt->getRunNumber() << endl << endl;
    isRecoOk=false;
  } 
  
  // Store truth hits
  std::vector<TBHit> TruthHitStore;    
  
  
  int nTruthHits = 0;
  if(isTruthOk) 
    nTruthHits = truthhitcol->getNumberOfElements();
  
  streamlog_out(MESSAGE2) << "Total of " << nTruthHits << " truth hits found" << endl;      
  
  for(int ihit=0; ihit< nTruthHits ; ++ihit)
  {
    
    TrackerHitImpl * thit = dynamic_cast<TrackerHitImpl*>( truthhitcol->getElementAt(ihit) ) ;
    TBHit TruthHit ( thit  );  

    // We have to find plane number of the hit 
    int daqid = TruthHit.GetDAQID();   
    int ipl = _detector.GetPlaneNumber(daqid);      
           
    if( dut.GetPlaneNumber() == ipl  )
    {   
      
      // Smear truth hit according to telescope error 
      TruthHit.GetCoord()[0][0] += double(RandGauss::shoot(0,_telerrorU));
      TruthHit.GetCoord()[1][0] += double(RandGauss::shoot(0,_telerrorV));
      
      streamlog_out(MESSAGE2) << " truth hit at position" 
                              << "   U [mm]= " << TruthHit.GetCoord()[0][0]
                              << "   V [mm]= " << TruthHit.GetCoord()[1][0] 
                              << endl;
       
      TruthHitStore.push_back( TruthHit );  
    }  
            
  } 
  
  _noOfTruthHits += TruthHitStore.size();   
  
  // Store reco hits
  std::vector<TBHit> RecoHitStore;     
  
  
  int nRecoHits = 0;
  if(isRecoOk) 
    nRecoHits = recohitcol->getNumberOfElements();
    
  streamlog_out(MESSAGE2) << "Total of " << nRecoHits << " reco hits found"  << endl; 
  
  for(int ihit=0; ihit< nRecoHits ; ihit++)
  {
    
    TrackerHitImpl * rhit = dynamic_cast<TrackerHitImpl*>( recohitcol->getElementAt(ihit) ) ;
    TBHit RecoHit ( rhit  );  
     
    // We have to find plane number of the hit 
    int daqid = RecoHit.GetDAQID(); 
    int ipl = _detector.GetPlaneNumber(daqid);         
           
    if( dut.GetPlaneNumber() == ipl  )
    {   
      streamlog_out(MESSAGE2) << " reco hit at position" 
                              << "   U [mm]= " << RecoHit.GetCoord()[0][0]
                              << "   V [mm]= " << RecoHit.GetCoord()[1][0] 
                              << endl;
       
      RecoHitStore.push_back( RecoHit );  
    }  
    
  }
  
  _noOfRecoHits += RecoHitStore.size(); 
  
  // 
  // Match reco hits with truth hits 
   
  int nMatch=0;	
  
  
  // Tagging each match in both ways
  vector<int> reco2truth(RecoHitStore.size(), -1);
  vector<int> truht2reco(TruthHitStore.size(), -1); 
  
  
  
  // Continue matching hits until all hits are matched 
  // or no match possible. 
  double distmin=numeric_limits<double >::max();
  
  do{
    int besttruthhit=-1;
    int bestrecohit=-1;
    
    distmin=numeric_limits<double >::max();
    
    // Find hit pairs with minimum (euclidean) distance. 
     
    for(int itruth=0;itruth<(int)TruthHitStore.size(); ++itruth)
    {
      // Check if matched to reco
      if (truht2reco[itruth] >= 0) continue;  
      
      for(int ireco=0; ireco< (int)RecoHitStore.size(); ++ireco)
      {
        // Check if matched to truth
        if (reco2truth[ireco] >= 0) continue;  
        
        TBHit& RecoHit = RecoHitStore[ireco];
        TBHit& TruthHit = TruthHitStore[itruth];
  
        // Compute the distance hit to track  
        double uhit = RecoHit.GetCoord()[0][0];
        double vhit = RecoHit.GetCoord()[1][0];
        double utrk = TruthHit.GetCoord()[0][0];
        double vtrk = TruthHit.GetCoord()[1][0];
        double hitdist = std::abs(utrk-uhit) + std::abs(vtrk-vhit);
        
        if(hitdist<distmin)
        {
          distmin=hitdist;
          besttruthhit=itruth;
          bestrecohit=ireco;
        }
      }
    }
    
    streamlog_out(MESSAGE2) << "In matching loop: best truth hit " << besttruthhit << " to best reco hit " << bestrecohit << endl; 
    streamlog_out(MESSAGE2) << "  distmin: " <<  distmin  << endl; 
    
    // Check iff best match is good enough
    if( distmin < _hitdistmax  )
    {   
      nMatch++;
      reco2truth[bestrecohit] = besttruthhit;
      truht2reco[besttruthhit] = bestrecohit;   
    } 
     
  } // End of matching loop 
  
  while( distmin < _hitdistmax );
  
  
  _noOfMatchedHits += nMatch; 
  streamlog_out(MESSAGE2) << nMatch << " Reco hits matched to truth hits" << endl;
   
  // Loop over reco hits and fill reco tree 
  for(int ihit=0;ihit<(int)RecoHitStore.size(); ++ihit)
  {
    _rootRunNumber = _runnumber;  
    _rootEventNumber = evt->getEventNumber();  
    _rootDetectorID = dut.GetDAQID();       
    _rootStartGate = 0;
    _rootNTelTracks = nTruthHits; 
    _rootNDUTHits = (int)RecoHitStore.size(); 
    
    TBHit& hit = RecoHitStore[ihit];
    
    _rootHitU = hit.GetCoord()[0][0];         
    _rootHitV = hit.GetCoord()[1][0];   
    
    _rootHitCol = dut.GetColumnFromCoord( _rootHitU, _rootHitV );  
    _rootHitRow = dut.GetRowFromCoord( _rootHitU, _rootHitV );  
 
    // Cluster shape variables   
    PixelCluster Cluster = hit.GetCluster();
       
    _rootHitQuality = 0; 
    _rootHitCharge = Cluster.getCharge() ; 
    _rootHitSeedCharge = Cluster.getSeedCharge() ; 
    _rootHitSize = Cluster.getSize();  
    _rootHitSizeCol = Cluster.getUSize();     
    _rootHitSizeRow = Cluster.getVSize();     
 
    // Add variables for matched track 
    if ( reco2truth[ihit] >= 0 ) {  
      
      _rootHitHasTrack = 0;  // matched        
      
      TBHit& thit = TruthHitStore[ reco2truth[ihit] ]; 
      double ufit = thit.GetCoord()[0][0];   
      double vfit = thit.GetCoord()[1][0];   
       
      // Get readout channels  
      int fitcol = dut.GetColumnFromCoord( ufit, vfit );     
      int fitrow = dut.GetRowFromCoord( ufit, vfit );           
       
      _rootHitFitdUdW = 0;     
      _rootHitFitdVdW = 0;    
      _rootHitFitU = ufit;      
      _rootHitFitV = vfit;          

      _rootHitFitErrorU  = _telerrorU;    
      _rootHitFitErrorV  = _telerrorV;  
      _rootHitPullResidualU = (hit.GetCoord()[0][0] - ufit) / TMath::Sqrt( _telerrorU*_telerrorU + hit.GetCov()[0][0] ) ;   
      _rootHitPullResidualV = (hit.GetCoord()[1][0] - vfit) / TMath::Sqrt( _telerrorV*_telerrorV + hit.GetCov()[1][1] ) ;  
                                         
      _rootHitFitCol = fitcol;      
      _rootHitFitRow = fitrow;    
      _rootHitFitPixU = dut.GetPixelCenterCoordU( fitrow, fitcol ); 
      _rootHitFitPixV = dut.GetPixelCenterCoordV( fitrow, fitcol );                                        
      _rootHitTrackChi2 = 0; 
      _rootHitTrackNDF =  8; // assume 6 pixel hits from a reference track
      _rootHitLocalChi2 = 0;   
      _rootStartGate = 0; 
      
    } else {
    
      _rootHitHasTrack = -1; // no match      
      // These are dummy values, always query hasTrack      
      _rootHitFitU = -1;           
      _rootHitFitV = -1;           
      _rootHitFitdUdW = -1;      
      _rootHitFitdVdW = -1;    
      _rootHitFitErrorU = -1;  
      _rootHitFitErrorV = -1;      
      _rootHitPullResidualU = -1;   
      _rootHitPullResidualV = -1;                        
      _rootHitFitCol = -1;      
      _rootHitFitRow = -1;      
      _rootHitFitPixU = -1;  
      _rootHitFitPixV = -1;                                       
      _rootHitTrackChi2 = -1;  
      _rootHitTrackNDF = -1; 
      _rootHitLocalChi2 = -1;  
      
    }
    
    // Fill tree with set variables 
    _rootFile->cd("");
    _rootHitTree->Fill();
  
  }
  
  for(int itrk=0;itrk<(int)TruthHitStore.size(); ++itrk)
  {
  
    _rootRunNumber = _runnumber;  
    _rootEventNumber = evt->getEventNumber();  
    _rootDetectorID = dut.GetDAQID();   
    _rootStartGate = 0;   
    _rootNTelTracks = nTruthHits; 
    _rootNDUTHits = (int)RecoHitStore.size();      
     
    // Read matched truth hit 
    TBHit& thit = TruthHitStore[itrk];
    
    double ufit = thit.GetCoord()[0][0];   
    double vfit = thit.GetCoord()[1][0]; 
           
    // Get readout channels  
    int fitcol = dut.GetColumnFromCoord( ufit, vfit );     
    int fitrow = dut.GetRowFromCoord( ufit, vfit );   
           
            
    _rootTrackFitdUdW = 0;     
    _rootTrackFitdVdW = 0;    
    _rootTrackFitU = ufit;           
    _rootTrackFitV = vfit;    
    _rootTrackFitCol = fitcol;      
    _rootTrackFitRow = fitrow;    
    _rootTrackFitPixU = dut.GetPixelCenterCoordU( fitrow, fitcol ); 
    _rootTrackFitPixV = dut.GetPixelCenterCoordV( fitrow, fitcol );                                        
    _rootTrackChi2 = 0; 
    _rootTrackNDF = 0;   
         
    _rootTrack1x1Charge = 0; 
    _rootTrack3x3Charge = 0;   
    _rootTrack1x1Quality = 0;    
    _rootTrack3x3Quality = 0;  
            
    if ( truht2reco[itrk] >= 0 ) {
      
      _rootTrackHasHit = 0;  // match   
      _rootTrackSeedCharge = -1;
       
      TBHit& hit = RecoHitStore[ truht2reco[itrk] ];  

      PixelCluster Cluster = hit.GetCluster();
      _rootTrackSeedCharge = Cluster.getSeedCharge() ; 
      
       
    } else {
      _rootTrackHasHit = -1; // no match
      _rootTrackSeedCharge = -1;
    }
     
    _rootStartGate = 0; 
    
    _rootFile->cd("");
    _rootTrackTree->Fill(); 
  
  }
  
  // Fill event tree
  _rootRunNumber = _runnumber;  
  _rootNDUTHits = (int)RecoHitStore.size(); 
  _rootEventNDUTTracks = (int) TruthHitStore.size(); 
  _rootNTelTracks = (int) TruthHitStore.size(); 
  _rootEventNMatched = nMatch; 
  _rootStartGate = 0;
  
  _rootFile->cd("");
  _rootEventTree->Fill();  
   
  return;
}


//
// Method called after each event to check the data processed
//
void DEPFETResoCalc::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void DEPFETResoCalc::end()
{

   // Print the summer
   streamlog_out(MESSAGE3) << endl << endl 
                           << "Total number of processed events:     " << setiosflags(ios::right) << _nEvt 
                           << resetiosflags(ios::right) << endl
                           << "Total number of matched hits:         " << setiosflags(ios::right) << _noOfMatchedHits
                           << resetiosflags(ios::right) << endl
                           << "Total number of truth hits:           " << setiosflags(ios::right) << _noOfTruthHits
                           << resetiosflags(ios::right) << endl
                           << "Total number of reco hits:            " << setiosflags(ios::right) << _noOfRecoHits
                           << resetiosflags(ios::right) << endl
                           << endl << endl; 
      
   
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

    
   // Close root  file
   _rootFile->Write();
   _rootFile->Close();
   delete _rootFile;

}



//
// Method printing processor parameters
//
void DEPFETResoCalc::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "DEPFETResoCalc Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}


 // ROOT_OUTPUT
void DEPFETResoCalc::bookHistos()
{   
   
   _rootFile = new TFile( _rootFileName.c_str(),"recreate");
   _rootFile->cd("");
   
   // 
   // Hit Tree  
   _rootHitTree = new TTree("Hit","Hit info");
   _rootHitTree->Branch("iRun"            ,&_rootRunNumber        ,"iRun/I");
   _rootHitTree->Branch("iEvt"            ,&_rootEventNumber      ,"iEvt/I");
   _rootHitTree->Branch("det"             ,&_rootDetectorID       ,"det/I");
   _rootHitTree->Branch("nTelTracks"      ,&_rootNTelTracks       ,"nTelTracks/I"); 
   _rootHitTree->Branch("startgate"       ,&_rootStartGate        ,"startgate/I");
   _rootHitTree->Branch("nduthits"        ,&_rootNDUTHits         ,"nduthits/I");
   
   
   _rootHitTree->Branch("clusterQuality"  ,&_rootHitQuality   ,"clusterQuality/I");
   _rootHitTree->Branch("u_hit"           ,&_rootHitU             ,"u_hit/D");
   _rootHitTree->Branch("v_hit"           ,&_rootHitV             ,"v_hit/D");     
   _rootHitTree->Branch("clusterCharge"   ,&_rootHitCharge    ,"clusterCharge/D");
   _rootHitTree->Branch("seedCharge"      ,&_rootHitSeedCharge       ,"seedCharge/D");
   _rootHitTree->Branch("sizeCol"         ,&_rootHitSizeCol   ,"sizeCol/I");
   _rootHitTree->Branch("sizeRow"         ,&_rootHitSizeRow   ,"sizeRow/I");
   _rootHitTree->Branch("size"            ,&_rootHitSize      ,"size/I");
   _rootHitTree->Branch("hasTrack"        ,&_rootHitHasTrack         ,"hasTrack/I");   
   _rootHitTree->Branch("u_fit"           ,&_rootHitFitU             ,"u_fit/D");
   _rootHitTree->Branch("v_fit"           ,&_rootHitFitV             ,"v_fit/D"); 
   _rootHitTree->Branch("dudw_fit"        ,&_rootHitFitdUdW          ,"dudw_fit/D");
   _rootHitTree->Branch("dvdw_fit"        ,&_rootHitFitdVdW          ,"dvdw_fit/D");         
   _rootHitTree->Branch("col_fit"         ,&_rootHitFitCol           ,"col_fit/I");
   _rootHitTree->Branch("row_fit"         ,&_rootHitFitRow           ,"row_fit/I");
   _rootHitTree->Branch("col_hit"         ,&_rootHitCol           ,"col_hit/I");
   _rootHitTree->Branch("row_hit"         ,&_rootHitRow           ,"row_hit/I");
   _rootHitTree->Branch("u_pixel"         ,&_rootHitFitPixU          ,"u_pixel/D");
   _rootHitTree->Branch("v_pixel"         ,&_rootHitFitPixV          ,"v_pixel/D");                                      
   _rootHitTree->Branch("chi2"            ,&_rootHitTrackChi2      ,"chi2/D");
   _rootHitTree->Branch("ndof"            ,&_rootHitTrackNDF       ,"ndof/I");
   _rootHitTree->Branch("chi2pred"        ,&_rootHitLocalChi2        ,"chi2pred/D");    
   _rootHitTree->Branch("u_fiterr"        ,&_rootHitFitErrorU        ,"u_fiterr/D");
   _rootHitTree->Branch("v_fiterr"        ,&_rootHitFitErrorV        ,"v_fiterr/D");   
   _rootHitTree->Branch("pull_resu"       ,&_rootHitPullResidualU    ,"pull_resu/D");
   _rootHitTree->Branch("pull_resv"       ,&_rootHitPullResidualV    ,"pull_resv/D");  
   
    
   // 
   // Track Tree 
   _rootTrackTree = new TTree("Track","Track info");
   _rootTrackTree->Branch("iRun"            ,&_rootRunNumber      ,"iRun/I");
   _rootTrackTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
   _rootTrackTree->Branch("det"             ,&_rootDetectorID     ,"det/I");
   _rootTrackTree->Branch("startgate"       ,&_rootStartGate      ,"startgate/I");
   _rootTrackTree->Branch("nTelTracks"      ,&_rootNTelTracks     ,"nTelTracks/I"); 
   _rootTrackTree->Branch("nduthits"           ,&_rootNDUTHits          ,"nduthits/I");
   
   _rootTrackTree->Branch("hasHit"          ,&_rootTrackHasHit         ,"hasHit/I");
   _rootTrackTree->Branch("seedCharge"      ,&_rootTrackSeedCharge     ,"seedCharge/D");                                                        
   _rootTrackTree->Branch("1x1Charge"       ,&_rootTrack1x1Charge      ,"1x1Charge/D");
   _rootTrackTree->Branch("3x3Charge"       ,&_rootTrack3x3Charge      ,"3x3Charge/D");
   _rootTrackTree->Branch("1x1Quality"      ,&_rootTrack1x1Quality     ,"1x1Quality/I");
   _rootTrackTree->Branch("3x3Quality"      ,&_rootTrack3x3Quality     ,"3x3Quality/I");
   _rootTrackTree->Branch("u_fit"           ,&_rootTrackFitU           ,"u_fit/D");
   _rootTrackTree->Branch("v_fit"           ,&_rootTrackFitV           ,"v_fit/D");
   _rootTrackTree->Branch("dudw_fit"        ,&_rootTrackFitdUdW        ,"dudw_fit/D");
   _rootTrackTree->Branch("dvdw_fit"        ,&_rootTrackFitdVdW        ,"dvdw_fit/D");
   _rootTrackTree->Branch("col_fit"         ,&_rootTrackFitCol         ,"col_fit/I");
   _rootTrackTree->Branch("row_fit"         ,&_rootTrackFitRow         ,"row_fit/I");
   _rootTrackTree->Branch("u_pixel"         ,&_rootTrackFitPixU        ,"u_pixel/D");
   _rootTrackTree->Branch("v_pixel"         ,&_rootTrackFitPixV        ,"v_pixel/D");
   _rootTrackTree->Branch("chi2"            ,&_rootTrackChi2           ,"chi2/D");
   _rootTrackTree->Branch("ndof"            ,&_rootTrackNDF            ,"ndof/I");
   
   // 
   // Event Summay Tree 
   _rootEventTree = new TTree("Event","Event info");
   _rootEventTree->Branch("iRun"            ,&_rootRunNumber      ,"iRun/I");
   _rootEventTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
   _rootEventTree->Branch("det"             ,&_rootDetectorID     ,"det/I");
   _rootEventTree->Branch("startgate"       ,&_rootStartGate      ,"startgate/I");       
   _rootEventTree->Branch("nTelTracks"      ,&_rootNTelTracks     ,"nTelTracks/I"); 
   _rootEventTree->Branch("nduthits"        ,&_rootNDUTHits          ,"nduthits/I");
   _rootEventTree->Branch("nDUTTracks"      ,&_rootEventNDUTTracks     ,"nDUTTracks/I");
   _rootEventTree->Branch("matched"         ,&_rootEventNMatched       ,"matched/I");
   
}

} // Namespace



