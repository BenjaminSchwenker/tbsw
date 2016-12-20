// MCTreeProducer Processor  
// 
// See MCTreeProducer.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "MCTreeProducer.h"

// TBTools includes
#include "TBHit.h"
#include "PixelCluster.h"
#include "Det.h"
#include "Utilities.h" 
#include "PhysicalConstants.h"

// Include basic C
#include <iostream>
#include <limits>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <string>
#include <algorithm>

// Include LCIO classes
#include "lcio.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>

// Include CLHEP classes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Vector/ThreeVector.h>

// Used namespaces
using namespace std; 
using namespace CLHEP; 
using namespace lcio ;
using namespace marlin ;

namespace depfet {

  //
  // Instantiate this object
  //
  MCTreeProducer aMCTreeProducer ;
  
  //
  // Constructor
  //
  MCTreeProducer::MCTreeProducer() : Processor("MCTreeProducer")
  {
    
    // Processor description
    _description = "MCTreeProducer: Matching SimtrackerHits to clusters and moving data into output TTree(s)" ;
    
    //
    // Input collections 
    
    registerInputCollection( LCIO::SIMTRACKERHIT,
                             "SimTrackerHitCollection" ,
                             "Name of collection with simulated hits"  ,
                             _simHitColName ,
                             std::string("SimTrackerHits") ) ;
    
    registerInputCollection( LCIO::TRACKERHIT,
                             "HitCollection" ,
                             "Name of hit collection"  ,
                             _hitColName ,
                             std::string("hit") ) ;
     
    
    // Processor parameters:
    
    registerProcessorParameter ("AlignmentDBFileName",
                                "This is the name of the LCIO file with the alignment constants (add .slcio)",
                                _alignmentDBFileName, static_cast< string > ( "alignmentDB.slcio" ) );     
    
    registerProcessorParameter ("DUTPlane",
                                "Plane number of DUT along the beam line",
                                _idut,  static_cast < int > (3));
       
    registerProcessorParameter ("MaxResidualU",
                                "Maximum u residual for matching simHits to hits [mm]. Put -1 to deactivate cut.",
                                _maxResidualU,  static_cast < double > (0.2));
    
    registerProcessorParameter ("MaxResidualV",
                                "Maximum v residual for matching simHits to hits [mm]. Put -1 to deactivate cut.",
                                _maxResidualV,  static_cast < double > (0.2));
    
    registerProcessorParameter( "RootFileName",
                                "Output root file name",
                                _rootFileName, std::string("histos.root"));  
    
    
  }
  
  //
  // Method called at the beginning of data processing
  //
  void MCTreeProducer::init() {
   
    // Initialize variables
    _nRun = 0 ;
    _nEvt = 0 ;
    _timeCPU = clock()/1000 ;
    
    // Print set parameters
    printProcessorParams();
    
    // Read detector constants from gear file
    _detector.ReadGearConfiguration();  
    
    // Read alignment data base file 
    _detector.ReadAlignmentDB( _alignmentDBFileName );    
    
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
  void MCTreeProducer::processRunHeader(LCRunHeader * run)
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
  void MCTreeProducer::processEvent(LCEvent * evt)
  {
    
    _nEvt ++ ;
     
    // Load DUT module    
    Det & dut = _detector.GetDet(_idut);   
        
    //
    // Get simhit collection
    //
  
    LCCollection* simHitCol = NULL;
    int nSimHit = 0;   
    try {
      simHitCol = evt->getCollection( _simHitColName ) ;
      nSimHit = simHitCol->getNumberOfElements();
    } catch (lcio::DataNotAvailableException& e) {
      streamlog_out(MESSAGE2) << "Not able to get collection "
                              << _simHitColName
                              << " from event " << evt->getEventNumber()
                              << " in run " << evt->getRunNumber()  << endl << endl;   
    }  
   
    streamlog_out(MESSAGE2) << "Total of " << nSimHit << " simHits in collection " << _simHitColName << endl; 
    
    // 
    // Get hit collection 
    // 
    
    LCCollection* hitCol = 0;
    int nHit = 0; 
    try {
      hitCol = evt->getCollection( _hitColName ) ;
      nHit = hitCol->getNumberOfElements();
    } catch (lcio::DataNotAvailableException& e) {
      streamlog_out(MESSAGE2) << "Not able to get collection "
                              << _hitColName
                              << " from event " << evt->getEventNumber()
                              << " in run " << evt->getRunNumber() << endl << endl;
    } 
    
    streamlog_out(MESSAGE2) << "Total of " << nHit << " hit(s) in collection " << _hitColName << endl;
    
    
    // Read hits and simHits 
    // ----------------------------------
    
    std::vector<SimTrackerHit*> SimHitStore;
    std::vector<TBHit> HitStore;    
 
    CellIDDecoder<SimTrackerHit> cellIDDec(simHitCol);
     
    for(int i=0; i< nSimHit ; i++)
    {  
      // Retrieve simtrackerhit 
      SimTrackerHit * simHit = dynamic_cast<SimTrackerHit*> (simHitCol->getElementAt(i));
      
      // Set current - layer ID, ladder ID and sensor ID
      int sensorID = cellIDDec(simHit)["sensorID"];
      int ipl = _detector.GetPlaneNumber(sensorID);
      
      // Store all hits on the DUT module  
      if( dut.GetPlaneNumber() == ipl )
      {         
        streamlog_out(MESSAGE2) << " SimHit at plane " << ipl << " at: (" << simHit->getPosition()[0] << ", " << simHit->getPosition()[1] << ")" 
                                << endl;
        
        SimHitStore.push_back(simHit);    
      }                            
    } 
    
    
    for(int i=0; i< nHit ; i++)
    {
      // Built a TBHit
      TrackerHitImpl * lciohit = dynamic_cast<TrackerHitImpl*>( hitCol->getElementAt(i) ) ;
      TBHit RecoHit ( lciohit  );        
         
      // We have to find plane number of the hit 
      int sensorID = RecoHit.GetDAQID();      
      int ipl = _detector.GetPlaneNumber(sensorID);  
      
      // Store all hits on the DUT module  
      if( dut.GetPlaneNumber() == ipl )
      {         
        streamlog_out(MESSAGE2) << " RecoHit at plane " << ipl << " at: (" << RecoHit.GetCoord()[0][0] << ", " << RecoHit.GetCoord()[1][0] << ")" 
                                << endl;
        
        HitStore.push_back( RecoHit );  
      }
    } 
    
    streamlog_out(MESSAGE2) << "Total of " << SimHitStore.size() << " simHit(s) on DUT" << endl; 
    streamlog_out(MESSAGE2) << "Total of " << HitStore.size() << " hit(s) on DUT" << endl; 
    
    
    
    // Record for each hit a matched simhit
    vector<int> hit2simhit(HitStore.size(), -1);
    
    // Continue matching until all hits are matched 
    // or no hit is close enough!!
    
    { 
      double distmin=numeric_limits<double >::max();
      int bestsimhit=-1;   
      int besthit=-1; 
      
      do{
        bestsimhit=-1;
        besthit=-1; 
        distmin=numeric_limits<double >::max();
        
        // Find hit/simhit pair with minimum chi2 distance.  
        for(int i=0; i< (int)HitStore.size() ; i++)
        {
            
          // If matched, skip hit 
          if (hit2simhit[i] >= 0) continue;
          
          for(int j=0;j<(int)SimHitStore.size(); j++)
          {
            SimTrackerHit * simHit = SimHitStore[j];
            double simHitPosU = simHit->getPosition()[0];
            double simHitPosV = simHit->getPosition()[1];
            
            TBHit& recoHit = HitStore[i];
            double hitPosU = recoHit.GetCoord()[0][0];
            double hitPosV = recoHit.GetCoord()[1][0];
             

            // Skip all DUT hits with too large residuum 
            if ( std::abs(hitPosU-simHitPosU) >= _maxResidualU && _maxResidualU > 0 ) continue;  
            if ( std::abs(hitPosV-simHitPosV) >= _maxResidualV && _maxResidualV > 0 ) continue; 
            
            // Finally, we will use a simple 2D distance to select best matching hit
            double hitdist = 0; 
            if ( _maxResidualU > 0 )  hitdist += std::abs(hitPosU-simHitPosU); 
            if ( _maxResidualV > 0 )  hitdist += std::abs(hitPosV-simHitPosV); 

                       
            if( hitdist<distmin )
            {
              distmin=hitdist;
              besthit=i;
              bestsimhit=j;
            }
          }
        }
        
        streamlog_out(MESSAGE2) << "In matching loop: best hit " << besthit << " to simhit " << bestsimhit << endl; 
        streamlog_out(MESSAGE2) << "  distmin: " <<  distmin  << endl; 
        
        // Check if a match was found
        if( bestsimhit>-1 &&  besthit>-1   )
        {   
          streamlog_out(MESSAGE2) << "  match found!!!"   << endl;
          hit2simhit[besthit] = bestsimhit;
        } 
  
      } // End of do loop of matching DUT hits to fitted positions
      while( bestsimhit>-1 &&  besthit>-1);
    }
    
    // Fill hit tree 
    
    _rootRunNumber = evt->getRunNumber();  
    _rootEventNumber = evt->getEventNumber();  
    _rootSensorID = dut.GetDAQID();   
    
    streamlog_out(MESSAGE2) << "Start fill hit tree" << endl; 
    
    for(int ihit=0;ihit<(int)HitStore.size(); ++ihit)
    {
      // Provide default values for all branches
      _rootClusterID = -1;              
      _rootClusterIDU = -1;             
      _rootClusterIDV = -1;             
      _rootClusterPosU = 0;          
      _rootClusterPosV = 0;          
      _rootClusterSigmaU = 0;        
      _rootClusterSigmaV = 0;        
      _rootClusterCharge = 0;        
      _rootSeedCharge = 0;           
      _rootClusterSize = 0;            
      _rootClusterSizeU = 0;            
      _rootClusterSizeV = 0;            
      _rootClusterStartU = 0;           
      _rootClusterStartV = 0;           
      _rootTrackMomentum = 0;         
      _rootTrackCharge = 0;        
      _rootTrackPosU = 0;           
      _rootTrackPosV = 0;                          
      _rootTrackdUdW = 0;           
      _rootTrackdVdW = 0;     
      
      if ( hit2simhit[ihit] >= 0 ) {        
        
        TBHit& hit = HitStore[ihit];
         
        _rootClusterPosU = hit.GetCoord()[0][0];         
        _rootClusterPosV = hit.GetCoord()[1][0];   
        _rootClusterSigmaU = hit.GetCov()[0][0];
        _rootClusterSigmaV = hit.GetCov()[1][1];
          
        PixelCluster Cluster = hit.GetCluster();
        
        _rootClusterCharge = Cluster.getCharge() ; 
        _rootSeedCharge = Cluster.getSeedCharge() ; 
        _rootClusterSize = Cluster.getSize();  
        _rootClusterSizeU = Cluster.getUSize();     
        _rootClusterSizeV = Cluster.getVSize();  
        _rootClusterStartU = Cluster.getUStart();
        _rootClusterStartV = Cluster.getVStart();   
        
        SimTrackerHit * simHit = SimHitStore[ hit2simhit[ihit] ]; 
        Hep3Vector momentum(simHit->getMomentum()[0],simHit->getMomentum()[1],simHit->getMomentum()[2]);
        
        std::cout << std::setiosflags(std::ios::fixed | std::ios::internal ) 
                  << std::setprecision(16)
                  << "momentum " << simHit->getMomentum()[0] << ", " << simHit->getMomentum()[1] << ", " << simHit->getMomentum()[2]
                  << std::resetiosflags(std::ios::showpos)
                  << std::setprecision(0) 
                  << std::endl;

        _rootTrackMomentum = momentum.mag();         
        _rootTrackCharge = simHit->getMCParticle()->getCharge();       
        _rootTrackPosU = simHit->getPosition()[0];          
        _rootTrackPosV = simHit->getPosition()[1];                           
        _rootTrackdUdW = momentum[0]/momentum[2];        
        _rootTrackdVdW = momentum[1]/momentum[2];    
         
        // Fill tree with set variables 
        _rootFile->cd("");
        _rootClusterTree->Fill();
      }
       
      
    }
    
    return;
  }


  //
  // Method called after each event to check the data processed
  //
  void MCTreeProducer::check( LCEvent * evt ) {}

  //
  // Method called after all data processing
  //
  void MCTreeProducer::end()
  {
    
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
  void MCTreeProducer::printProcessorParams() const 
  {
    
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "MCTreeProducer Development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
    
  }
  
  void MCTreeProducer::bookHistos()
  {   
     
    _rootFile = new TFile( _rootFileName.c_str(),"recreate");
    _rootFile->cd("");
     
    // 
    // Cluster Tree  
    _rootClusterTree = new TTree("Cluster","Cluster info");
    _rootClusterTree->Branch("iRun"            ,&_rootRunNumber              ,"iRun/I");
    _rootClusterTree->Branch("iEvt"            ,&_rootEventNumber            ,"iEvt/I");
    _rootClusterTree->Branch("sensorID"        ,&_rootSensorID               ,"sensorID/I");
    _rootClusterTree->Branch("clusterID"       ,&_rootClusterID              ,"clusterID/I"); 
    _rootClusterTree->Branch("clusterIDU"      ,&_rootClusterIDU             ,"clusterIDU/I"); 
    _rootClusterTree->Branch("clusterIDV"      ,&_rootClusterIDV             ,"clusterIDV/I"); 
    _rootClusterTree->Branch("clusterPosU"     ,&_rootClusterPosU            ,"clusterPosU/D");
    _rootClusterTree->Branch("clusterPosV"     ,&_rootClusterPosV            ,"clusterPosV/D");  
    _rootClusterTree->Branch("clusterSigmaU"   ,&_rootClusterSigmaU          ,"clusterSigmaU/D");
    _rootClusterTree->Branch("clusterSigmaV"   ,&_rootClusterSigmaV          ,"clusterSigmaV/D");   
    _rootClusterTree->Branch("clusterCharge"   ,&_rootClusterCharge          ,"clusterCharge/D");
    _rootClusterTree->Branch("seedCharge"      ,&_rootSeedCharge             ,"seedCharge/D");
    _rootClusterTree->Branch("clusterSizeU"    ,&_rootClusterSizeU           ,"clusterSizeU/I");
    _rootClusterTree->Branch("clusterSizeV"    ,&_rootClusterSizeV           ,"clusterSizeV/I");
    _rootClusterTree->Branch("clusterSize"     ,&_rootClusterSize            ,"clusterSize/I");
    _rootClusterTree->Branch("clusterStartU"   ,&_rootClusterStartU          ,"clusterSizeV/I");
    _rootClusterTree->Branch("clusterStartV"   ,&_rootClusterStartV          ,"clusterSize/I");
    _rootClusterTree->Branch("trackMomentum"   ,&_rootTrackMomentum          ,"trackMomentum/D"); 
    _rootClusterTree->Branch("trackCharge"     ,&_rootTrackCharge            ,"trackCharge/D"); 
    _rootClusterTree->Branch("trackPosU"       ,&_rootTrackPosU              ,"trackPosU/D");
    _rootClusterTree->Branch("trackPosV"       ,&_rootTrackPosV              ,"trackPosV/D"); 
    _rootClusterTree->Branch("trackDuDw"       ,&_rootTrackdUdW              ,"trackDuDw/D");
    _rootClusterTree->Branch("trackDvDw"       ,&_rootTrackdVdW              ,"trackDvDw/D");    
    
     
  }

} // Namespace



