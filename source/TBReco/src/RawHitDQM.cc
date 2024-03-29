// RawHitDQM - Marlin Processor
// 
// Produces DQM plots 
//
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 
#include "TBDetector.h"

#include "RawHitDQM.h"
#include "TBHit.h"
#include "Det.h"
#include "Utilities.h"
#include "HitFactory.h"
#include "PixelCluster.h"

// Include basic C
#include <limits>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <algorithm>

// Include LCIO classes
#include <lcio.h>

#include <TMath.h>

// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;


namespace depfet {

//
// Instantiate this object
//
RawHitDQM aRawHitDQM ;

//
// Constructor
//
RawHitDQM::RawHitDQM() : Processor("RawHitDQM")
{

// Processor description
   _description = "RawHitDQM: Producing DQM histos from hits";
   

//   
// First of all, we need to register the input/output collections
    
   vector< string > inputHitCollectionNameVecExample;
   inputHitCollectionNameVecExample.push_back( "hit" );
   
   registerInputCollections (LCIO::TRACKERHIT, "InputHitCollectionNameVec",
                            "Hit collection names",
                            _inputHitCollectionNameVec, inputHitCollectionNameVecExample );
   
   registerProcessorParameter("RootFileName", "Output root file name",
                               _rootFileName, std::string("tb_hits.root"));
   
                                 
}

//
// Method called at the beginning of data processing
//
void RawHitDQM::init() {
    
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   _timeCPU = clock()/1000;
                 
   // Print set parameters
   printProcessorParams();
   
   
   
   // Book root histos
   bookHistos();   
}

//
// Method called for each run
//
void RawHitDQM::processRunHeader(LCRunHeader * run)
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
void RawHitDQM::processEvent(LCEvent * evt)
{
   
  // Print event number
  if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                   << (evt->getEventNumber())
                                                                   << std::endl << std::endl;
   
    
  _nEvt ++ ;
   
  // Read input hit collections
  // ===========================
  
  // Hit factory sorts hits according to plane    
  HitFactory HitStore( TBDetector::GetInstance() );   
    
  for ( size_t iCol = 0 ; iCol < _inputHitCollectionNameVec.size(); ++iCol ) {
     
    try {
       
       LCCollectionVec * hitcol = dynamic_cast < LCCollectionVec * > (evt->getCollection( _inputHitCollectionNameVec.at( iCol ) ) );
              
       for ( int ihit = 0 ; ihit < (int) hitcol->size() ; ++ihit ) {
             
         // Built a TBHit     
         TrackerHitImpl * lciohit = dynamic_cast< TrackerHitImpl* > ( hitcol->getElementAt( ihit ) );
         TBHit RecoHit ( lciohit ); 
         HitStore.AddRecoHit(RecoHit);
                      
       } // End for hit loop      
          
    } catch (lcio::DataNotAvailableException& e) {
       streamlog_out ( MESSAGE2 ) << "Not able to get collection "
                                  << _inputHitCollectionNameVec.at( iCol )
                                  << "\nfrom event " << evt->getEventNumber()
                                  << " in run " << evt->getRunNumber()  << endl;
       
    }  
  }   
  
  streamlog_out ( MESSAGE2 ) << "Total of " << HitStore.GetNHits() << " hits." << endl;
  
  // Fill DQM histos 
  // ===============
  
  int nhits = 0; 
  
  for (int iplane = 0; iplane < TBDetector::GetInstance().GetNSensors(); ++iplane ) {
    
    std::string histoName;
    std::string planeNumber=(to_string( iplane ));
    nhits += HitStore.GetNHits(iplane); 
    
    histoName = "hnhits_sensor"+planeNumber;
    _histoMap[ histoName ]->Fill(HitStore.GetNHits(iplane));
    
    for (int ihit = 0; ihit < HitStore.GetNHits(iplane); ++ihit ) {
      
      const TBHit & Hit = HitStore.GetRecoHitFromID(ihit, iplane); 
      
      double hit_u = Hit.GetCoord()[0]; 
      double hit_v = Hit.GetCoord()[1]; 
      
      double hit_sigma_u = TMath::Sqrt(Hit.GetCov()(0,0)); 
      double hit_sigma_v = TMath::Sqrt(Hit.GetCov()(1,1)); 
      double hit_corr_uv = Hit.GetCov()(0,1)/(hit_sigma_u*hit_sigma_v); 
        
      histoName = "hsigma_hit_u_det"+planeNumber;
      _histoMap[ histoName  ]->Fill(hit_sigma_u);
      
      histoName = "hsigma_hit_v_det"+planeNumber;
      _histoMap[ histoName  ]->Fill(hit_sigma_v); 
      
      histoName = "hcorr_uv_det"+planeNumber;
      _histoMap[ histoName  ]->Fill(hit_corr_uv);  
      
      histoName = "hhit_u_det"+planeNumber;
      _histoMap[ histoName  ]->Fill(hit_u); 
      
      histoName = "hhit_v_det"+planeNumber;
      _histoMap[ histoName  ]->Fill(hit_v);
      
      histoName = "hhitmap"+planeNumber;
      _histoMap2D[ histoName  ]->Fill(hit_u,hit_v);   
      
      PixelCluster Cluster = Hit.GetCluster(); 
      
      histoName = "hdigitmap"+planeNumber;
      _histoMap2D[ histoName  ]->Fill(Cluster.getUSeed(),Cluster.getVSeed()); 
      
      histoName = "hcls_charge_sensor"+planeNumber;
      _histoMap[ histoName  ]->Fill(Cluster.getCharge()); 
      
      histoName = "hseed_charge_sensor"+planeNumber;
      _histoMap[ histoName  ]->Fill(Cluster.getSeedCharge()); 
      
      histoName = "hcls_type_sensor"+planeNumber;
      _histoMap[ histoName  ]->Fill(0); 
       
      histoName = "hsize_sensor"+planeNumber;
      _histoMap[ histoName  ]->Fill(Cluster.getSize()); 
      
      histoName = "hsizeU_sensor"+planeNumber;
      _histoMap[ histoName  ]->Fill(Cluster.getUSize()); 
      
      histoName = "hsizeV_sensor"+planeNumber;
      _histoMap[ histoName  ]->Fill(Cluster.getVSize()); 
      
    }
  }

}


//
// Method called after each event to check the data processed
//
void RawHitDQM::check( LCEvent * )
{
}

//
// Method called after all data processing
//
void RawHitDQM::end()
{
   
   // CPU time end
   _timeCPU = clock()/1000 - _timeCPU;

   for (int ipl = 0; ipl < TBDetector::GetInstance().GetNSensors(); ++ipl ) {
     
     string histoName = "hnhits_sensor"+to_string( ipl );
     double nhits = _histoMap[ histoName ]->GetMean();
     
     histoName = "hnhits";
     _histoMap[ histoName ]->SetBinContent(ipl+1, nhits);  
   }
   
   _rootFile->Write();
   _rootFile->Close();
   
   delete _rootFile;

}


//
// Method printing processor parameters
//
void RawHitDQM::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "RawHitDQM Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}


void RawHitDQM::bookHistos() {
  
  
  // BOOK HISTOS -------------------------- 
   
  // Get number of detectors
  //------------------------
  int nDet = TBDetector::GetInstance().GetNSensors();
  
  // Write ROOT Output
  _rootFile = new TFile(_rootFileName.c_str(),"recreate");
   
  
  std::string dirName; 
  std::string histoName;
  std::string histoTitle;
  
  histoName = "hnhits";
  _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", nDet, 0, nDet);
  _histoMap[ histoName ]->SetXTitle("plane number"); 
  _histoMap[ histoName ]->SetYTitle("mean hits per event");   
  
  // Create subdirs for detectors
  
  for (int ipl=0 ; ipl < nDet; ipl++) {
    std::string dirName = "Det"+to_string( ipl );
    _rootFile->mkdir(dirName.c_str());     
  }      
  
  // Detector histograms 
  for (int ipl=0 ; ipl < nDet; ipl++) {
    
    dirName = "/Det"+to_string(ipl)+"/";
    _rootFile->cd(dirName.c_str());
    
    
    double max; 
    double safetyFactor = 1.1;
          
    // This is a alignable detector   
    const Det & adet = TBDetector::Get(ipl); 
                
    histoName = "hsigma_hit_u_det"+to_string( ipl );
    histoTitle ="Cluster sigma u"; 
    max = 10*adet.GetSensitiveMaxU()/(adet.GetMaxUCell()+2); 
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 400, 0, max); 
    _histoMap[ histoName  ]->SetXTitle("cluster sigma u [mm]"); 
    _histoMap[ histoName  ]->SetYTitle("number of cluster");
    
    histoName = "hsigma_hit_v_det"+to_string( ipl );
    histoTitle ="Cluster sigma v"; 
    max = 10*adet.GetSensitiveMaxV()/(adet.GetMaxVCell()+2); 
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 400, 0, max);
    _histoMap[ histoName  ]->SetXTitle("cluster sigma v [mm]"); 
    _histoMap[ histoName  ]->SetYTitle("number of cluster"); 

    histoName = "hcorr_uv_det"+to_string( ipl );
    histoTitle ="Cluster correlation uv"; 
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 100, -1, +1);
    _histoMap[ histoName  ]->SetXTitle("uv correlation coefficient"); 
    _histoMap[ histoName  ]->SetYTitle("number of cluster"); 
   
    int minUCell = adet.GetMinUCell();
    int maxUCell = adet.GetMaxUCell();
    int minVCell = adet.GetMinVCell();
    int maxVCell = adet.GetMaxVCell();
    double  uMin = safetyFactor * adet.GetSensitiveMinU();
    double  uMax = safetyFactor * adet.GetSensitiveMaxU();
    int uBins = maxUCell-minUCell+1;             
    
    histoName = "hhit_u_det"+to_string( ipl );
    histoTitle ="Cluster u"; 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), histoTitle.c_str(), uBins ,  uMin, uMax ); 
    _histoMap[ histoName ]->SetXTitle("cluster u [mm]"); 
    _histoMap[ histoName ]->SetYTitle("number of clusters"); 

    double  vMin = safetyFactor * adet.GetSensitiveMinV();
    double  vMax = safetyFactor * adet.GetSensitiveMaxV();
    int vBins = maxVCell-minVCell+1;
    
    histoName = "hhit_v_det"+to_string( ipl );
    histoTitle ="Cluster v"; 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), histoTitle.c_str(), vBins, vMin, vMax);
    _histoMap[ histoName ]->SetXTitle("cluster v [mm]"); 
    _histoMap[ histoName ]->SetYTitle("number of clusters"); 
    
    histoName = "hhitmap"+to_string( ipl );
    histoTitle ="Hitmap for plane " +to_string( ipl )+" SensorID " + to_string( adet.GetSensorID()  );
    _histoMap2D[ histoName] = new TH2D(histoName.c_str(), histoTitle.c_str(),uBins, uMin, uMax, vBins, vMin, vMax);
    _histoMap2D[histoName]->SetXTitle("cluster u [mm]"); 
    _histoMap2D[histoName]->SetYTitle("cluster v [mm]");    
    _histoMap2D[histoName]->SetStats( false );  
    
    histoName = "hdigitmap"+to_string( ipl );
    histoTitle ="Hitmap for plane " +to_string( ipl )+" SensorID " + to_string( adet.GetSensorID() );
    _histoMap2D[ histoName] = new TH2D(histoName.c_str(), histoTitle.c_str(),uBins, minUCell, maxUCell+1, vBins, minVCell, maxVCell+1);// +1 for max with same bincount before, max is not included upper bin edge
    _histoMap2D[histoName]->SetXTitle("cluster u seed [cellID]"); 
    _histoMap2D[histoName]->SetYTitle("cluster v seed [cellID]");    
    _histoMap2D[histoName]->SetStats( false ); 
   
    histoName = "hcls_charge_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 255, 0, 255);
    _histoMap[ histoName ]->SetXTitle(" cluster charge [ADU]"); 
    _histoMap[ histoName ]->SetYTitle(" number of clusters");  

    histoName = "hseed_charge_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 255, 0, 255);
    _histoMap[ histoName ]->SetXTitle(" seed charge [ADU]"); 
    _histoMap[ histoName ]->SetYTitle(" number of clusters");      

    histoName = "hcls_type_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 10, 0, 10);
    _histoMap[ histoName ]->SetXTitle(" cluster type"); 
    _histoMap[ histoName ]->SetYTitle(" number of clusters");   

    histoName = "hsize_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 10, 0, 10);
    _histoMap[ histoName ]->SetXTitle(" cluster size [pixels]"); 
    _histoMap[ histoName ]->SetYTitle(" number of clusters");    
    
    histoName = "hsizeU_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 6, 0, 6);
    _histoMap[ histoName ]->SetXTitle(" cluster size u [cells]"); 
    _histoMap[ histoName ]->SetYTitle(" number of clusters");    
    
    histoName = "hsizeV_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 6, 0, 6);
    _histoMap[ histoName ]->SetXTitle(" cluster size v [cells]"); 
    _histoMap[ histoName ]->SetYTitle(" number of clusters");  
          
    histoName = "hnhits_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, 0, 200);
    _histoMap[ histoName ]->SetXTitle("hits per event"); 
    _histoMap[ histoName ]->SetYTitle("number of events");   
     
  } 
        
}


} // Namespace

