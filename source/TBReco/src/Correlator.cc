// Correlator processor   
//
// See Correlator.h for full documentation of this processor. 
//		
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "Correlator.h"

// Include TBTools header files
#include "TBDetector.h"
#include "TBHit.h"
#include "TBTrack.h"
#include "Det.h"
#include "Utilities.h"
#include "HitFactory.h"
#include "SeedGenerator.h"
#include "GenericTrackFitter.h"

// Include basic C
#include <cstdlib>
#include <iostream>
#include <limits>
#include <algorithm>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>

#include <Eigen/Core>

// ROOT includes
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>

using Eigen::Vector3d;

// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;

namespace depfet {

//
// Instantiate this object
//
Correlator aCorrelator ;

//
// Constructor
//
Correlator::Correlator() : Processor("Correlator")
{ 

// Processor description
   _description = "Correlator: Pre- alignment in xy of tracking telescopes";
     
//   
// First of all, we need to register the input/output collections
    
   vector< string > inputHitCollectionNameVecExample;
   inputHitCollectionNameVecExample.push_back( "hit" );
   
   registerInputCollections (LCIO::TRACKERHIT, "InputHitCollectionNameVec",
                            "Hit collection names",
                            _inputHitCollectionNameVec, inputHitCollectionNameVecExample );
   
//   
// Define processor parameters 
  
   registerProcessorParameter ("OutputRootFileName",
                              "This is the name of the output root file",
                              _rootFileName, string("XCorrelator.root"));
   
   registerProcessorParameter ("ReferencePlane",
                              "Reference sensor plane number, counted along beam line",
                              _refPlane, static_cast <int> (0) );
   
   registerProcessorParameter ("ParticleMomentum", "Particle momentum [GeV]",
                              _momentum,  static_cast < double > (4.0));
   
   registerProcessorParameter ("ParticleMass", "Particle mass [GeV]",
                              _mass,  static_cast < double > (0.000511));
   
   registerProcessorParameter ("ParticleCharge", "Particle charge [e]",
                              _charge,  static_cast < double > (-1));
     
}

//
// Method called at the beginning of data processing
//
void Correlator::init() {
  
  // Initialize variables
  _iRun = 0 ;
  _iEvt = 0 ;
   
  // Initialize new random generator with unique seed 
  gRandom = new TRandom3();
  gRandom->SetSeed(0);

  // Print random generator info
  streamlog_out(MESSAGE3) << "Random Generator setup with seed: "
                          << (gRandom->GetSeed())
                          << std::endl << std::endl;
  
  // Book correlation histograms   
  bookHistos();   
  
}

//
// Method called for each run
//
void Correlator::processRunHeader(LCRunHeader * run)
{
     
  // Print run number
  streamlog_out(MESSAGE3) << "Processing run: "
                          << (run->getRunNumber())
                          << std::endl << std::endl;
  _iRun++ ;

 
 
}

//
// Method called for each event
//
void Correlator::processEvent(LCEvent * evt)
{
     
  // Print event number
  if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                 << (evt->getEventNumber())
                                                                 << std::endl << std::endl;
   
  ++_iEvt;
          
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
         
         // Only use good hits 
         if ( RecoHit.GetQuality() == 0  ) { 
           HitStore.AddRecoHit(RecoHit);
         }
                      
       } // End for hit loop      
          
    } catch (lcio::DataNotAvailableException& e) {
       streamlog_out ( MESSAGE2 ) << "Not able to get collection "
                                  << _inputHitCollectionNameVec.at( iCol )
                                  << "\nfrom event " << evt->getEventNumber()
                                  << " in run " << evt->getRunNumber()  << endl;
       
    }  
  }   
  
  streamlog_out ( MESSAGE2 ) << "Total of " << HitStore.GetNHits() << " good hits." << endl;
      
  // Fill correlation histos 
  // ============================
  // 
  // Beam constrained tracks (BT) are formed from each hit 
  // on the reference detector assuming beam parallel to 
  // telescope z axis.
  //  
  // Then, BT's are extrapolated to local sensors using the 
  // nominal alignment constants. Hit pairs (HP) are formed 
  // by merging a BT with a local hit. For each HP, we record
  // BT impact position and HP residuals. 
       
  // Configure track fitter 
  GenericTrackFitter TrackFitter(TBDetector::GetInstance());
  TrackFitter.SetNumIterations(1);   
  
  // Configue seed track generator  
  SeedGenerator TrackSeeder(_charge, _momentum);
  
  for (int iref = 0; iref < HitStore.GetNHits(_refPlane); ++iref ) {
      
    //    
    // Create a beam constrained track 
      
    TBTrack rectrack(TBDetector::GetInstance());
    
    // Set particle hypothesis 
    rectrack.SetMass( _mass );
    rectrack.SetCharge( _charge );
    rectrack.SetMomentum( _momentum ); 
       
    // Set seed paramaters   
    const TBHit& refhit = HitStore.GetRecoHitFromID(iref, _refPlane);  
    TBTrackState Seed = TrackSeeder.CreateSeedTrack(refhit, TBDetector::GetInstance());   
    rectrack.SetReferenceState(Seed);
    
    // Extrapolate seed to all planes
    bool trkerr = TrackFitter.ExtrapolateSeed(rectrack);
    if ( trkerr ) { // just skip the track 
      continue;
    }  
    
    //                        
    // Loop over detector planes 
     
    for(int ipl=0;ipl<TBDetector::GetInstance().GetNSensors();++ipl)  { 
       
      // Check iff sensor is  reference sensor
      if ( ipl == _refPlane ) continue;    
      
      // Skip track if it does not intersect
      if (!rectrack.GetTE(ipl).IsCrossed()) continue;

      // Get extrapolated intersection coordinates
      double u = rectrack.GetTE(ipl).GetState().GetPars()[2]; 
      double v = rectrack.GetTE(ipl).GetState().GetPars()[3];  

      const Det & det = TBDetector::GetInstance().GetDet(ipl); 

	  double  PitchU = ( det.GetSensitiveMaxU() - det.GetSensitiveMinU() ) /(det.GetMaxUCell()-det.GetMinUCell()+1); 
	  double  PitchV = ( det.GetSensitiveMaxV() - det.GetSensitiveMinV() ) /(det.GetMaxVCell()-det.GetMinVCell()+1);     
           
      // Loop over all hits on this detector 
      for (int ihit = 0; ihit < HitStore.GetNHits(ipl); ++ihit ) 
      {      
        
        const TBHit & anyhit = HitStore.GetRecoHitFromID(ihit, ipl);
           
        // Measured hit coordinates
        Vector3d anypos = anyhit.GetLocalSpacePoint();
        double um = anypos[0]; 
        double vm = anypos[1];  

		// Smear out residuals with a Gaussian to remove peaks from residual distribution
		double random_U = gRandom->Gaus(0, TMath::Sqrt( PitchU/2.0 ));
		double random_V = gRandom->Gaus(0, TMath::Sqrt( PitchU/2.0 ));
        
        _hitUCorrelationMatrix[ ipl ] -> Fill ( u, um + random_U) ;
        _hitVCorrelationMatrix[ ipl ] -> Fill ( v, vm + random_V) ;
                     
        _hitUResidualHisto[ ipl ]->Fill( um - u + random_U);
        _hitVResidualHisto[ ipl ]->Fill( vm - v + random_V);
              
      }  // End inner hit loop 
           
    } // End plane loop 
                  
  } // End outer hit loop 
                   
}


//
// Method called after each event to check the data processed
//
void Correlator::check( LCEvent * )
{
}

//
// Method called after all data processing
//
void Correlator::end()
{
  
  ////////////////////////////////////////////////////////////
  // Try to create a better aligned detector  
   
  streamlog_out(MESSAGE3) << "Calculating alignment corrections  ... Correlation band method" << endl << endl;        
  

  for ( int ipl = 0 ; ipl < TBDetector::GetInstance().GetNSensors() ; ++ipl ) 
  {           
       
    if (ipl == _refPlane) continue; 
       
    
    // Find position of residual peak in U
    
    float peakFrequencyU = 0.;
    float peakPositionU = 0.;
    for(int ibin = 1; ibin <= _hitUResidualHisto[ ipl ]->GetXaxis()->GetNbins(); ibin++)
    {
      double shift = _hitUResidualHisto[ ipl ]->GetXaxis()->GetBinCenter(ibin);
      double frequency = _hitUResidualHisto[ ipl ]->GetBinContent( ibin ); 
                  
      if( frequency > peakFrequencyU )
      {
        peakPositionU = shift; 
        peakFrequencyU = frequency;
      }
    }
        
    // Find position of residual peak in U       
               
    float peakFrequencyV = 0.;
    float peakPositionV = 0.;
    for(int ibin = 1; ibin <= _hitVResidualHisto[ ipl ]->GetXaxis()->GetNbins(); ibin++)
    {
      double shift = _hitVResidualHisto[ ipl ]->GetXaxis()->GetBinCenter(ibin);
      double frequency = _hitVResidualHisto[ ipl ]->GetBinContent( ibin ); 
                  
      if( frequency > peakFrequencyV )
      {
        peakPositionV = shift; 
        peakFrequencyV = frequency;
      }
    }
    
    // Use alignment offsets to shift the residual peaks to zero     
    double offsetU = peakPositionU; 
    double offsetV = peakPositionV; 
    
    // Compute alignment corrections  
             
    streamlog_out(MESSAGE3) << "Plane number : " << ipl << endl  
                            << "   U offset is " << offsetU << " mm"
                            << endl
                            << "   V offset is " << offsetV << " mm"
                            << endl << endl;  
      
    
    Det & adet = TBDetector::GetInstance().GetDet(ipl);
      
    // We have calculated offset in local coord.
    Vector3d local_offset;
    local_offset << -offsetU, -offsetV, 0;
      
    // And transform to offsets in global coord.  
    Vector3d global_offset = adet.GetNominal().GetRotation().transpose() * local_offset ;
    double dx = global_offset[0]; 
    double dy = global_offset[1];      
    double dz = global_offset[2];    
          
    // Compute a 'delta' frame from corrections - no tilt corrections  
    ReferenceFrame deltaFrame = ReferenceFrame::create_karimaki_delta(dx,dy,dz,0,0,0); 
      
    // Merge nominal frame and delta frame 
    ReferenceFrame alignFrame = ReferenceFrame::combine_karimaki(adet.GetNominal(), deltaFrame); 
      
    // Update nominal sensor reference frame
    adet.SetNominalFrame(alignFrame); 
  
  }
       
  streamlog_out ( MESSAGE3 ) << "Detector geometry after alignment procedure: " << endl;   
  TBDetector::GetInstance().Print();   
       
  // Print message
  streamlog_out(MESSAGE3) << std::endl
                          << "Processor succesfully finished!"
                          << std::endl;
   
  // ROOT Output 
  _rootFile->Write();   
  _rootFile->Close(); 
  
  delete _rootFile; 
    
}


void Correlator::bookHistos() {
  
  streamlog_out ( MESSAGE4 ) <<  "Booking histograms" << endl;
  
  // ROOT Output 
  _rootFile = new TFile(_rootFileName.c_str(),"recreate");
     
  // Create all the directories first
  vector< string > dirNames;
  
  // This is the reference detector   
  const Det & refdet = TBDetector::GetInstance().GetDet(_refPlane); 
    
  dirNames.push_back ("HitU");
  dirNames.push_back ("HitV");
  dirNames.push_back ("HitUShift");
  dirNames.push_back ("HitVShift");
     
  for ( size_t iPos = 0 ; iPos < dirNames.size() ; iPos++ ) {
    _rootFile->mkdir( dirNames[iPos].c_str()  ); 
  }
    
  string tempHistoName;
  string tempHistoTitle;
  string tempXAxis;
  string tempYAxis; 
  
  for(int ipl=0; ipl < TBDetector::GetInstance().GetNSensors(); ++ipl)    
  {
        
    if ( ipl == _refPlane ) continue; 
                
    // This is a alignable detector   
    const Det & det = TBDetector::GetInstance().GetDet(ipl); 
    
    double safetyFactor = 2.0;  // 2 should be enough because it
                                // means that the sensor is wrong
                                // by all its size.
        
    
    /////////////////////////////////////////////////
    // book U histos
              
    _rootFile->cd("/HitU/");
         
    double  uMin = safetyFactor  * det.GetSensitiveMinU();               
    double  uMax = + safetyFactor * det.GetSensitiveMaxU(); 
    double  PitchU = (det.GetSensitiveMaxU() - det.GetSensitiveMinU() )/(det.GetMaxUCell()-det.GetMinUCell()+1); 
    int     uBins = static_cast<int>( (uMax - uMin)/(2.01*PitchU) );     
    
    // avoid too many bins 
    if (uBins > 2000 ) uBins = 2000;            
    
    double  uMinRef = safetyFactor * refdet.GetSensitiveMinU();               
    double  uMaxRef = safetyFactor * refdet.GetSensitiveMaxU(); 
    double  refPitchU = (refdet.GetSensitiveMaxU()  - refdet.GetSensitiveMinU() )/(refdet.GetMaxUCell()-refdet.GetMinUCell()+1); 
    int     uBinsRef = static_cast<int>( (uMaxRef - uMinRef)/(2.01*refPitchU) );    
    
    // avoid too many bins 
    if (uBinsRef > 2000 ) uBinsRef = 2000;   

           
    tempHistoName  = "HitUCorrelationHisto_d" + to_string( ipl );
    tempHistoTitle = "HitUCorrelationHisto_d" + to_string( ipl );

    _hitUCorrelationMatrix[ ipl  ] = new TH2D(tempHistoName.c_str(),tempHistoTitle.c_str(),uBinsRef, uMinRef, uMaxRef, uBins, uMin, uMax);
               
    tempXAxis = "Detector " + to_string( _refPlane ) + " Hit U [mm]";
    tempYAxis = "Detector " + to_string( ipl ) + " Hit U [mm]";
    _hitUCorrelationMatrix[ ipl  ]->SetXTitle(tempXAxis.c_str()); 
    _hitUCorrelationMatrix[ ipl  ]->SetYTitle(tempYAxis.c_str());   
    _hitUCorrelationMatrix[ ipl  ]->SetStats( false );    
     
     

    /////////////////////////////////////////////////
    // book V
                
    _rootFile->cd("/HitV/");
     
    double  vMin = safetyFactor * det.GetSensitiveMinV();               
    double  vMax = safetyFactor * det.GetSensitiveMaxV(); 
    double  PitchV = ( det.GetSensitiveMaxV() - det.GetSensitiveMinV() ) /(det.GetMaxVCell()-det.GetMinVCell()+1); 
    int     vBins = static_cast<int>( (vMax - vMin)/(2.01*PitchV) );     
        
    // avoid too many bins 
    if (vBins > 2000 ) vBins = 2000;   

    double  vMinRef = safetyFactor * refdet.GetSensitiveMinV();               
    double  vMaxRef = safetyFactor * refdet.GetSensitiveMaxV(); 
    double  refPitchV = ( refdet.GetSensitiveMaxV()- refdet.GetSensitiveMinV() )/(refdet.GetMaxVCell()-refdet.GetMinVCell()+1); 
    int     vBinsRef = static_cast<int>( (vMaxRef - vMinRef)/(2.01*refPitchV) );     
        
    // avoid too many bins 
    if (vBinsRef > 2000 ) vBinsRef = 2000;   
    
    tempHistoName =  "HitVCorrelationHisto_d" + to_string( ipl );
    tempHistoTitle = "HitVCorrelationHisto_d" + to_string( ipl );
    
    _hitVCorrelationMatrix[ ipl  ] = new TH2D(tempHistoName.c_str(),tempHistoTitle.c_str(),vBinsRef, vMinRef, vMaxRef, vBins, vMin, vMax);
             
    tempXAxis = "Detector " + to_string( _refPlane ) + " Hit V [mm]";
    tempYAxis = "Detector " + to_string( ipl ) + " Hit V [mm]";
    _hitVCorrelationMatrix[ ipl  ]->SetXTitle(tempXAxis.c_str()); 
    _hitVCorrelationMatrix[ ipl  ]->SetYTitle(tempYAxis.c_str());    
    _hitVCorrelationMatrix[ ipl  ]->SetStats( false );  
        
    /////////////////////////////////////////////////
    // book U residuals 
            
    _rootFile->cd("/HitUShift/");
     
    tempHistoName =  "HitUResidualsHisto_d" + to_string( ipl );
    tempHistoTitle =  "HitUResidualsHisto_d" + to_string( ipl );
               
    _hitUResidualHisto[ ipl ] = new TH1D(tempHistoName.c_str(), tempHistoTitle.c_str(), std::max(uBinsRef,uBins), std::min(uMinRef,uMin), std::max(uMaxRef,uMax));
              
    tempXAxis = "Detector " + to_string( ipl ) + " Residual U [mm]";
    tempYAxis = "# Hits";
    _hitUResidualHisto[ ipl ]->SetXTitle(tempXAxis.c_str()); 
    _hitUResidualHisto[ ipl ]->SetYTitle(tempYAxis.c_str());         
                         
    /////////////////////////////////////////////////
    // book V residuals 
    
    _rootFile->cd("/HitVShift/");
      
    tempHistoName = "HitVResidualsHisto_d" + to_string( ipl ) ;
    tempHistoTitle = "HitVResidualsHisto_d" + to_string( ipl ) ;
                          
    _hitVResidualHisto[ ipl ] = new TH1D(tempHistoName.c_str(), tempHistoTitle.c_str(), std::max(vBinsRef,vBins), std::min(vMinRef,vMin), std::max(vMaxRef,vMax));
             
    tempXAxis = "Detector " + to_string( ipl ) + " Residual V [mm]";
    tempYAxis = "# Hits";
    _hitVResidualHisto[ ipl ]->SetXTitle(tempXAxis.c_str()); 
    _hitVResidualHisto[ ipl ]->SetYTitle(tempYAxis.c_str());            
                    
                   
  } 
        
}


} // Namespace

