// TriplettCorrelator    
//
// See TriplettCorrelator.h for full documentation of this processor. 
//		
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "TriplettCorrelator.h"

// Include DEPFETTrackTools header files
#include "TBHit.h"
#include "TBTrack.h"
#include "Det.h"
#include "Utilities.h"
#include "HitFactory.h"
#include "GenericTrackFitter.h"
#include "TrackInputProvider.h"

// Include basic C
#include <cstdlib>
#include <iostream>
#include <limits>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>

// Include CLHEP classes
#include <CLHEP/Matrix/Vector.h>

// Used namespaces
using namespace std; 
using namespace lcio;
using namespace CLHEP; 
using namespace marlin;

namespace depfet {

//
// Instantiate this object
//
TriplettCorrelator aTriplettCorrelator ;

//
// Constructor
//
TriplettCorrelator::TriplettCorrelator() : Processor("TriplettCorrelator")
{ 

// Processor description
   _description = "TriplettCorrelator: Pre- alignment in xy of tracking telescopes with trackletts";
     
//   
// Define processor parameters 
    
   vector< string > inputHitCollectionNameVecExample;
   inputHitCollectionNameVecExample.push_back( "hit" );
   
   registerInputCollections (LCIO::TRACKERHIT, "InputHitCollectionNameVec",
                            "Hit collection names",
                            _inputHitCollectionNameVec, inputHitCollectionNameVecExample );

   registerInputCollection(LCIO::TRACK,"TrackCollectionName",
                          "Tracklett collection for alignment",
                          _trackCollectionName,std::string("aligntracks"));
      
   registerProcessorParameter ("AlignmentDBFileName",
                             "This is the name of the file with the alignment constants (add .root)",
                             _alignmentDBFileName, static_cast< string > ( "alignmentDB.root" ) );   
                   
   registerProcessorParameter ("UpdateAlignment",
                              "Update alignment DB using offset corrections (true/false)?",
                              _updateAlignment, static_cast <bool> (false) ); 
 
   registerProcessorParameter ("NewAlignment",
                              "Start alignment from scratch (true/false)?",
                              _newAlignment, static_cast <bool> (false) ); 
   
   registerProcessorParameter ("OutputRootFileName",
                              "This is the name of the output root file",
                              _rootFileName, string("XCorrelator.root"));
   
   
   
   
     
}

//
// Method called at the beginning of data processing
//
void TriplettCorrelator::init() {
  
  // Initialize variables
  _iRun = 0 ;
  _iEvt = 0 ;
   
  // Read detector constants from gear file
  _detector.ReadGearConfiguration();    
  
  // Read alignment data base file 
  if(!_newAlignment) _detector.ReadAlignmentDB( _alignmentDBFileName );
  // This is needed, because if the AlignmentDB is not read, the detector construct doesn't know the alignmentDB name
  else  _detector.SetAlignmentDBName( _alignmentDBFileName );
  
  // Apply pre-alignment only to active sensors
  _isActive.resize(_detector.GetNSensors(), true);   
   
  // Book correlation histograms   
  bookHistos();   
  
}

//
// Method called for each run
//
void TriplettCorrelator::processRunHeader(LCRunHeader * run)
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
void TriplettCorrelator::processEvent(LCEvent * evt)
{
     
  // Print event number
  if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                 << (evt->getEventNumber())
                                                                 << std::endl << std::endl;
   
  ++_iEvt;
          
  // Read input hit collections
  // ===========================
  
  // Hit factory sorts hits according to plane    
  HitFactory HitStore( _detector );   
    
  for ( size_t iCol = 0 ; iCol < _inputHitCollectionNameVec.size(); ++iCol ) {
     
    try {
       
       LCCollectionVec * hitcol = dynamic_cast < LCCollectionVec * > (evt->getCollection( _inputHitCollectionNameVec.at( iCol ) ) );
              
       for ( int ihit = 0 ; ihit < (int) hitcol->size() ; ++ihit ) {
             
         // Built a TBHit     
         TrackerHitImpl * lciohit = dynamic_cast< TrackerHitImpl* > ( hitcol->getElementAt( ihit ) );
         TBHit RecoHit ( lciohit ); 
         
         // Only use good hits 
         if ( RecoHit.GetQuality() == 0  ) { 
           // Data for correlation bands 
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
      
  LCCollection* trackcol;
  try {
      trackcol = evt->getCollection(_trackCollectionName);
  } catch (DataNotAvailableException& e) {
      throw SkipEventException(this);
  }


  // Fill correlation histos 
  // ============================
  // 
  // Loop over tracks and refit them. Then loop over  
  // all planes and fill residuals histos between the
  // extrapolated track and all detector hits. 
       
  // Main loop over all tracks
   
  int nTracks = trackcol->getNumberOfElements(); 
  TrackInputProvider TrackLCIOReader;  

  // Configure track fitter 
  GenericTrackFitter TrackFitter(_detector);
  TrackFitter.SetNumIterations(1);    

  for (int itrk = 0; itrk < nTracks; itrk++) {
    
    // Retrieve track from LCIO 
    Track * lciotrk = dynamic_cast<Track*> (trackcol->getElementAt(itrk));
    
    // Convert LCIO -> TB track  
    TBTrack rectrack = TrackLCIOReader.MakeTBTrack( lciotrk, _detector );  
 
    // Fit the track. This also extrapolates to all known detector planes
    bool trkerr = TrackFitter.Fit(rectrack);
    if ( trkerr ) { // just skip the track 
      continue;
    }  
    
    for(int ipl=0;ipl<_detector.GetNSensors();++ipl)  { 
       
      // Skip track if it does not intersect
      if (!rectrack.GetTE(ipl).IsCrossed()) continue;
      
      // Do not touch sensors contributing hits in trackletts 
      if (rectrack.GetTE(ipl).HasHit() ) _isActive[ipl]=false; 

      // Get extrapolated intersection coordinates
      double u = rectrack.GetTE(ipl).GetState().GetPars()[2][0]; 
      double v = rectrack.GetTE(ipl).GetState().GetPars()[3][0];           
           
      

      // Loop over all hits on this detector 
      for (int ihit = 0; ihit < HitStore.GetNHits(ipl); ++ihit ) 
      {      
        TBHit & anyhit = HitStore.GetRecoHitFromID(ihit, ipl);
           
        // Measured hit coordinates
        HepVector anypos = anyhit.GetLocalSpacePoint();
        double um = anypos[0]; 
        double vm = anypos[1];   
        
        _hitUCorrelationMatrix[ ipl ] -> Fill ( u, um ) ;
        _hitVCorrelationMatrix[ ipl ] -> Fill ( v, vm ) ;
                     
        _hitUResidualHisto[ ipl ]->Fill( um - u);
        _hitVResidualHisto[ ipl ]->Fill( vm - v);

        
              
      }  // End hit loop 
           
    } // End ipl loop   
                               
  } // End loop over all tracks 
                 
}


//
// Method called after each event to check the data processed
//
void TriplettCorrelator::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void TriplettCorrelator::end()
{
  
  ////////////////////////////////////////////////////////////
  // Try to create a better aligned detector  
  
  TBDetector tmp_detector = _detector;
   
  streamlog_out(MESSAGE3) << "Calculating alignment corrections  ... Correlation band method" << endl << endl;        
  

  for ( int ipl = 0 ; ipl < tmp_detector.GetNSensors() ; ++ipl ) 
  {           
       
    // Only correct position of active sensors
    if (!_isActive[ipl]) continue; 
       
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

    Det & adet = tmp_detector.GetDet(ipl);
      
    // We have calculated offset in local coord.
    HepVector local_offset(3);
    local_offset[0] = -offsetU; 
    local_offset[1] = -offsetV;      
    local_offset[2] = 0;  
      
    // And transform to offsets in global coord.  
    HepVector global_offset = adet.GetNominal().GetRotation().T() * local_offset ;  
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
      
  //////////////////////////////////////////////////////////////////////  
  // Save alignment DB
      
  _detector = tmp_detector; 
   
  if ( _updateAlignment ) { 
    
    ////////////////////////////////////////////////////////////
    // Print detector geometrie after alignment   
    
    streamlog_out ( MESSAGE3 ) << "Detector geometry after alignment procedure: " << endl;   
    _detector.Print();   
       
    ////////////////////////////////////////////////////////////
    // Overwrite the alignment data base with new geometry 
    _detector.WriteAlignmentDB( ); 
    
  } else {
    streamlog_out ( MESSAGE3 ) << endl;
    streamlog_out ( MESSAGE3 ) << "NO UPDATE OF ALIGNMENT DB" << endl; 
  }
  
  // Print message
  streamlog_out(MESSAGE3) << std::endl
                          << "Processor succesfully finished!"
                          << std::endl;
   
  // ROOT Output 
  _rootFile->Write();   
  _rootFile->Close(); 
  
  delete _rootFile; 
    
}


void TriplettCorrelator::bookHistos() {
  
  streamlog_out ( MESSAGE4 ) <<  "Booking histograms" << endl;
  
  // ROOT Output 
  _rootFile = new TFile(_rootFileName.c_str(),"recreate");
     
  // Create all the directories first
  vector< string > dirNames;
    
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
  
  
  for(int ipl=0; ipl < _detector.GetNSensors(); ++ipl)    
  {
                        
    // This is a alignable detector   
    Det & det = _detector.GetDet(ipl); 
             
    double safetyFactor = 1.4;  // 2 should be enough because it
                                // means that the sensor is wrong
                                // by all its size.
        
    /////////////////////////////////////////////////
    // book U histos
              
    _rootFile->cd("/HitU/");
         
    double  uMin = - safetyFactor * 0.5 * det.GetSensitiveSizeU();               
    double  uMax = + safetyFactor * 0.5 * det.GetSensitiveSizeU(); 
    int     uBins = static_cast<int>( (uMax - uMin)/(2*det.GetPitchU()) );     
    
    // avoid too many bins 
    if (uBins > 2000 ) uBins = 2000;            
               
    tempHistoName  = "HitUCorrelationHisto_d" + to_string( ipl );
    tempHistoTitle = "HitUCorrelationHisto_d" + to_string( ipl );

    _hitUCorrelationMatrix[ ipl  ] = new TH2D(tempHistoName.c_str(),tempHistoTitle.c_str(),uBins, uMin, uMax, uBins, uMin, uMax);
               
    tempXAxis = "Tracklett U [mm]";
    tempYAxis = "Detector " + to_string( ipl ) + " Hit U [mm]";
    _hitUCorrelationMatrix[ ipl  ]->SetXTitle(tempXAxis.c_str()); 
    _hitUCorrelationMatrix[ ipl  ]->SetYTitle(tempYAxis.c_str());   
    _hitUCorrelationMatrix[ ipl  ]->SetStats( false );  
         
    /////////////////////////////////////////////////
    // book V
                
    _rootFile->cd("/HitV/");
     
    double  vMin = - safetyFactor * 0.5 * det.GetSensitiveSizeV();               
    double  vMax = + safetyFactor * 0.5 * det.GetSensitiveSizeV(); 
    int     vBins = static_cast<int>( (vMax - vMin)/(2*det.GetPitchV()) );     
        
    // avoid too many bins 
    if (vBins > 2000 ) vBins = 2000;   
    
    tempHistoName =  "HitVCorrelationHisto_d" + to_string( ipl );
    tempHistoTitle = "HitVCorrelationHisto_d" + to_string( ipl );

    _hitVCorrelationMatrix[ ipl  ] = new TH2D(tempHistoName.c_str(),tempHistoTitle.c_str(),vBins, vMin, vMax, vBins, vMin, vMax);
             
    tempXAxis = "Tracklett V [mm]";
    tempYAxis = "Detector " + to_string( ipl ) + " Hit V [mm]";
    _hitVCorrelationMatrix[ ipl  ]->SetXTitle(tempXAxis.c_str()); 
    _hitVCorrelationMatrix[ ipl  ]->SetYTitle(tempYAxis.c_str());    
    _hitVCorrelationMatrix[ ipl  ]->SetStats( false );        
  

    /////////////////////////////////////////////////
    // book U residuals 
            
    _rootFile->cd("/HitUShift/");
     
    tempHistoName =  "HitUResidualsHisto_d" + to_string( ipl );
    tempHistoTitle =  "HitUResidualsHisto_d" + to_string( ipl );
               
                   
    _hitUResidualHisto[ ipl ] = new TH1D(tempHistoName.c_str(), tempHistoTitle.c_str(), uBins, uMin, uMax);
              
    tempXAxis = "Detector " + to_string( ipl ) + " Residual U [mm]";
    tempYAxis = "# Hits";
    _hitUResidualHisto[ ipl ]->SetXTitle(tempXAxis.c_str()); 
    _hitUResidualHisto[ ipl ]->SetYTitle(tempYAxis.c_str());         
                         
    /////////////////////////////////////////////////
    // book V residuals 
    
    _rootFile->cd("/HitVShift/");
      
    tempHistoName = "HitVResidualsHisto_d" + to_string( ipl ) ;
    tempHistoTitle = "HitVResidualsHisto_d" + to_string( ipl ) ;
                          
    _hitVResidualHisto[ ipl ] = new TH1D(tempHistoName.c_str(), tempHistoTitle.c_str(), vBins, vMin, vMax);
             
    tempXAxis = "Detector " + to_string( ipl ) + " Residual V [mm]";
    tempYAxis = "# Hits";
    _hitVResidualHisto[ ipl ]->SetXTitle(tempXAxis.c_str()); 
    _hitVResidualHisto[ ipl ]->SetYTitle(tempYAxis.c_str());                    
                   
  } 
        
}


} // Namespace

