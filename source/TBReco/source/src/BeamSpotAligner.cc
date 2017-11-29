// BeamSpotAligner    
//
// See BeamSpotAligner.h for full documentation of this processor. 
//		
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "BeamSpotAligner.h"

// Include TBTools header files
#include "TBHit.h"
#include "TBTrack.h"
#include "TrackInputProvider.h"

// Include basic C
#include <iostream>
#include <limits>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/LCCollectionVec.h>

// ROOT includes
#include "TMath.h"
#include "TF2.h"
#include "TVirtualFitter.h"


// Used namespaces
using namespace std; 
using namespace lcio;
using namespace CLHEP; 
using namespace marlin;

namespace depfet {

//
// Instantiate this object
//
BeamSpotAligner aBeamSpotAligner ;

//
// Constructor
//
BeamSpotAligner::BeamSpotAligner() : Processor("BeamSpotAligner")
{ 

// Processor description
   _description = "BeamSpotAligner: Alignment of sensor XY position relative to the particle beam";
     
//   
// First of all, we need to register the triplett tracks 
    
   registerInputCollection(LCIO::TRACK,"TrackCollectionName",
                          "Tracklett collection for alignment",
                          _trackCollectionName,std::string("tracks"));
   
//   
// Define processor parameters 
   
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
                              _rootFileName, string("BeamSpot.root"));
    
}

//
// Method called at the beginning of data processing
//
void BeamSpotAligner::init() {
  
  // Initialize variables
  _nRun = 0 ;
  _nEvt = 0 ;
   
  // Read detector constants from gear file
  _detector.ReadGearConfiguration();    
  
  // Read alignment data base file 
  if(!_newAlignment) _detector.ReadAlignmentDB( _alignmentDBFileName );
  // This is needed, because if the AlignmentDB is not read, the detector construct doesn't know the alignmentDB name
  else  _detector.SetAlignmentDBName( _alignmentDBFileName );   
  
  // Book all needed histograms 
  bookHistos();
  
}

//
// Method called for each run
//
void BeamSpotAligner::processRunHeader(LCRunHeader * run)
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
void BeamSpotAligner::processEvent(LCEvent * evt)
{
  
  ++_nEvt;
     
  // Print event number
  if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                 << (evt->getEventNumber())
                                                                 << std::endl << std::endl;
   
  
  // Read input track collections
  // ===========================
   
  LCCollection* trackcol;
  try {
      trackcol = evt->getCollection(_trackCollectionName);
  } catch (DataNotAvailableException& e) {
      throw SkipEventException(this);
  }
  
  // Main loop over all tracks
   
  int nTracks = trackcol->getNumberOfElements(); 
  TrackInputProvider TrackLCIOReader;  
  
  for (int itrk = 0; itrk < nTracks; itrk++) {
    
    // Retrieve track from LCIO 
    Track * lciotrk = dynamic_cast<Track*> (trackcol->getElementAt(itrk));
    
    // Convert LCIO -> TB track  
    TBTrack rectrack = TrackLCIOReader.MakeTBTrack( lciotrk, _detector );  

    for(int ipl=0;ipl<_detector.GetNSensors();++ipl)  { 
      
      if ( rectrack.GetTE(ipl).HasHit() ) {
      
        double um = rectrack.GetTE(ipl).GetHit().GetCoord()[0][0];
        double vm = rectrack.GetTE(ipl).GetHit().GetCoord()[1][0];
        
        string histoName = "hhitmap_sensor"+to_string( ipl );
        _histoMap2D[ histoName  ]->Fill(um,vm); 
      }
    }
                                     
  }  
                     
}


//
// Method called after each event to check the data processed
//
void BeamSpotAligner::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void BeamSpotAligner::end()
{
  
  ////////////////////////////////////////////////////////////
  // Try to create a better aligned detector  
  
  TBDetector tmp_detector = _detector;
  
  streamlog_out(MESSAGE3) << "Calculating alignment corrections  ... Beam spot method" << endl << endl;   
    
  for ( int ipl = 0 ; ipl < tmp_detector.GetNSensors() ; ++ipl ) 
  {           
       
    if (ipl != 0) continue; 

    string histoName = "hhitmap_sensor"+to_string( ipl );
  
    // First the center of the beam spot is fitted using the local hit 
    // map. We assume a 2d Gaussian beam density. 

    double sizeU = _detector.GetDet(ipl).GetSensitiveSizeU();  
    double sizeV = _detector.GetDet(ipl).GetSensitiveSizeV(); 
    
    // Get range for beam spot fitting 
    double MinU = -0.5*sizeU;  
    double MaxU = +0.5*sizeU;  
    double MinV = -0.5*sizeV;  
    double MaxV = +0.5*sizeV;  
     
    // Define the fit function with ranges and initial parameters 
    const int npar = 5;  
    double f2params[npar] = {100,0,6,0,6,};
    TF2 *f2 = new TF2("f2",beamspotfitter::gauss2Dsimple,MinU,MaxU,MinV,MaxV, npar);
    f2->SetParameters(f2params);
    f2->SetParNames("Constant","MeanU","SigmaU","MeanV","SigmaV");    
     
    //Fit beam spot histogram 
    _histoMap2D[ histoName  ]->Fit("f2","R");
    f2->Draw("cont1 same");
      
    // Read fit results 
    double offsetU = f2->GetParameter(1); 
    double offsetV = f2->GetParameter(3);
    
    // Discard fit result, if beam spot outside sensor area
    if ( ! _detector.GetDet(ipl).SensitiveCrossed(offsetU, offsetV)) continue;
    
    streamlog_out(MESSAGE3) << "Plane number : " << ipl << endl  
                            << "   U offset is " << offsetU << " mm"
                            << endl
                            << "   V offset is " << offsetV << " mm"
                            << endl << endl;  
     
    // Alignment shifts in x/y are computed to move the fitted center on 
    // the z axis. 
    
    // We have calculated offset in local coord  
    HepVector local_offset(3);
    local_offset[0] = -offsetU; 
    local_offset[1] = -offsetV;      
    local_offset[2] = 0;  
      
    // Transform local offset to global offset   
    HepVector global_offset = _detector.GetDet(ipl).GetNominal().GetRotation().T() * local_offset ;  
    double dx = global_offset[0]; 
    double dy = global_offset[1];      
    double dz = global_offset[2];   
             
    // Compute a 'delta' frame from corrections - no tilt corrections  
    ReferenceFrame deltaFrame = ReferenceFrame::create_karimaki_delta(dx,dy,dz,0,0,0); 
      
    // Merge nominal frame and delta frame 
    ReferenceFrame alignFrame = ReferenceFrame::combine_karimaki(_detector.GetDet(ipl).GetNominal(), deltaFrame); 
      
    // Update nominal sensor reference frame
    _detector.GetDet(ipl).SetNominalFrame(alignFrame); 
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
    streamlog_out ( MESSAGE3 ) << "NO UPDATE OF ALIGNMENT DB LCIO FILE" << endl; 
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

void BeamSpotAligner::bookHistos() {


  // ROOT Output 
  _rootFile = new TFile(_rootFileName.c_str(),"recreate");

  string histoName;

  for(int ipl=0; ipl < _detector.GetNSensors(); ++ipl)    
  {
    
    // Get handle to sensor data
    Det & Sensor = _detector.GetDet(ipl); 
    
    double  uBox = 1.0 * 0.5 * Sensor.GetSensitiveSizeU();   
    int uBins = Sensor.GetNColumns()/12;    
    
    double  vBox = 1.0 * 0.5 * Sensor.GetSensitiveSizeV();
    int vBins = Sensor.GetNRows()/12; 

    histoName = "hhitmap_sensor"+to_string( ipl );
    _histoMap2D[ histoName] = new TH2D(histoName.c_str(), "" ,uBins, -uBox, +uBox, vBins, -vBox, +vBox);
    _histoMap2D[histoName]->SetXTitle("u [mm]"); 
    _histoMap2D[histoName]->SetYTitle("v [mm]");    
    _histoMap2D[histoName]->SetStats( false );  
  }

}


} // depfet Namespace


/** Bivariate gaussian fit model in sensor u,v coordinates 
 *  Summary of model parameters: 
 *
 * par[0] :	amplitude factor 
 * par[1] : 	mean in u in mm
 * par[2] :  	standard deviation in u in mm
 * par[3] : 	mean in v in mm
 * par[4] :   	standard deviation in v in mm
 * par[5] : 	correlation coefficient
 *
 * Summary of variables 
 *
 * x[0]   :   	u coordinate in mm
 * x[1]   :	v coordinate in mm
 */ 
double beamspotfitter::gauss2D(double *x, double *par) {
   
   double z1 = double((x[0]-par[1])/par[2]);
   double z2 = double((x[1]-par[3])/par[4]);
   double r  = double(par[5]); 
   return par[0]*exp(-(z1*z1+z2*z2-2*r*z1*z2)/(2*(1-r*r)));
}   

/** Bivariate gaussian fit model in sensor u,v coordinates.
 *  Simplified model without correlations.  
 *  Summary of model parameters: 
 *
 * par[0] :	amplitude factor 
 * par[1] : 	mean in u in mm
 * par[2] :  	standard deviation in u in mm
 * par[3] : 	mean in v in mm
 * par[4] :   	standard deviation in v in mm
 *
 * Summary of variables 
 *
 * x[0]   :   	u coordinate in mm
 * x[1]   :	v coordinate in mm
 */ 
double beamspotfitter::gauss2Dsimple(double *x, double *par) {
   
   double r1 = double((x[0]-par[1])/par[2]);
   double r2 = double((x[1]-par[3])/par[4]);
   return par[0]*exp(-0.5*(r1*r1+r2*r2));
}   



