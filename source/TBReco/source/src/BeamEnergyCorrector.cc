// BeamEnergyCorrector implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// Local includes 
#include "BeamEnergyCorrector.h"

// DEPFETTrackTools includes
#include "TBTrack.h"
#include "TrackInputProvider.h"
#include "GenericTrackFitter.h"

// C++ includes
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <Exceptions.h>
#include <IMPL/LCFlagImpl.h>




// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;
using namespace CLHEP;

namespace depfet {

//
// Instantiate this object
//
BeamEnergyCorrector aBeamEnergyCorrector ;

//
// Constructor
//
BeamEnergyCorrector::BeamEnergyCorrector() : Processor("BeamEnergyCorrector")
{
   
// Processor description
  _description = "BeamEnergyCorrector: Applies corrections to the beam energy profile";
   

//
// Input collections  
  registerInputCollection(LCIO::TRACK,"InputTrackCollectionName",
                          "Track input collection",
                          _inputTrackCollectionName,std::string("tracks"));

//
// Output collections  
  registerOutputCollection(LCIO::TRACK,"OutputTrackCollectionName",
                           "Collection name for corrected tracks",
                           _outputTrackCollectionName, string ("corrtracks"));
   
// 
// Processor parameters
  
  registerProcessorParameter ("AlignmentDBFileName",
                             "This is the name of the LCIO file with the alignment constants (add .slcio)",
                             _alignmentDBFileName, static_cast< string > ( "alignmentDB.slcio" ) ); 
  
  registerProcessorParameter ("ReferencePlane",
                              "Plane number along the beam line",
                              _iref,  static_cast < int > (3));
  
  
  registerProcessorParameter ("CentralMomentum",
                              "Momentum at centre of reference plane [GeV]",
                              _M0,  static_cast < double > (3.));

  registerProcessorParameter ("SlopeUMomentum",
                              "Momentum slope along u [GeV/mm]",
                              _Mu,  static_cast < double > (0.0));
   
  registerProcessorParameter ("SlopeVMomentum",
                              "Momentum slope along v [GeV/mm]",
                              _Mv,  static_cast < double > (0.0));
                                 
}

//
// Method called at the beginning of data processing
//
void BeamEnergyCorrector::init() {
  
  // Initialize variables
  _nRun = 0 ;
  _nEvt = 0 ;
  
  // Print set parameters
  printProcessorParams();
  
  // CPU time start
  _timeCPU = clock()/1000;
  
  // Read detector constants from gear file
  _detector.ReadGearConfiguration();    
  
  // Read alignment data base file 
  _detector.ReadAlignmentDB( _alignmentDBFileName );       
     
}

//
// Method called for each run
//
void BeamEnergyCorrector::processRunHeader(LCRunHeader * run)
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
void BeamEnergyCorrector::processEvent(LCEvent * evt)
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
  
  GenericTrackFitter TrackFitter(_detector);
  TrackFitter.SetNumIterations(2); 
   
  LCCollection* inputCollection;
  try {
      inputCollection = evt->getCollection(_inputTrackCollectionName);
  } catch (DataNotAvailableException& e) {
      throw SkipEventException(this);
  }

  // Create output track collection
  LCCollectionVec * outputCollection = new LCCollectionVec(LCIO::TRACK);
    
  // Set flag for storing track hits in track collection
  LCFlagImpl flag(outputCollection->getFlag());
  flag.setBit( LCIO::TRBIT_HITS );
  outputCollection->setFlag(flag.getFlag());
  
  
  // Main loop over all tracks
  int nTracks = inputCollection->getNumberOfElements(); 
  
  for (int itrk = 0; itrk < nTracks; itrk++) {
    
    // Retrieve track from LCIO 
    Track * inputtrack = dynamic_cast<Track*> (inputCollection->getElementAt(itrk));
    
    // Convert LCIO -> TB track  
    TBTrack track = TrackIO.MakeTBTrack( inputtrack, _detector );  
    
    // ReFit track 
    bool trkerr = TrackFitter.Fit(track);
    if ( trkerr ) {
      continue;
    }  
    
    // Get fit results at reference plane 
    TBTrackElement& TE = track.GetTE(_iref);
    if ( !TE.IsCrossed() ) continue;
    
    double u = TE.GetState().GetPars()[2][0];
    double v = TE.GetState().GetPars()[3][0];
     
    // Compute corrected momentum
    double corrmom = _M0 + _Mu * u + _Mv * v;     
    track.SetMomentum( corrmom ); 
    track.GetReferenceState().Pars[4][0] = track.GetCharge() / corrmom ;
     
    // Convert TBTrack to LCIO::Track  
    TrackImpl* outputtrack = TrackIO.MakeLCIOTrack( track );
       
    // Add lcio track to lcio collection 
    outputCollection->addElement(outputtrack);
    
    
                               
  } // End loop over all tracks 
        
  evt->addCollection(outputCollection, _outputTrackCollectionName); 

}

//
// Method called after each event to check the data processed
//
void BeamEnergyCorrector::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void BeamEnergyCorrector::end()
{
   
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
void BeamEnergyCorrector::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "BeamEnergyCorrector Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}

} // Namespace
