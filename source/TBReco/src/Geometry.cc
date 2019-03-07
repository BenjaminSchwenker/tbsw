// Geometry Processor 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// TBTools includes 
#include "Geometry.h"

// C++ includes
#include <iostream>
#include <iomanip>


// Used namespaces
using namespace std; 
using namespace marlin;

namespace depfet {

//
// Instantiate this object
//
Geometry aGeometry ;

//
// Constructor
//
Geometry::Geometry() : Processor("Geometry")
{
   
  // Processor description
  _description = "Geometry: ";
   
  // 
  // Processor parameters
  registerProcessorParameter ("AlignmentDBFileName",
                             "This is the name of the file with the alignment constants (add .root)",
                             _alignmentDBFileName, static_cast< string > ( "alignmentDB.root" ) ); 
  
  registerProcessorParameter ("UpdateAlignment",
                              "Update lcio alignmentDB using alignment results (true/false)?",
                              _updateAlignment, static_cast <bool> (true) ); 

  registerProcessorParameter ("NewAlignment",
                              "Start alignment from scratch (true/false)?",
                              _newAlignment, static_cast <bool> (false) ); 
                                
}

//
// Method called at the beginning of data processing
//
void Geometry::init() {
  
  // Initialize variables
  _nRun = 0 ;
  _nEvt = 0 ;
  
  // Print set parameters
  printProcessorParams();
  
  // CPU time start
  _timeCPU = clock()/1000;
  
  std::string gearFile = Global::parameters->getStringVal("GearXMLFile" ) ;
  
  streamlog_out( MESSAGE3 ) << "Print before reading gear file " << endl;
  _detector.Print();

  // Read detector constants from gear file OLD WAY
  //_detector.ReadGearConfiguration();    
   
  // Read detector constants from gear file NEW WAY
  _detector.ReadGearConfiguration(gearFile);
  
  streamlog_out( MESSAGE3 ) << "Print after reading gear " << endl;
  _detector.Print();
  
  /*
  // Read alignment data base file 
  if(!_newAlignment) _detector.ReadAlignmentDB( _alignmentDBFileName );
  // This is needed, because if the AlignmentDB is not read, the detector construct doesn't know the alignmentDB name
  else  _detector.SetAlignmentDBName( _alignmentDBFileName ); 
  */     
}

//
// Method called for each run
//
void Geometry::processRunHeader(LCRunHeader * run)
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
void Geometry::processEvent(LCEvent * evt)
{
    
  //////////////////////////////////////////////////////////////////////  
  // Process next event
  ++_nEvt;
   
  if ( _nEvt % 1000 == 0 ) {
    streamlog_out( MESSAGE3 ) << "Processing event "
                              << evt->getEventNumber() << " in run "
                              << evt->getRunNumber() << endl; 
  }      
}

//
// Method called after each event to check the data processed
//
void Geometry::check( LCEvent * )
{
}

//
// Method called after all data processing
//
void Geometry::end()
{
  // TODO This will break in the new code 
  //TBDetector tmp_detector = _detector;
  
  /* 
  KalmanAlignmentAlgorithm2 Aligner;
  AlignableDet reco_const = Aligner.Fit(tmp_detector, alignment_data, AlignState, _maxTracks, _annealingTracks, _annealingFactor,  _pValueCut, _deviationCut, _useBC, _logLevel );
  
  bool error_fim = Aligner.AlignDetector(tmp_detector, reco_const);
  if ( error_fim ) {
    streamlog_out ( MESSAGE3 ) << "Alignment failed!" << endl;
  } 
  */
  
  /*
  // Create aligned detector   
  // TODO This will break in the new code 
  //_detector = tmp_detector; 
   
  if ( _updateAlignment ) { 
    _detector.WriteAlignmentDB( ); 
  } 
  */   
       
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
void Geometry::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "Geometry Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}

} // Namespace
