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
  _description = "Geometry: Setup geometry description for test beam telescope.";
   
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
   
  // Print set parameters
  printProcessorParams();
  
  std::string gearFile = Global::parameters->getStringVal("GearXMLFile" ) ;
   
  // Read detector constants from gear file NEW WAY
  // This creates an instance of TBDetector which is globally available
  TBDetector::GetInstance().ReadGearConfiguration(gearFile);
  
  streamlog_out( MESSAGE3 ) << "Print after reading gear " << endl;
  TBDetector::GetInstance().Print();
  
  // Read alignment data base file 
  if(!_newAlignment) TBDetector::GetInstance().ReadAlignmentDB( _alignmentDBFileName );
  // This is needed, because if the AlignmentDB is not read, the detector construct doesn't know the alignmentDB name
  else  TBDetector::GetInstance().SetAlignmentDBName( _alignmentDBFileName ); 
       
}


//
// Method called after all data processing
//
void Geometry::end()
{
  
  if ( _updateAlignment ) { 
    TBDetector::GetInstance().WriteAlignmentDB( ); 
  } 
           
  // Print message
  streamlog_out(MESSAGE3) << std::endl
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
