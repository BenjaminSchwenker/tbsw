// Geometry Processor 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

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
  registerProcessorParameter ("AlignmentDBFilePath",
                             "This is the path to the alignment constants (add .root)",
                             m_alignmentDBFilePath, static_cast< string > ( "alignmentDB.root" ) ); 
  
  registerProcessorParameter ("OverrideAlignment",
                              "Override alignmentDB (true/false)?",
                              m_overrideAlignment, static_cast <bool> (true) ); 

  registerProcessorParameter ("ApplyAlignment",
                              "Apply corrections from alignmentDB (true/false)?",
                              m_applyAlignment, static_cast <bool> (true) ); 
                                
}

//
// Method called at the beginning of data processing
//
void Geometry::init() {
   
  // Print set parameters
  printProcessorParams();
  
  std::string gearFile = Global::parameters->getStringVal("GearXMLFile" ) ;
   
  // Read detector constants from gear file and create a global instance
  // of TBDetector 
  TBDetector::GetInstance().ReadGearConfiguration(gearFile);
  
  // Set the path to alignmentDB file
  TBDetector::GetInstance().SetAlignmentDBPath( m_alignmentDBFilePath );
  
  // Read and apply alignment data base file 
  if(m_applyAlignment) TBDetector::GetInstance().ApplyAlignmentDB();
  
}


//
// Method called after all data processing
//
void Geometry::end()
{
  
  if ( m_overrideAlignment ) { 
    // Print message
    streamlog_out(MESSAGE3) << std::endl
                            << "Override the alignmentDB in path " << m_alignmentDBFilePath
                            << std::endl;
    
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
                            << "Geometry Processor development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}

} // Namespace
