// ////////////////////////////////////////////////////////////////////////// //
//                                                                            //
//    Geometry - Marlin Processor                                             //
// ////////////////////////////////////////////////////////////////////////// //

#ifndef GEOMETRYPROCESSOR_H
#define GEOMETRYPROCESSOR_H 1

// TBTools includes   
#include "TBDetector.h"

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include basic C
#include <string>


namespace depfet {

//! Geometry processor 
/*! 
 *  Setup geometry description for test beam telescope 
 *   
 *  Author: Benjamin Schwenker, Göttingen University 
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */


class Geometry : public marlin::Processor {
  
 public:
  
   //!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new Geometry ; }
   
   //!Constructor - set processor description and register processor parameters
   Geometry();
   
   //!Method called at the beginning of data processing - used for initialization
   virtual void init();
   
   //!Method called after all data processing
   virtual void end();
   
   //!Method printing processor parameters
   void printProcessorParams() const;
      
 protected:
   
   //! Processor Parameters 
   
   //! AlignmentDB file path
   std::string m_alignmentDBFilePath;
      
   //! Override alignmentDB file   
   bool m_overrideAlignment;
   
   //! Apply alignmentDB corrections   
   bool m_applyAlignment;
    
}; // Class

} // Namespace

#endif 

