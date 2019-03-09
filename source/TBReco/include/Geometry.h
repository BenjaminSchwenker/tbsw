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

   //! Path for the geometry in the parameter space
   std::string m_geometryPath;

   /** Whether or not this module will raise an error if the geometry is
     * already present. This can be used to add the geometry multiple times if
     * it's not clear if it's already present in another path */
   bool m_ignoreIfPresent{false};

   /** If true we need to create a payload */
   bool m_createGeometryPayload{false};
      
   //! AlignmentDB file name 
   std::string _alignmentDBFileName;
      
   //! Update alignmentDB   
   bool _updateAlignment;
   
   //! New alignment  
   /*! Don't use current alignment data base, but start from scratch   
    */
   bool _newAlignment;
    
}; // Class

} // Namespace

#endif 

