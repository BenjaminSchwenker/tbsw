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
 * 
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
   
   //!Method called for each run - used for run header processing
   virtual void processRunHeader(LCRunHeader * run);
   
   //!Method called for each event - used for event data processing
   virtual void processEvent(LCEvent * evt);
   
   //!Method called after each event - used for data checking
   virtual void check(LCEvent * evt);
   
   //!Method called after all data processing
   virtual void end();
   
   //!Method printing processor parameters
   void printProcessorParams() const;
      
 protected:
   
   //! Processor Parameters 
      
   //! AlignmentDB file name 
   std::string _alignmentDBFileName;
      
   //! Update alignmentDB   
   bool _updateAlignment;
   
   //! New alignment  
   /*! Don't use current alignment data base, but start from scratch   
    */
   bool _newAlignment;
   
 private:
   
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
    
   // Handle to detector data 
   TBDetector  _detector;        
}; // Class

} // Namespace

#endif 

