// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    ClusterShapeToAsciiPrinter - Marlin Processor                                         //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef ClusterShapeToAsciiPrinter_H
#define ClusterShapeToAsciiPrinter_H 1

// TBTools includes
#include "TBDetector.h"

// Include basic C
#include <vector>
#include <string>
#include <fstream>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"


namespace depfet {

  //! ClusterShapeToAsciiPrinter  
  /*! 
   *  The task of this processor is compute cluster shapes from sensor 
   *  clusters and print them to an ascii file. 
   *  
   *  Author: B.Schwenker, Universität Göttingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
   
  class ClusterShapeToAsciiPrinter : public marlin::Processor {
    
   public:
    
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new ClusterShapeToAsciiPrinter ; }
    
    //!Constructor - set processor description and register processor parameters
    ClusterShapeToAsciiPrinter();
    
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
     
   protected:
    
    //!Method printing processor parameters
    void printProcessorParams() const;
    
    // Processor Parameters
    
    //! Input track collection name
    std::string _inputTrackCollectionName;
    
    //! AlignmentDB file name 
    std::string _alignmentDBFileName;
       
    //! Output file name  
    std::string _asciiFileName;    
    
    //! Ignore clusters from these sensorIDs 
    std::vector<int >  _ignoreIDVec;
        
   private:
          
    // Handle to detector data 
    TBDetector  _detector;    

    // Handle to ascii file
    std::ofstream _outfile;
     
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
     
  }; // Class

} // Namespace

#endif 



