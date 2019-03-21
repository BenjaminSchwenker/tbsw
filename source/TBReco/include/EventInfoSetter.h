#ifndef EventInfoSetter_H
#define EventInfoSetter_H 1

// personal includes ".h"

// marlin includes ".h"
#include "marlin/DataSourceProcessor.h"

// lcio includes <.h>

// system includes <>


namespace depfet
{
  
  /** The EventInfoSetter Processor
   * Creates an empty lcio::Event and sets the event number. 
   * 
   * Author: Benjamin Schwenker, Göttingen University 
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
  
  class EventInfoSetter:public marlin::DataSourceProcessor 
  {
   public:
     
    //! Default constructor
    EventInfoSetter();
    virtual EventInfoSetter * newProcessor ();
    virtual void readDataSource (int Ntrig);
    virtual void init ();
    virtual void end ();
    
   protected:
     
    // Processor parameters 
    int m_runNumber;
    std::string m_detectorName;
     
  };
  
  EventInfoSetter gEventInfoSetter;
 
}            
#endif
