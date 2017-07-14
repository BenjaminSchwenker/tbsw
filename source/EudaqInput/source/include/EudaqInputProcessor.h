#ifndef EudaqInputProcessor_H
#define EudaqInputProcessor_H 1

// personal includes ".h"

// marlin includes ".h"
#include <marlin/DataSourceProcessor.h>
#include <marlin/Exceptions.h>

// lcio includes <.h>

// system includes <>


namespace eudaqinput
{
  
  /** The EudaqInputProcessor Processor
   *
   * The input processor reads .raw files produced by EUDAQ event by 
   * event. A new lcio::LCEvent is created and the detector raw data
   * gets written into lcio::LCCollection objects. Subsquent processors
   * can be used to unpack and reconstruct the raw data collections. 
   * 
   * Author: Benjamin Schwenker, Göttingen University 
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
  
  class EudaqInputProcessor:public marlin::DataSourceProcessor 
  {
   public:
     
    //! Default constructor
    EudaqInputProcessor();
    virtual EudaqInputProcessor * newProcessor ();
    virtual void readDataSource (int Ntrig);
    virtual void init ();
    virtual void end ();
    
   protected:
     
    // Processor parameters 
    std::string m_filename;
    std::string m_detectorName;
    
    
   private:
    unsigned m_ndata;
    unsigned m_nbore; 
    unsigned m_neore;       
    
     
  };
  
  EudaqInputProcessor gEudaqInputProcessor;
 
}            
#endif
