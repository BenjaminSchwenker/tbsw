#ifndef AsciiInputProcessor_H
#define AsciiInputProcessor_H 1


// marlin includes ".h"
#include <marlin/DataSourceProcessor.h>
#include <marlin/Exceptions.h>

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>
#include <vector>
#include <string>

namespace depfet
{
  
  /** The AsciiInputProcessor Processor
   *
   * The input processor reads one or more ascii files produced containing 
   * Hits tables. It reads hits one by one from file, builds 
   * events (lcio::LCEvent) and fills a lcio::LCCollection with hits.  
   * 
   * Author: Benjamin Schwenker, Göttingen University 
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
  
  class AsciiInputProcessor:public marlin::DataSourceProcessor 
  {
   public:
     
    //! Default constructor
    AsciiInputProcessor();
    virtual AsciiInputProcessor * newProcessor ();
    virtual void readDataSource (int Ntrig);
    virtual void init ();
    virtual void end ();
    
    // read one event from file
    int getEventFromFile(std::string fileName, int currEvt,  std::vector< std::vector<int> >& hits);
    
    // read all event from file
    int getAllEventFromFile(std::string fileName, std::vector< std::vector<int> >& hits, std::vector< std::vector<int> >& events );
    
   protected:
     
    // Processor parameters 
    std::vector< std::string >  m_fileNameVec;  
    std::string m_detectorName;
    
   private:
    CellIDEncodeConstructHelper _outputEncoderHelper;


    
     
  };
  
  AsciiInputProcessor gAsciiInputProcessor;
 
}            
#endif
