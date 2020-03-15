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
   * The input processor reads one or more ascii files containing 
   * Hits tables. It reads all hits one by one from file, builds 
   * lcio::LCEvent and fills a lcio::LCCollection with hits.  
   * 
   * The input files contain comma seperated values with one hit 
   * per row. There should five columns containing 'event_number'
   * 'layerID', 'column', 'row' and 'charge'.  
   * 
   * The 'layerID' must be unique identifier for a detector layer
   * described in the gear file. 
   * 
   * The 'event_number' must be unique across all files and serves
   * to merge hits from different input files into global events.
   * 
   * This processor provides a simple means to feed data from 
   * different subdetectors into tbsw.  
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
