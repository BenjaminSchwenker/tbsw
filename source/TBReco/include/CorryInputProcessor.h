
#ifndef CorryInputProcessor_H
#define CorryInputProcessor_H 1


// marlin includes ".h"
#include <marlin/DataSourceProcessor.h>
#include <marlin/Exceptions.h>

// system includes <>
#include <vector>
#include <string>
#include <map>
#include <fstream>

namespace depfet
{
  
  /** The CorryInputProcessor Processor
   *
   * The input processor reads one or more ascii files containing 
   * hit tables. It reads all hits one by one from file, builds 
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
   * In addition, the 'event_number' must non negative and increase
   * monotonically. 
   * 
   * This processor provides a simple means to feed data from 
   * different subdetectors into tbsw.  
   * 
   * Author: Benjamin Schwenker, Göttingen University 
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
  
  class CorryInputProcessor:public marlin::DataSourceProcessor 
  {
   public:
     
    //! Default constructor
    CorryInputProcessor();
    virtual CorryInputProcessor * newProcessor ();
    virtual void readDataSource (int Ntrig);
    virtual void init ();
    virtual void end ();
    
   
    // read all event from file
    bool getAllEventFromFile(std::fstream& fin, std::vector< std::vector<int> >& hits, std::vector< std::vector<int> >& events, int max_events, int max_event_number);
    
   protected:
     
    // Processor parameters 

    //! Input hit files 
    std::vector< std::string >  m_fileNameVec;

    //! List of sensor names 
    std::vector< std::string >  m_sensorNamesVec;

    //! List of sensorIDs
    std::vector< int >  m_sensorIDsVec;

    //! Gear detector name   
    std::string m_detectorName;

    //! Output hit collection name
    std::string m_rawHitCollectionName;

    //! Run number 
    int m_runNumber;

   private:  
    //! Number of events per chunk
    int m_chunkSize;
    //! Number of processed triggers (events)
    unsigned m_ntriggers;
    //! Map from name of senor to sensorID 
    std::map<std::string, int> m_sensor_lookup;
  };
  
  CorryInputProcessor gCorryInputProcessor;
 
}            
#endif
