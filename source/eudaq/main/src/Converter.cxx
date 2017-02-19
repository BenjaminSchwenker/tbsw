//!  Converter - EUDAQ raw -> LCIO data conversion 
/*!  This program reads the data produced by the EUDAQ for the 
 *   EUDET telescope in its native format and converts data to 
 *   the LCIO format. The LCIO format has the appropriate data 
 *   model for event reconstruction.
 *   
 *   The Converter uses a plugin mechanism to handle data from 
 *   different types of pixel detectors. A typical telescope 
 *   setup will consist of a number of reference planes and at 
 *   least one DUT sensor. Every detector group should supply a 
 *   DataConverterPlugin to convert their data into approbriate 
 *   LCIO collections.   
 *   
 *   @author  Benjamin Schwenker, Universitaet Goettingen
 *   <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
  
  
// eudaq includes 
#include "eudaq/FileReader.hh"
#include "eudaq/DetectorEvent.hh"
#include "eudaq/RawDataEvent.hh"
#include "eudaq/PluginManager.hh"
#include "eudaq/OptionParser.hh"
#include "eudaq/counted_ptr.hh"
#include "eudaq/Utils.hh"
#include "eudaq/FileNamer.hh"
#include "eudaq/Exception.hh"

// stl includes 
#include <iostream>

// lcio includes 
#include "lcio.h"
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCFlagImpl.h>
#include "IO/ILCFactory.h"
#include "EVENT/LCIO.h"
#include "Exceptions.h"
#include "IMPL/LCTOOLS.h"
#include "IO/LCWriter.h"


using namespace std;
using namespace eudaq;
using namespace lcio;

namespace {
    static const unsigned TLUID = Event::str2id("_TLU");
    static const unsigned IDMASK = 0x7fff;
  }

std::vector<unsigned> parsenumbers(const std::string & s) {
  std::vector<unsigned> result;
  std::vector<std::string> ranges = split(s, ",");
  for (size_t i = 0; i < ranges.size(); ++i) {
    size_t j = ranges[i].find('-');
    if (j == std::string::npos) {
      unsigned v = from_string(ranges[i], 0);
      result.push_back(v);
    } else {
      long min = from_string(ranges[i].substr(0, j), 0);
      long max = from_string(ranges[i].substr(j+1), 0);
      if (j == 0 && max == 1) {
        result.push_back((unsigned)-1);
      } else if (j == 0 || j == ranges[i].length()-1 || min < 0 || max < min) {
        std::cout << "Bad range" << std::endl ;
      } else {
        for (long n = min; n <= max; ++n) {
          result.push_back(n);
        }
      }
    }
  }
  return result;
}

int main(int /*argc*/, char ** argv) {
   
  eudaq::OptionParser op("EUDAQ File Converter", "2.0", 
                         "A command-line tool for converting a raw data file to lcio file",
                         1);
  eudaq::Option<std::string> events(op, "e", "events", "", "numbers", "Event numbers to convert (eg. '1-10,99' default is all)");
  eudaq::Option<std::string> ipat(op, "i", "inpattern", "../data/run$6R.raw", "string", "Input filename pattern");
  eudaq::Option<std::string> opat(op, "o", "outpattern", "run$6R$X", "string", "Output filename pattern");
  eudaq::Option<std::string> level(op, "l", "log-level", "INFO", "level", "The minimum level for displaying log messages locally");
  
  bool debug = true; 
  
  try {
    op.Parse(argv);
    std::vector<unsigned> numbers = parsenumbers(events.Value());
    
    // Loop over data files 
    
    for (size_t i = 0; i < op.NumArgs(); ++i) {
                     
      //eudaq::FileReader reader(op.GetArg(i), ipat.Value(), true);
      eudaq::FileReader reader(op.GetArg(i), ipat.Value(), false);
      std::cout << "Reading: " << reader.Filename() << std::endl ;
      
      string m_filepattern = opat.Value();  
      int m_runnumber = reader.RunNumber();    
            
      // Create lcio file writer here
      LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
      
      // Open a new file
      try {
        lcWriter->open(FileNamer(m_filepattern).Set('X', ".slcio").Set('R', m_runnumber),lcio::LCIO::WRITE_NEW) ;
      } catch(const lcio::IOException & e) {
        std::cout << e.what() << std::endl ;
      }
      
      unsigned ndata = 0, ndatalast = 0, nnondet = 0, nbore = 0, neore = 0;       
      
      do {
        const eudaq::Event & ev = reader.GetEvent();
        
        if (ev.IsBORE()) {
          nbore++;
          if (nbore > 1) {
            std::cout << "Multiple BOREs " << to_string(nbore) << std::endl ;
          }
          
          
          // Process BORE event 
          if (const eudaq::DetectorEvent * dev = dynamic_cast<const eudaq::DetectorEvent *>(&ev)) {
             
            // Initialize pugin manager
            eudaq::PluginManager::Initialize(*dev);
            
            // Write LCIO run header here 
            IMPL::LCRunHeaderImpl* lcHeader = new IMPL::LCRunHeaderImpl;
            lcHeader->setDescription(" Reading from file " + reader.Filename());
            lcHeader->setRunNumber((*dev).GetRunNumber());
            lcHeader->setDetectorName("EUTelescope");
            
            try { 
              lcWriter->writeRunHeader(lcHeader);
            } catch(const lcio::IOException & e) {
              std::cout << e.what() << std::endl ;
            }
            delete lcHeader;  
          }
          
           
        } else if (ev.IsEORE()) {
          neore++;
          if (neore > 1) {
            std::cout << "Warning: Multiple EOREs: " << to_string(neore) << std::endl;
          }
          
        } else {
          ndata++;
                 
           
           
          if (const eudaq::DetectorEvent * dev = dynamic_cast<const eudaq::DetectorEvent *>(&ev)) {
                     

            
            for (size_t i = 0; i < (*dev).NumEvents(); ++i) {
              const eudaq::Event* subev = (*dev).GetEvent(i);           
               
              cout << "  TrigID  " << (eudaq::PluginManager::GetTriggerID(*subev) & IDMASK) 
                   << "  (" << subev->GetSubType() << ")"         
                   << endl;  
            } 
              
            
            LCEvent * lcEvent = eudaq::PluginManager::ConvertToLCIO(*dev);
            
            if ( lcEvent == NULL ) {
              std::cout << "The eudaq plugin manager is not able to create a valid LCEvent" << endl;
              exit(-1);
            }
             
            try { 
                lcWriter->writeEvent( lcEvent );
            } catch(const lcio::IOException & e) {
                std::cout << e.what() << std::endl ;
            }
            
            delete lcEvent;
            
          } 
          
          
        }
          
      } while (reader.NextEvent());
      
       
      try { 
        lcWriter->close(); 
      } catch(const lcio::IOException & e) {
        std::cout << e.what() << std::endl ;
      }        

      
      //std::cout << "Number of data events: " << to_string(ndata) << std::endl;
      if (!nbore) {
        std::cout << "Warning: No BORE found" << std::endl;
      } else if (nbore > 1) {
        std::cout << "Warning: Multiple BOREs found: " << nbore << std::endl;
      }
            
      if (!neore) {
        std::cout << "Warning: No EORE found, possibly truncated file." << std::endl;
      } else if (neore > 1) {
        std::cout << "Warning: Multiple EOREs found: " << nbore << std::endl;
      }
            
      if (nnondet || (nbore != 1) || (neore > 1)) std::cout << "Probably corrupt file." << std::endl;
        
    }
    
  } catch (...) {
    return op.HandleMainException();
  }
  return 0;
}
