#include <iostream>

// LCIO headers
#include "lcio.h"
#include "UTIL/LCTOOLS.h"
#include "UTIL/LCTime.h"

// LCCD headers
#include "lccd.h"
#include "lccd/DataFileHandler.hh"


#include "CalibrationConstant.hh"


using namespace std ;
using namespace lcio;
//using namespace lccd ;

/** Test program that reads a data file and treats the given collection as 
 *  conditions data, i.e. adds it to every event unless a new version 
 *  exists in the event.
 * 
 * @author F.Gaede, DESY
 * @version $Id: readdatafile.cc,v 1.2 2005-02-22 18:55:22 gaede Exp $
 */

int main(int argc, char** argv ) {
  
  // read file name and collection name from command line 
  if( argc < 3) {
    cout << " usage: readdatafile <input-file> <collection name>" << endl ;
    exit(1) ;
  }
  
  string fileName( argv[1] ) ;
  string colName( argv[2] ) ;
  
  
  
  // create a conditions handler
  // NB: the rest of the code uses the abstract interface !
  
  lccd::IConditionsHandler* conData = new lccd::DataFileHandler( colName ) ;
  
  // now read the data file
  LCReader* lcReader = LCFactory::getInstance()->createLCReader() ;
  lcReader->open( fileName ) ;

  LCEvent* evt ;
  int nEvents = 0 ;

  //----------- the event loop -----------
  while( (evt = lcReader->readNextEvent()) != 0 ) {
    
    conData->updateEvent( evt ) ;
    
    LCTOOLS::dumpEvent( evt ) ;

    nEvents ++ ;
  }
  // -------- end of event loop -----------


  // clean up
  delete conData ;
  delete lcReader ;

}
