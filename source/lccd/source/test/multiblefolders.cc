#ifdef USE_CONDDB

// -- no CondDB headers needed !

// -- LCIO headers
#include "lcio.h"
#include "lccd/SimpleFileHandler.hh"

#include "UTIL/LCTOOLS.h"


// -- LCCD headers
#include "lccd.h"
#include "CalibrationConstant.hh"
#include "lccd/StreamerMgr.hh"
#include "lccd/VCollectionStreamer.hh"
#include "lccd/DBInterface.hh"

// -- C++ headers 
#include <iostream>
#include <sstream>
#include <cstdlib>

int NFolders = 10 ;

using namespace std ;
using namespace lcio;
//using namespace lccd ;


/** Test program that reads some calibration from a file
 *  and writes it to multible folders in the database testing 
 *  the reference counting for the DB connections.
 *  
 * @author F.Gaede, DESY
 * @version $Id: multiblefolders.cc,v 1.1 2006-11-14 17:34:27 gaede Exp $
 */

int main(int argc, char** argv ) {
  
  // read file name and collection name from command line 
  if( argc < 6) {
    cout << " usage: multiblefolders <input-file> <input collection name> <db collection name> <since> <till>" << endl ;
    exit(1) ;
  }
  
  string fileName( argv[1] ) ;
  string colName( argv[2] ) ;
  string colNameDB( argv[3] ) ;
  string sinceStr( argv[4] ) ;
  string tillStr( argv[5] ) ;
  
  // create a conditions handler
  // NB: the rest of the code uses the abstract interface !
 
  lccd::IConditionsHandler* conData = new lccd::SimpleFileHandler( fileName, colNameDB, colName ) ;
  

  // read data from file - timestamp is ignored for SimpleFileHandler
  conData->update( lccd::LCCDTimeStamp( 123456789) ) ;


  cout << endl 
       << " -- calibration data has been read from the file " 
       << fileName << " in collection " << colName 
       << endl 
       << endl ;

  LCCollection* col = conData->currentCollection() ;
  
  LCTOOLS::printLCGenericObjects( col ) ;
  
  lccd::DBInterface* db[NFolders] ;
  //------------------ write constants to database -----------------------
  for(int i=0 ; i< NFolders ; i++ ){

    stringstream str ;
    
    str << "/lccd/" <<  colNameDB << "_" << i ;

    string folder = str.str() ;
 
    db[i] = new lccd::DBInterface( folder , true ) ;
    
    lccd::LCCDTimeStamp since = std::atoll( sinceStr.c_str()  )  ;
    lccd::LCCDTimeStamp till  = std::atoll( tillStr.c_str()  )  ;
    
    db[i]->storeCollection( since, till, col , "test cond data" )  ;
  }

   for(int i=0 ; i< NFolders ; i++ ){
     delete db[i] ;
   }

  //------------------ end write constants to database -----------------------
  


  // clean up
  delete conData ;

}

#endif
