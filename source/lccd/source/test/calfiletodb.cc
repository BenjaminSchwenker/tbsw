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
#include <cstdlib>

using namespace std ;
using namespace lcio;
//using namespace lccd ;


/** Test program that reads some calibration from a file
 *  and writes it to the database.
 * 
 * @author F.Gaede, DESY
 * @version $Id: calfiletodb.cc,v 1.4 2005-02-21 13:01:06 gaede Exp $
 */

int main(int argc, char** argv ) {
  
  // read file name and collection name from command line 
  if( argc < 6) {
    cout << " usage: calfiletodb <input-file> <input collection name> <db collection name> <since> <till>" << endl ;
    exit(1) ;
  }
  
  string fileName( argv[1] ) ;
  string colName( argv[2] ) ;
  string colNameDB( argv[3] ) ;
  string sinceStr( argv[4] ) ;
  string tillStr( argv[5] ) ;
  
  string folder( "/lccd/" + colNameDB ) ;


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
  //   // or use the class CalibrationConstant to print, i.e. access the data
  //   for(int i=0 ; i< col->getNumberOfElements() ; i++) {
  //     CalibrationConstant cc( col->getElementAt(i) ) ;
  //     cc.print( std::cout ) ; 
  //   }
  //   cout << endl 
  //        << " vailid since " <<  conData->validSince() 
  //        << "  till " <<  conData->validTill() 
  //        << endl 
  //        << endl ;
  
  
  //------------------ write constants to database -----------------------
  
  lccd::DBInterface db( folder , true ) ;
//   lccd::DBInterface db( lccd::getDBInitString() , folder , true ) ;
// lccd::DBInterface db("localhost:lccd_test:calvin:hobbes" , folder , true ) ;
  
  lccd::LCCDTimeStamp since = std::atoll( sinceStr.c_str()  )  ;
  lccd::LCCDTimeStamp till  = std::atoll( tillStr.c_str()  )  ;
  
  db.storeCollection( since, till, col , "test cond data" )  ;
  
  //------------------ end write constants to database -----------------------
  


  // clean up
  delete conData ;

}

#endif
