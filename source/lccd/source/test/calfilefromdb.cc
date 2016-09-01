#ifdef USE_CONDDB

// -- CondDB headers
#include "ConditionsDB/ICondDBMgr.h"
#include "ConditionsDB/CondDBMySQLMgrFactory.h"
#include "ConditionsDB/CondDBObjFactory.h"
//#include "exampleObject.h"

// -- LCIO headers
#include "lcio.h"
#include "lccd/SimpleFileHandler.hh"
#include "IO/LCWriter.h"
#include "IMPL/LCRunHeaderImpl.h"
#include "IMPL/LCEventImpl.h"
//
#include "UTIL/LCTOOLS.h"
#include "UTIL/LCTime.h"


// -- LCCD headers
#include "lccd.h"
#include "lccd/StreamerMgr.hh"
#include "lccd/VCollectionStreamer.hh"
#include "lccd/DBCondHandler.hh"
#include "lccd/DBInterface.hh"
#include "lccd/ConditionsMap.hh"

#include "CalibrationConstant.hh"

// -- C++ headers 
#include <iostream>
#include <algorithm>

using namespace std ;
using namespace lcio;
//using namespace lccd ;


/** Test program that reads some calibration from a data base and writes 
 *  it to an LCIO file. Demonstrates the usage of lccd::DBCondHandler and
 *  lccd::ConditionsMap. <br>
 *  To simply create a file with conditions data
 *  for a certain point in time use lccd::DBInterface::createSimpleFile() 
 *  see commented out code at the end.
 * 
 * @author F.Gaede, DESY
 * @version $Id: calfilefromdb.cc,v 1.8 2005-04-22 14:26:51 gaede Exp $
 */

int main(int argc, char** argv ) {
  
  // enable LCIO exception handling (usually automatically done when Reader/Writer exists)
  HANDLE_LCIO_EXCEPTIONS ;

  // read file name and collection name from command line 
  if( argc < 4) {
    cout << " usage: calfilefromdb <db collection name> <timestamp> <output-file> [<tag>]" << endl ;
    exit(1) ;
  }

  string colName( argv[1] ) ;
  string timeStr( argv[2] ) ;
  string fileName( argv[3] ) ;

  lccd::LCCDTimeStamp  timeStamp ( atol( timeStr.c_str() ) )   ;  
  
  string folder( "/lccd/" + colName ) ;

  string tag("") ;

  if( argc > 4 )
    tag = argv[4] ;
  

  // testing: dump all collections in tag :
  //   lccd::DBInterface db("localhost:lccd_test:calvin:hobbes" , folder , false ) ;
  
  //   lccd::ColVec colVec ;
  //   db.findCollections( colVec , tag ) ;
  
  //   std::cout << " ---- Collections defined for tag: " << tag ;
  //   std::for_each(  colVec.begin() , colVec.end() , LCTOOLS::printLCGenericObjects ) ;
  //-----------  end dump all ----------------------------------


  // ---- use the DBCondHandler -----------

  lccd::IConditionsHandler* conData = 
    new lccd::DBCondHandler( lccd::getDBInitString() , folder, colName, tag ) ;
//     new lccd::DBCondHandler( "localhost:lccd_test:calvin:hobbes", folder, colName, tag ) ;
  

  // ------ testing: create a calibration map ------------------
  typedef lccd::ConditionsMap<int,CalibrationConstant> CalMap ;
  
  CalMap calMap( &CalibrationConstant::getCellID )   ;
  
  conData->registerChangeListener(  &calMap )  ;
  
  conData->update( timeStamp ) ;
  
  calMap.print( std::cout ) ;
  //--------------- end calibration map -----------------------



  lcio::LCTime t0 ( conData->validSince()  ) ;
  lcio::LCTime t1 ( conData->validTill()  ) ;

  cout << endl 
       << " -- calibration data has been read from data base folder" << folder 
       << endl 
       << " ---- valid from: " <<  t0.getDateString()  
       << " [" <<  conData->validSince()  << "] " 
       << endl 
       << " ----       till: " <<  t1.getDateString()  
       << " [" <<  conData->validTill()  << "] "
       << endl 
       << endl ;
   
  LCCollection* col = conData->currentCollection() ;
  
  //------------------ write constants to file -----------------------
  
  LCWriter* wrt = LCFactory::getInstance()->createLCWriter() ;
  
  wrt->open( fileName , LCIO::WRITE_NEW )  ;
  
  
  // create and add a run header to the file 
  // this isn't needed by LCCD but allows to use the 'dumpevent' tool
  LCRunHeader* rHdr = new LCRunHeaderImpl ;
  wrt->writeRunHeader( rHdr ) ;
  
  LCEvent* evt = new LCEventImpl ;
  evt->addCollection(  col , colName  ) ;  

  LCTOOLS::dumpEventDetailed( evt ) ; 
  wrt->writeEvent( evt ) ;
  
  // this is needed because the event should not own the collection when beeing deleted 
  // FIXME: not very nice, is it ?
  evt->takeCollection(  colName  ) ;  

  wrt->close() ;
  //----------------------------------------------------------------

  // clean up
  delete evt ; 
  delete wrt ;
  delete rHdr ;
  delete conData ;


  //---------------------------------------------------------------------------

  // while the above demonstrates and tests the use of DBCondHandler and ConditionsMap
  // we could have created the LCIO file much easier :)

  lccd::DBInterface db(  folder )  ;
  db.createSimpleFile( timeStamp , tag , true ) ;

  //---------------------------------------------------------------------------


}
#endif
