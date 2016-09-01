#include <iostream>

// LCIO headers
#include "lcio.h"
#include "IMPL/LCEventImpl.h"
#include "UTIL/LCTOOLS.h"
#include "UTIL/LCTime.h"

// LCCD headers
#include "lccd.h"
#include "lccd/DBFileHandler.hh"
#include "lccd/IConditionsChangeListener.hh"


#include "CalibrationConstant.hh"


using namespace std ;
using namespace lcio;
//using namespace lccd ;



/** Test class that simply dumps the conditions data collection whenever it changed*/
class TestListener : public lccd::IConditionsChangeListener {
  
public:
  TestListener( LCEventImpl** evt ) : _evt( evt ) {} 
  
  void conditionsChanged( LCCollection* col ) {
    
    std::string colName = col->getParameters().getStringVal("CollectionName") ;

    std::cout << " TestListener::conditionsChanged for event time stamp : "
	      << (*_evt)->getTimeStamp() 
	      << " collection name : "  << colName
	      <<  std::endl ;
    
    LCTOOLS::printLCGenericObjects( col ) ;
  }

  // need pointer to pointer to retireve current time stamp
  LCEventImpl** _evt ;
} ;


/** Test program that reads calibration constants from a db-file.
 *  Note: This assumes that the time stamps are not real but reasonably
 *  small numbers, e.g 0-1000.
 *   
 * 
 * @author F.Gaede, DESY
 * @version $Id: testdbfile.cc,v 1.2 2005-02-22 18:55:22 gaede Exp $
 */

int main(int argc, char** argv ) {
  

  // read file name and collection name from command line 
  if( argc < 5) {
    cout << " usage: testdbfile <input-file> <collection name> <timestamp0> <timestamp1>" << endl ;
    exit(1) ;
  }
  
  string fileName( argv[1] ) ;
  string colName( argv[2] ) ;

  lccd::LCCDTimeStamp t0 = atoll( argv[3] ) ;
  lccd::LCCDTimeStamp t1 = atoll( argv[4] ) ;


  // create a conditions handler
  // NB: the rest of the code uses the abstract interface !
 
  lccd::IConditionsHandler* conData = new lccd::DBFileHandler( fileName, colName ) ;
  

  LCEventImpl*  evt ;
  TestListener tl( & evt ) ;

  conData->registerChangeListener( &tl ) ;
 

  // simulate events in time intervall - one per ns (:-)
  for( lcio::long64 i=t0 ; i<t1 ; i++){

    evt = new LCEventImpl  ;
//     *evtptr = evt ;


    std::cout << " processing timestamp : " << i << std::endl ; 

    evt->setTimeStamp( i ) ;

    conData->updateEvent(  evt  )  ;

    delete evt ;
  }

  // clean up
  delete conData ;

}
