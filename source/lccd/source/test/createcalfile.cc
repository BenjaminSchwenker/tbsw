#include <iostream>

// --- lcio headers 
#include "lcio.h"
#include "IO/LCWriter.h"
#include "IMPL/LCRunHeaderImpl.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCCollectionVec.h"
// ------------------

#include "CalibrationConstant.hh"


static int NCELL = 3 ;

using namespace std ;
using namespace lcio;


/** Test program that writes some calibration constants to an LCIO file.<br>
 *  NB: This uses LCIO only - no LCCD code !!!
 * 
 * @author F.Gaede, DESY
 * @version $Id: createcalfile.cc,v 1.1.1.1 2005-02-10 12:17:54 gaede Exp $
 */

int main(int argc, char** argv ) {
  
  // read file name and collection name from command line 
  if( argc < 3) {
    cout << " usage: createcalfile <output-file> <collection name>" << endl ;
    exit(1) ;
  }
  
  string fileName( argv[1] ) ;
  string colName( argv[2] ) ;
  
  
  LCWriter* wrt = LCFactory::getInstance()->createLCWriter() ;
  
  wrt->open( fileName , LCIO::WRITE_NEW )  ;


  // create and add a run header to the file 
  // this isn't needed by LCCD but allows to use the 'dumpevent' tool

  LCRunHeader* rHdr = new LCRunHeaderImpl ;
  wrt->writeRunHeader( rHdr ) ;

  //--------- create event with one collection --------------

  
  LCEvent* evt = new LCEventImpl ;
  
  LCCollectionVec* col = new LCCollectionVec( LCIO::LCGENERICOBJECT ) ;

  
  for(int i=0; i< NCELL ; i++ ) {
    
    int cellId   = i<<16 ;
    float offset =  i/100. ;
    float gain   = 3.141592 / (i+.0000001 )  ;
    
    CalibrationConstant* calCon = new CalibrationConstant( cellId, offset, gain  ) ;
    
    col->addElement( calCon ) ;

  }

  evt->addCollection(  col , colName ) ;
  
  evt->parameters().setValue("isCalibration", 1 ) ;

  wrt->writeEvent( evt ) ; 

  // --- clean up ---

  delete evt ;

  wrt->close() ;

  //----------------

}
