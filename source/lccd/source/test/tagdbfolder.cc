#ifdef USE_CONDDB

// -- CondDB headers
#include "ConditionsDB/ICondDBMgr.h"
#include "ConditionsDB/CondDBMySQLMgrFactory.h"
#include "ConditionsDB/CondDBObjFactory.h"
//#include "exampleObject.h"

// -- LCIO headers
#include "lcio.h"
#include "UTIL/LCTime.h"

// -- LCCD headers
#include "lccd.h"
#include "lccd/DBInterface.hh"

// -- C++ headers 
#include <iostream>

using namespace std ;
using namespace lcio;


/** Test program that tags an existing folder in the database.
 * 
 * @author F.Gaede, DESY
 * @version $Id: tagdbfolder.cc,v 1.2 2005-02-18 16:51:55 gaede Exp $
 */

int main(int argc, char** argv ) {
  
  // enable LCIO exception handling (usually automatically done when Reader/Writer exists)
  HANDLE_LCIO_EXCEPTIONS ;

  // read file name and collection name from command line 
  if( argc < 4) {
    cout << " usage: tagdbfolder <folder> <tag> <description> [<dbInit>]" << endl ;
    exit(1) ;
  }
  
  string folder     ( argv[1]  ) ;
  string tag        ( argv[2]  ) ;
  string description( argv[3]  ) ;
  


  string dbInit("") ;
  if( argc > 4 ) 
    dbInit = argv[4] ;
  else
    dbInit = lccd::getDBInitString() ; // "localhost:lccd_test:calvin:hobbes"  ;
  
  lccd::DBInterface db( dbInit, folder , true ) ;
  
  db.tagFolder( tag, description )  ;
  
  lcio::LCTime t  ;
  std::cout << " Tagged folder " << folder << " as [" << tag << "] " 
	    << " on " << t.getDateString() 
	    << std::endl ; 

}

#endif

