#ifdef USE_CONDDB

// -- LCCD headers
#include "lccd.h"
#include "lccd/DBInterface.hh"

// -- C++ headers 
#include <iostream>

using namespace std ;
using namespace lcio;
//using namespace lccd ;


/** Test program that creates an LCIO fole with all conditions data in a database folder 
 *  for a given tag.
 * 
 * @author F.Gaede, DESY
 * @version $Id: createdbfile.cc,v 1.2 2005-02-18 16:51:55 gaede Exp $
 */

int main(int argc, char** argv ) {
  
  // enable LCIO exception handling (usually automatically done when Reader/Writer exists)
  HANDLE_LCIO_EXCEPTIONS ;
  
  // read file name and collection name from command line 
  if( argc < 2) {
    cout << " usage: createdbfile <db collection name> [<tag>]" << endl ;
    exit(1) ;
  }
  
  string colName( argv[1] ) ;
  string folder( "/lccd/" + colName ) ;
  string tag("") ;
  
  if( argc > 2 )
    tag = argv[2] ;
  

  lccd::DBInterface db( folder , false ) ;

  db.createDBFile( tag ) ; 


}
#endif
