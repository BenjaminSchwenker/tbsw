#include "gearxml/tinyxml.h"
#include "gearxml/MergeXML.h"
#include "gear/GEAR.h"


/** Program to merge two xml files. This can be used if a gear xml 
 * (that has been automatically created with MokkaGear) is lacking 
 * some information. XML Elements that are present in only one of the input
 * files will be merged into the target file. In case the same element
 * exists in either input file the attribute values from the primary file
 * are taken. 
 * 
 * @author R.Lippe, DESY
 * @version $Id: mergeXML.cc,v 1.3 2007-11-16 17:10:30 gaede Exp $
 */


int main( int argc, char *argv[] ) {

  try{

    if( argc < 4 ) {

      std::cout << std::endl 
		<< " Program to merge two xml files. This can be used if a gear xml" << std::endl
		<< " (that has been automatically created with MokkaGear) is lacking" << std::endl
		<< " some information. XML Elements that are present in only one of the input" << std::endl
		<< " files will be merged into the target file. In case the same element"  << std::endl
		<< " exists in either input file the attribute values from the primary file" << std::endl
		<< " are taken." << std::endl
		<< std::endl
		<< "usage: " << std::endl
		<< "  mergeXML  <secondary>  <primary>  <target> " << std::endl
		<< "   <secondary>   filename for base-file that will be overwritten" << std::endl
		<< "   <primary>     filename for file that overrules secondary file" << std::endl
		<< "   <target>      filename for output file " << std::endl
		<< std::endl ;
      exit(1) ;
    }

    
    std::string fileName1   = argv[1] ;
    std::string fileName2   = argv[2] ;
    std::string fileNameNew = argv[3] ;
    
    gear::MergeXML mergeXML ;
    
    if ( mergeXML.setFile1 (fileName1) ) {
      std::cout << "secondary   : ok" << std::endl ;
    }
    else {
      std::cout << "secondary   : error reading file " << fileName1 << std::endl ;
      exit (2) ;
    }

    if ( mergeXML.setFile2 (fileName2) ) {
      std::cout << "primary: ok" << std::endl ;
    }
    else {
      std::cout << "primary: error reading file " << fileName2 << std::endl ;
      exit (2) ;
    }

    if ( mergeXML.mergeFiles(fileNameNew) ) {
      std::cout << "target  : ok" << std::endl ;
    }
    else {
      std::cout << "target: error merging files in " << fileNameNew << std::endl ;
      exit(1) ;
    }
    
    return 0;

  }
  catch( gear::Exception& e) {
    std::cout << " An Exception occured : " << e.what() << std::endl ;
  }

}
  

