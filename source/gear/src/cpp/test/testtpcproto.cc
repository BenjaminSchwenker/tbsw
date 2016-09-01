#include "gearimpl/RectangularPadRowLayout.h"
#include "gearimpl/Util.h"
#include "gearxml/GearXML.h"
#include "gear/GearMgr.h"
#include "gear/GEAR.h"

#include <iostream>
#include <assert.h>

#include <exception>
#include <cstdlib>

#ifdef GEAR_USE_AIDA  
#include "AIDA/AIDA.h"
using namespace AIDA ;
#endif


using namespace gear ;

void testRectangularPadRowLayout( const RectangularPadRowLayout& pl ) ;


void gear_unexpected(){

  try {

    throw ;

  } catch( std::exception& e) {
    
    std::cout << " A runtime error has occured : "
              << e.what()
              << std::endl
              << " the program will have to be terminated - sorry." << std::endl ;
    exit(1) ;

  }

  catch(...) {
    
    std::cout << " An unknown runtime error has occured : "
              << std::endl
              << " the program will have to be terminated - sorry." << std::endl ;
    exit(1) ;
  }
}


/** Testprogram for TPC prototype description. If linked with an AIDA implementation a hit map is created that
 *  shows the pad layout.
 */

int main(int argc, char**argv){
  
  
  std::set_unexpected( gear_unexpected ) ;
  std::set_terminate( gear_unexpected ) ;
  
  if( argc < 2 ) {
    std::cout << " testgear:  Testprogram for TPC prototype description " << std::endl 
	      << " usage:     testtpcproto input.xml  " << std::endl ;
    exit(1) ;
   }
  
  std::string fileName( argv[1] ) ;
  
  GearXML gearXML( fileName ) ;
  
  GearMgr* gearMgr = gearXML.createGearMgr() ;
  
  
  std::cout << " testtpcproto - instantiated  GearMgr from file " << fileName 
	    << std::endl ;

  std::cout  << *gearMgr  << std::endl ;


  try{
    const RectangularPadRowLayout& pl = 
      dynamic_cast<const RectangularPadRowLayout&>( gearMgr->getTPCParameters().getPadLayout() );
    
    testRectangularPadRowLayout( pl ) ;
  }
  catch(std::bad_cast& e){
    std::cout << "  wrong type of layout - expected RectangularPadRowLayout ! " << std::endl ;
  }
  
  // --- test writing of XML file ---------

  GearXML::createXMLFile( gearMgr, "testgear_out.xml" ) ;

}



void testRectangularPadRowLayout( const  RectangularPadRowLayout& pl ) {
  

#ifdef GEAR_USE_AIDA  

  std::string storeName("aidaroot.root");

  IAnalysisFactory * myaida = AIDA_createAnalysisFactory() ;

  ITreeFactory * mytreefactory = myaida->createTreeFactory() ;

  ITree * mytree = mytreefactory->create(storeName,"root",false,true,"none") ;
  //              mytreefactory->create(storeName,storeType,readOnly,createNew,options)

  IHistogramFactory * myhistofactory = myaida->createHistogramFactory( *mytree ) ;

#endif

  double xmin = pl.getPlaneExtent()[0] ;
  double xmax = pl.getPlaneExtent()[1] ;
  double ymin = pl.getPlaneExtent()[2] ;
  double ymax = pl.getPlaneExtent()[3] ;


  int nRow = pl.getNRows() ;
  
  std::cout << "   RectangularPadRowLayout :  " << std::endl
 	    << "         xMin:      " << xmin   << std::endl
 	    << "         xMax:      " << xmax   << std::endl
 	    << "         yMin:      " << ymin   << std::endl
 	    << "         yMax:      " << ymax   << std::endl
 	    << "         padHeight: " <<  pl.getPadHeight(0)  << std::endl 
 	    << "         padWidth:  " <<  pl.getPadWidth(0)  << std::endl
 	    << "         nRows :    " <<  nRow << std::endl 
  	    << std::endl 
 	    << std::endl ;
  
#ifdef GEAR_USE_AIDA  
  IHistogram2D *padMap = myhistofactory->createHistogram2D("padMap",
							   500, xmin, xmax,
							   500, ymin, ymax );
#endif

  int nPadTotal = 0 ;

  for( int i = 0 ; i < nRow ; i++) {
    
    if( i==0 || i == nRow-1 ) 
      std::cout << " --------- row : " << i << std::endl ;
    
    const  std::vector<int>& pads = pl.getPadsInRow( i ) ;
    
    int nPad = pads.size() ;
    nPadTotal += nPad ;
    
    for( int j = 0 ; j < nPad ; j++) {
      
      int iRow = pl.getRowNumber( pads[j] ) ;
      int iPad = pl.getPadNumber( pads[j] ) ;
      
      if( j != 0 ) {

 	int ln = pl.getLeftNeighbour(  pl.getPadIndex( iRow , iPad ) ) ; 

 	assert(  pl.getPadNumber( ln ) ==  iPad - 1 ) ;
      }

      if( j != nPad-1 ) {

	int rn = pl.getRightNeighbour(  pl.getPadIndex( iRow , iPad ) ) ; 

 	assert(  pl.getPadNumber( rn )  ==  iPad + 1 ) ;
      }


      Vector2D p = pl.getPadCenter( pads[j] ) ;


#ifdef GEAR_USE_AIDA  
      for(unsigned k=0;k<1000;++k){
	double x = xmin + ( xmax - xmin ) * rand()/ double(RAND_MAX ) ;
	double y = ymin + ( ymax - ymin ) * rand()/ double(RAND_MAX ) ;
	
	if(  pl.isInsidePad( x , y ) )
	  padMap->fill(  x , y  ) ;

      }
#endif

      if( (i==0 && j < 10 ) || ( i == nRow-1 && j > nPad-9 ) ) {
	
	std::cout << "         pad: "  
		  << " [" << iRow << "," << iPad << "] "
		  << " - ( " << p[0] << " , " << p[1]  << ") "
		  << std::endl ;
      }
      
      assert(  pl.getNearestPad(  p[0] , p[1] ) == pads[j]  ) ;
      assert(  pl.isInsidePad(  p[0] , p[1] , pads[j] ) ) ;
    }
  }

  assert( nPadTotal ==  pl.getNPads() ) ;

#ifdef GEAR_USE_AIDA  
  mytree->commit();
  //  mytree->close();
#endif


  //---------------------------------
}
  
