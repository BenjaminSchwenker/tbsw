
#include "gearimpl/Util.h"
#include "gearxml/GearXML.h"
#include "gear/GearMgr.h"
#include "gear/GEAR.h"
#include "gear/VXDParameters.h"

#ifdef CGA
#include "gearcga/CGAGearDistanceProperties.h"
#include "gearcga/CGAGearPointProperties.h"
#endif

#include <iostream>
#include <assert.h>

#include <exception>
#include <cstdlib>

#include <sstream>
#include <fstream>

#include <math.h>

#ifdef GEAR_USE_AIDA  
#include "AIDA/AIDA.h"
using namespace AIDA ;
#endif



using namespace gear ;

void testVXD( const VXDParameters& vxdParams ) ;

void createCheckPlots( const VXDParameters& vxdParams ) ;

void testVXDPoint( const Vector3D p , const VXDParameters& vxdParams ) ;

void testVXDDist( const Vector3D p , const VXDParameters& vxdParams ) ;

void testAllVXDPoint( const VXDParameters& vxdParams ) ;

void testAllVXDDist( const VXDParameters& vxdParams ) ;

void testVXDIntersection( const VXDParameters& vxdParams ) ;

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
}


#ifdef GEAR_USE_AIDA  
// global variables for AIDA histograms and tuples

  IAnalysisFactory * myaida ;
  ITreeFactory * mytreefactory ;
  ITree * mytree ;
  IHistogramFactory * myhistofactory ;
  ITupleFactory * mytuplefactory ;

#endif





/** Testprogram for gear classes.
 * 
 */

int main(int argc, char**argv){
  
  
  std::set_unexpected( gear_unexpected ) ;
  std::set_terminate( gear_unexpected ) ;
  
  if( argc < 2 ) {
    std::cout << " testgear:  Testprogram for vertex detecto in gear. " << std::endl 
	      << " usage:     testgear input.xml  " << std::endl ;
    exit(1) ;
   }
  
  std::string fileName( argv[1] ) ;
  
  GearXML gearXML( fileName ) ;
  
  GearMgr* gearMgr = gearXML.createGearMgr() ;
  
  
  std::cout << " testVTXgear - instantiated  GearMgr from file " << fileName 
	    << std::endl ;

  const VXDParameters& vp = 
    dynamic_cast<const VXDParameters&>( gearMgr->getVXDParameters() ) ;
 

  std::cout << " VXD parameters : " << std::endl 
	    << vp 
	    << std::endl ;



#ifdef GEAR_USE_AIDA  

  std::string storeName("vxdCheckPlot.root");

  myaida = AIDA_createAnalysisFactory() ;

  mytreefactory = myaida->createTreeFactory() ;

  mytree = mytreefactory->create(storeName,"root",false,true,"none") ;
  //              mytreefactory->create(storeName,storeType,readOnly,createNew,options)

  myhistofactory = myaida->createHistogramFactory( *mytree ) ;
  mytuplefactory = myaida->createTupleFactory( *mytree ) ;
  
//   const VXDParameters& vp = 
//     dynamic_cast<const VXDParameters&>( gearMgr->getVXDParameters() ) ;

  //testVXDIntersection( vp ) ;
  
  createCheckPlots( vp ) ;
  
#else
  
  std::cout << " Compile with GEAR_USE_AIDA==1 to get some histograms   ! " 
	    << std::endl; 
  
#endif
  

  testVXD( vp ) ;

}

void testVXD( const VXDParameters& vxdParams ) {

  // Methods for special debugging

  Vector3D p0 ( 0, 0 , 0  ) ;
  Vector3D p1 ( 12 , 2 , 0 ) ;
//   Vector3D p1 ( -36.8 , 7 , 10 ) ;
  Vector3D p2 ( -15 ,-8, .5 ) ;
  Vector3D p3 ( -30 , -2 , .5 ) ;
  Vector3D p4 ( -25 ,-25.0 , 0 ) ;
 
//   testVXDPoint( p0 , vxdParams ) ;
//   testVXDPoint( p1 , vxdParams ) ;
//   testVXDPoint( p2 , vxdParams ) ;
//   testVXDPoint( p3 , vxdParams ) ;
//   testVXDPoint( p4 , vxdParams ) ;
 

  testVXDDist( p0, vxdParams ) ;
  testVXDDist( p1, vxdParams ) ;
  testVXDDist( p2, vxdParams ) ;
  testVXDDist( p3, vxdParams ) ;
  testVXDDist( p4, vxdParams ) ;
  
  // testAllVXDDist( vxdParams ) ;
  testAllVXDPoint( vxdParams ) ;
}

#ifdef GEAR_USE_AIDA

void createCheckPlots( const VXDParameters& vxdParams ) {
  
  // creates sophisticated checkplots
  
  // values for LoRes Overview
  double testStart = -65 ;
  double testEnd = 65 ;
  double testStep = .05 ;
  double zStart = 10 ;
  double zEnd = 11 ;
  double zStep = 1 ;

  // Every _skipDist_ Step the distance test is performed
  int skipDist = 1000 ;
  
  // values for HiRes Part
  double testDetailStart = 2 ;
  double testDetailEnd = 18 ;
  double testDetailStep = 0.005 ;

  float myStatusP = 0, myStatusL = 0;
  int nLadder = 0 , nSensitive = 0 ;

  // ladder and sensitive Map
  std::string mapName = "vxdMap" ;
  std::string mapTitle = "Map of the VTX" ;
  std::vector<std::string> columnNames ;
  std::vector<std::string> columnTypes ;

  columnNames.push_back( "x" )         ; columnTypes.push_back( "double" ) ; 
  columnNames.push_back( "y" )         ; columnTypes.push_back( "double" ) ; 
  columnNames.push_back( "z" )         ; columnTypes.push_back( "double" ) ; 
  columnNames.push_back( "sensitive" ) ; columnTypes.push_back( "bool" ) ;

  ITuple *vxdMap = mytuplefactory->create( mapName , mapTitle , columnNames , columnTypes ) ;

  // ladder and sensitive Detail
  mapName = "vxdDetailMap" ;
  mapTitle = "Detailed Map of part of the VTX" ;

  ITuple *vxdDetailMap = mytuplefactory->create( mapName , mapTitle , columnNames , columnTypes ) ;

  // map of points pointed on by VXDParametes::distanceToLadder()
  mapName = "vxdPointsOnSurface" ;
  mapTitle = "Map of Points On Surface" ;

  ITuple *vxdSurfaceMap = mytuplefactory->create( mapName , mapTitle , columnNames , columnTypes ) ;

  // map of point returned from intersection function between line and vxd
  mapName = "vxdIntersectionPoints" ;
  mapTitle = "Map of Intersection Points" ;

  ITuple *vxdIntersectionMap = mytuplefactory->create( mapName , mapTitle , columnNames , columnTypes ) ;

  // Histogram of Distances Surface to Ladder
  IHistogram1D *distanceHisto = myhistofactory->createHistogram1D("distance",
								  100, 0, 0.00000000000001 ) ;
  
  // number of point to be tested
  int nPoints = int( (testEnd-testStart)*(testEnd-testStart)*(zEnd-zStart) 
    / ( testStep * testStep * zStep ) 
    +  (testDetailEnd-testDetailStart)*(testDetailEnd-testDetailStart)*(zEnd-zStart) 
    / (testDetailStep * testDetailStep * zStep ) ) ;

  std::cout << "Testing " << nPoints << " Points to AIDA" << std::endl ;

  //counter for skipping
  int skipped = 0 ;
  
  // first run over normal map
  // x
  for ( double x = testStart ; x < testEnd ; x += testStep ) {
    
    // Status
    float curStatus = (x - testStart) / (testEnd - testStart) * 50 ;
    if ( curStatus > ( myStatusP + 0.25 ) ) {
      myStatusP = curStatus ;
      std::cout << "." << std::flush ;
    }
    if ( curStatus > ( myStatusL + 10 ) ) {
      myStatusL = curStatus ;
      std::cout << floor(myStatusL) << "% finished" << std::endl ;
    }

   	 
    // y
    for ( double y=testStart ; y < testEnd ; y += testStep ) {
      
      // z
      for ( double z = zStart ; z < zEnd ; z += zStep ) {
	
	Vector3D p ( x , y , z ) ;
	bool isPoint = vxdParams.isPointInLadder( p ) ;
	bool isSensitive = vxdParams.isPointInSensitive( p ) ;
	
	if( isPoint ) nLadder ++;
	if( isSensitive ) nSensitive ++;

	if( isPoint || isSensitive ) {

	  // regular map
	  vxdMap->fill( 0 , p[0] ) ;
	  vxdMap->fill( 1 , p[1] ) ;
	  vxdMap->fill( 2 , p[2] ) ;
	  vxdMap->fill( 3 , isSensitive ) ;
	  vxdMap->addRow() ;

	}
	
	if( skipped >= skipDist ) {
	  // points pointed on
	  // ladder
	  Vector3D v = vxdParams.distanceToNearestLadder( p ) ;
	  Vector3D np (  p[0] + v[0] , p[1] + v[1] , p[2] + v[2] ) ;
	  // sensitive
	  v = vxdParams.distanceToNearestSensitive( p ) ;
	  Vector3D nps (  p[0] + v[0] , p[1] + v[1] , p[2] + v[2] ) ;
	  
	  vxdSurfaceMap->fill( 0 , np[0] ) ;
	  vxdSurfaceMap->fill( 1 , np[1] ) ;
	  vxdSurfaceMap->fill( 2 , np[2] ) ;
	  vxdSurfaceMap->fill( 3 , false ) ;
	  vxdSurfaceMap->addRow() ;
	  
	  vxdSurfaceMap->fill( 0 , nps[0] ) ;
	  vxdSurfaceMap->fill( 1 , nps[1] ) ;
	  vxdSurfaceMap->fill( 2 , nps[2] ) ;
	  vxdSurfaceMap->fill( 3 , true ) ;
	  vxdSurfaceMap->addRow() ;
	  
	  // distance histogram
	  v = vxdParams.distanceToNearestLadder( np ) ;
	  distanceHisto->fill( sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] ) ) ;
	  
	  // reset counter
	  skipped = 0 ;
	}
	else {
	  skipped ++ ;
	}
	
      } // z
    } // y
  } //x


  // second run over detailed map
  
  // x
  for ( double x = testDetailStart ; x < testDetailEnd ; x += testDetailStep ) {
    
    // Status
    float curStatus = (x - testDetailStart) / (testDetailEnd - testDetailStart) * 50 + 50 ;
    if ( curStatus > ( myStatusP + 0.25 ) ) {
      myStatusP = curStatus ;
      std::cout << "." << std::flush ;
    }
    if ( curStatus > ( myStatusL + 10 ) ) {
      myStatusL = curStatus ;
      std::cout << floor(myStatusL) << "% finished" << std::endl ;
    }

   	 
    // y
    for ( double y=testDetailStart ; y < testDetailEnd ; y += testDetailStep ) {
      
      // z
      for ( double z = zStart ; z < zEnd ; z += zStep ) {
	
	Vector3D p ( x , y , z ) ;
	bool isPoint = vxdParams.isPointInLadder( p ) ;
	bool isSensitive = vxdParams.isPointInSensitive( p ) ;
	
	if( isPoint ) nLadder ++;
	if( isSensitive ) nSensitive ++;

	if( isPoint || isSensitive ) {

	  // regular map
	  vxdDetailMap->fill( 0 , p[0] ) ;
	  vxdDetailMap->fill( 1 , p[1] ) ;
	  vxdDetailMap->fill( 2 , p[2] ) ;
	  vxdDetailMap->fill( 3 , isSensitive ) ;
	  vxdDetailMap->addRow() ;

	}

      } // z
    } // y
  } //x
 

  // testing intersectionVXD
 
  std::cout << "Testing intersection " << std::flush ;

  Vector3D p0 (0. , 0. , 10. ) ;
  Vector3D p1 (75. , 75. , 10. ) ;
  
  Vector3D lineVec ;

  for( double phi=-M_PI ; phi < M_PI ; phi += (M_PI/180) ) {
    lineVec = Vector3D ( sin(phi) , cos(phi) , 0 ) ;
    Vector3D iPoint0 = vxdParams.intersectionLadder( p0 , lineVec ) ;
    Vector3D iPoint1 = vxdParams.intersectionLadder( p1 , lineVec ) ;
    vxdIntersectionMap->fill( 0, iPoint0[0] ) ;
    vxdIntersectionMap->fill( 1, iPoint0[1] ) ;
    vxdIntersectionMap->fill( 2, iPoint0[2] ) ;
    vxdIntersectionMap->fill( 3, false ) ;
    vxdIntersectionMap->addRow() ;  
    vxdIntersectionMap->fill( 0, iPoint1[0] ) ;
    vxdIntersectionMap->fill( 1, iPoint1[1] ) ;
    vxdIntersectionMap->fill( 2, iPoint1[2] ) ;
    vxdIntersectionMap->fill( 3, false ) ;
    vxdIntersectionMap->addRow() ;
  }

  mytree->commit() ;
  
  std::cout << "\ntotal tested points: " << nPoints << "\n" 
	    << "          in ladder: " << nLadder << " (" << nLadder/nPoints*100 << ") \n"
	    << "       in sensitive: " << nSensitive << " (" << nSensitive/nPoints*100 << ")\n"
	    << "Points finished " << std::endl ;
  

  // start visualizing
  
  

}

#endif


void testVXDPoint( const Vector3D p , const VXDParameters& vxdParams ) {
  
  std::cout <<  p[0] <<" " <<  p[1] <<" "<<  p[2] ;
  bool isPoint = vxdParams.isPointInLadder( p ) ;
  
  if( isPoint ) {
    std::cout <<  p[0] <<" " <<  p[1] <<" "<<  p[2] << " inside";
  }
  
  std::cout << std::endl ;
  
}

void testVXDDist( const Vector3D p , const VXDParameters& vxdParams ) {
  
  Vector3D v = vxdParams.distanceToNearestLadder( p ) ;

  std::cout << "\nPoint  (" << p[0] << "," << p[1] << "," << p[2] << ") "
	    << "\ndistance Vector " << v << std::endl;
// 	    << "\ndistance Vector (" << v[0] << "," << v[1] << "," << v[2] << ") " << std::endl;


  Vector3D paL = p +  v ;

//   //FG: test the point found:
//   std::cout << " is point in ladder : " << vxdParams.isPointInLadder( paL )  << std::endl ;
  
  // how far are we still from the ladder ?
  Vector3D v2 = vxdParams.distanceToNearestLadder( paL ) ;

  std::cout << " distance to ladder : " << v2.r()  << std::endl ;
  



//   bool isFound = false ;

//   double myStep =1 ;
//   double lBound = 1 ;
//   double uBound = 1 ;

//   for(double r1 = lBound ; r1 <= uBound ; r1+= myStep ) {
//     for(double r2 = lBound ; r2 <= uBound ; r2+= myStep ) {
//       for(double r3 = lBound ; r3 <= uBound ; r3+= myStep ) {

// 	Vector3D iP ( p[0] + r1* v[0] , p[1] + r2* v[1] , p[2] + r3*v[2] ) ;

// 	bool isPoint = vxdParams.isPointInLadder( iP ) ;
// 	if( isPoint ) {
// 	  std::cout << "Vector point withhin Ladder at r =" << r1 <<" , " << r2 <<" , " << r3 << std::endl ;
// 	  std::cout << "                                  " << iP[0] << " , " << iP[1] << " , " << iP[2] << std::endl ;
// 	  isFound = true ;
// 	  break ;
// 	}
// 	else{
// 	  std::cout << "ERROR Vector NOT pointin in Ladder." << std::endl ;
// 	}
//       }
//       if (isFound ) break ;
//     }
//     if (isFound ) break ;
//   }


}

void testAllVXDPoint( const VXDParameters& vxdParams ) {

  double testStart = -65 ;
  double testEnd = 65 ;
  double testStep = 0.05 ;
  double zStart = 10 ;
  double zEnd = 11 ;
  double zStep = 1 ;

  float myStatusP = 0, myStatusL = 0;
  int nLadder = 0 , nSensitive = 0 ;

  std::ofstream myfile;
  const char *fileName = "./gearVXDTest.txt" ;
  myfile.open (fileName, std::ios::out | std::ios::trunc );

  
#ifdef GEAR_USE_AIDA  
// IHistogram2D *ladderMap = myhistofactory->createHistogram2D("ladderMap",
// 							   500, testStart, testEnd,
//							   500, testStart, testEnd );
  std::string mapName = "vxdMap" ;
  std::string mapTitle = "Map of the VTX" ;
  std::vector<std::string> columnNames ;
  std::vector<std::string> columnTypes ;

  columnNames.push_back( "x" )         ; columnTypes.push_back( "double" ) ; 
  columnNames.push_back( "y" )         ; columnTypes.push_back( "double" ) ; 
  columnNames.push_back( "z" )         ; columnTypes.push_back( "double" ) ; 
  columnNames.push_back( "sensitive" ) ; columnTypes.push_back( "bool" ) ;

  ITuple *ladderMap = mytuplefactory->create( mapName , mapTitle , columnNames , columnTypes ) ;

  
#endif

  int nPoints = int( (testEnd-testStart)*(testEnd-testStart)*(zEnd-zStart) 
		     / ( testStep * testStep * zStep ) ) ;
  std::cout << "Testing " << nPoints << " Points to AIDA" << std::endl ;
  
  // x
  for ( double x=testStart ; x < testEnd ; x += testStep ) {
    
    // Status
    float curStatus = (x - testStart) / (testEnd - testStart) * 100 ;
    if ( curStatus > ( myStatusP + 0.25 ) ) {
      myStatusP = curStatus ;
      std::cout << "." << std::flush ;
    }
    if ( curStatus > ( myStatusL + 10 ) ) {
      myStatusL = curStatus ;
      std::cout << floor(myStatusL) << "% finished" << std::endl ;
    }
	 
    // y
    for ( double y=testStart ; y < testEnd ; y += testStep ) {
      
      // z
      for ( double z = zStart ; z < zEnd ; z += zStep ) {
	
	Vector3D p ( x , y , z ) ;
	bool isPoint = vxdParams.isPointInLadder( p ) ;
	bool isSensitive = vxdParams.isPointInSensitive( p ) ;

	if( isPoint ) nLadder ++;
	if( isSensitive ) nSensitive ++;

	if( isPoint || isSensitive ) {
	  myfile <<  p[0] <<" " <<  p[1] <<" "<<  p[2] << " " << isSensitive << "\n" ;

#ifdef GEAR_USE_AIDA  
	  ladderMap->fill( 0 , p[0] ) ;
	  ladderMap->fill( 1 , p[1] ) ;
	  ladderMap->fill( 2 , p[2] ) ;
	  ladderMap->fill( 3 , isSensitive ) ;
	  ladderMap->addRow() ;
#endif
	}
  
      } // z
    } // y
  } //x
  
  myfile.close() ;
  
  std::cout << "\ntotal tested points: " << nPoints << "\n" 
	    << "          in ladder: " << nLadder << " (" << nLadder/nPoints*100 << ") \n"
	    << "       in sensitive: " << nSensitive << " (" << nSensitive/nPoints*100 << ")\n"
	    << "Points finished " << std::endl ;

}


void testAllVXDDist( const VXDParameters& vxdParams ) {
  
  double testStart = -65 ;
  double testEnd = 65 ;
  double testStep = 1 ;
  double zStart = -120 ;
  double zEnd = 120 ;
  double zStep = 10 ;
    
  float myStatusP = 0, myStatusL = 0;
  
  std::ofstream myfile, myHisto;
  const char *fileName     = "./gearVXDTest_out.txt" ;
  const char *histoName  = "./gearVXDTest_Histogram.txt" ;
  myfile.open (fileName, std::ios::out | std::ios::trunc );
  myHisto.open( histoName, std::ios::out | std::ios::trunc ) ;
  
  int nPoints = int ( (testEnd-testStart)*(testEnd-testStart)*(zEnd-zStart) 
		      / ( testStep * testStep * zStep ) ) ;
  std::cout << "Testing " << nPoints << " Points to file : " << fileName <<std::endl ;
  
  // x
  for ( double x=testStart ; x < testEnd ; x += testStep ) {
    
    // Status
    float curStatus = (x - testStart) / (testEnd - testStart) * 100 ;
    if ( curStatus > ( myStatusP + 0.25 ) ) {
      myStatusP = curStatus ;
      std::cout << "." << std::flush ;
    }
    if ( curStatus > ( myStatusL + 10 ) ) {
	myStatusL = curStatus ;
	std::cout << floor(myStatusL) << "% finished" << std::endl ;
    }
    
    // y
    for ( double y=testStart ; y < testEnd ; y += testStep ) {
      
      // z
      for ( double z = zStart ; z < zEnd ; z += zStep ) {
	
	Vector3D p ( x , y , z ) ;
	Vector3D nv = vxdParams.distanceToNearestLadder( p ) ;
	Vector3D np ( p[0] + nv[0] , p[1] + nv[1] , p[2] + nv[2] ) ;
	  
	//bool isPoint = vxdParams.isPointInLadder( np ) ;
	
	//if( !isPoint ) {
	  Vector3D nv2 = vxdParams.distanceToNearestLadder( np ) ;
	  myHisto << sqrt( nv2[0]*nv2[0] + nv2[1]*nv2[1] + nv2[2]*nv2[2] ) << "\n" ;
	  myfile << np[0] <<" " <<  np[1] <<" "<<  np[2] << 0 << "\n" ;
	  //}
	
      } // z
    } // y
  } //x
  
  myfile.close() ;
  myHisto.close() ;
    
  std::cout << "Points finished " << std::endl ;
}
