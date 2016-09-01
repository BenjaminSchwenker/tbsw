#ifdef USE_CONDDB
#include "lccd/VCollectionStreamer.hh"

#include "lccd/Xdr.hh"

#include "EVENT/LCParameters.h"

using namespace lcio ;


namespace lccd {


  void VCollectionStreamer::writeParameters( std::string& s ) const {
    
    const LCParameters& params = _col->getParameters() ;
    
    StringVec intKeys ;
    int nIntParameters = params.getIntKeys( intKeys ).size() ;
    
    Xdr::tostream( s , &nIntParameters , 1 ) ;

    for(int i=0; i< nIntParameters ; i++ ){

      IntVec intVec ;

      params.getIntVals(  intKeys[i], intVec ) ;

      int nInt  = intVec.size()  ;     // = params.getNInt( intKeys[i] ) ;

      if( nInt > 0 ){ 
	Xdr::tostream( s , &intKeys[i] , 1 ) ;
	Xdr::tostream( s , &nInt, 1 ) ;
	Xdr::tostream( s , &intVec[0], nInt  ) ;
      }
    }

    StringVec floatKeys ;
    int nFloatParameters = params.getFloatKeys( floatKeys ).size() ;
    
    Xdr::tostream( s , &nFloatParameters , 1 ) ;

    for(int i=0; i< nFloatParameters ; i++ ){

      FloatVec floatVec ;

      params.getFloatVals(  floatKeys[i], floatVec ) ;

      int nFloat  = floatVec.size()  ;     // = params.getNFloat( floatKeys[i] ) ;

      if( nFloat > 0 ){ 
	Xdr::tostream( s , &floatKeys[i] , 1 ) ;
	Xdr::tostream( s , &nFloat, 1 ) ;
	Xdr::tostream( s , &floatVec[0], nFloat  ) ;
      }
    }
    
    StringVec stringKeys ;
    int nStringParameters = params.getStringKeys( stringKeys ).size() ;
    
    Xdr::tostream( s , &nStringParameters , 1 ) ;

    for(int i=0; i< nStringParameters ; i++ ){
      
      StringVec stringVec ;
      
      params.getStringVals(  stringKeys[i], stringVec ) ;
      
      int nString  = stringVec.size()  ;     // = params.getNString( stringKeys[i] ) ;
      
      if( nString > 0 ){ 
	Xdr::tostream( s , &stringKeys[i] , 1 ) ;
	Xdr::tostream( s , &nString, 1 ) ;
       	
	for(int j=0; j< nString ; j++ ){
	  Xdr::tostream( s , &stringVec[j] , 1 ) ;
	}
      }
    }


  }  

  unsigned VCollectionStreamer::readParameters( const std::string& s , unsigned index ){
    
    LCParameters& params = _col->parameters() ;

    int nIntParameters ;

    index = Xdr::fromstream( s , index , &nIntParameters , 1 ) ;

    for(int i=0; i< nIntParameters ; i++ ){

      std::string key ;
      index = Xdr::fromstream( s , index , &key , 1 ) ;


      int nInt  ;
      index = Xdr::fromstream( s , index , &nInt , 1 ) ;
      
      IntVec intVec( nInt ) ;
      index = Xdr::fromstream( s , index , &intVec[0] , nInt ) ;
     
      params.setValues( key , intVec ) ;
    }

    int nFloatParameters ;
    index = Xdr::fromstream( s , index , &nFloatParameters , 1 ) ;

    for(int i=0; i< nFloatParameters ; i++ ){
      
      std::string key ;
      index = Xdr::fromstream( s , index , &key , 1 ) ;
      
      int nFloat  ;
      index = Xdr::fromstream( s , index , &nFloat , 1 ) ;
      
      FloatVec floatVec(nFloat) ;
      index = Xdr::fromstream( s , index , &floatVec[0] , nFloat ) ;
      
      
      params.setValues( key , floatVec ) ;
    }

    int nStringParameters ;
    index = Xdr::fromstream( s , index , &nStringParameters , 1 ) ;

    for(int i=0; i< nStringParameters ; i++ ){

      std::string key ;
      index = Xdr::fromstream( s , index , &key , 1 ) ;
      
      int nString  ;
      index = Xdr::fromstream( s , index , &nString , 1 ) ;
      StringVec stringVec(nString) ;

      for(int j=0; j< nString ; j++ ){
	index = Xdr::fromstream( s , index , &stringVec[j] , 1 ) ;
      }

      params.setValues( key , stringVec ) ;
    }

    return index ;
  }


} 
#endif
