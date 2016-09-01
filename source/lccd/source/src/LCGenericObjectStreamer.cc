#ifdef USE_CONDDB
#include "lccd/LCGenericObjectStreamer.hh"
#include "lccd/Xdr.hh"


#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCGenericObjectImpl.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/LCParametersImpl.h"
#include "UTIL/LCTOOLS.h"

#include <vector>

using namespace lcio ;

namespace lccd {


  LCGenericObjectStreamer::LCGenericObjectStreamer() {

//     std::cout << " LCGenericObjectStreamer::LCGenericObjectStreamer() : " 
// 	      << this << std::endl ; 
 

  }
  
  void LCGenericObjectStreamer::storeToDB(std::string &s) const {
    
    
//     std::cout << " LCGenericObjectStreamer::storeToDB() : " 
// 	      << this 
// 	      << " _col : " << _col 
// 	      << std::endl ; 

    bool isFixedSize = true ;

    unsigned int nObj = _col->getNumberOfElements() ;
    
    // create a vector with pointers
    std::vector<LCGenericObject*> objects ;
    objects.resize( nObj ) ;
    for( unsigned int i=0 ; i < nObj ; i++ ){
      objects[i] =  dynamic_cast< LCGenericObject* >( _col->getElementAt( i ) ) ; 
    }
    
    // need to check whether we have fixed size objects only
    for( unsigned int i=0 ; i < nObj ; i++ ){
      if( !  objects[i]->isFixedSize() ){ 
	isFixedSize = false ;
	break ;
      }
      
    }
    LCFlagImpl flag( _col->getFlag() ) ;
    
    if(  _col->parameters().getStringVal( "TypeName" ).size() ==  0 ){
      if (nObj>0) 
	_col->parameters().setValue( "TypeName", objects[0]->getTypeName() ) ;
      else
	_col->parameters().setValue( "TypeName", "LCGenericObject" ) ;
    }
    if( isFixedSize ) {
      flag.setBit( LCIO::GOBIT_FIXED ) ;

      if(  _col->parameters().getStringVal( "DataDescription" ).size() ==  0 ){
	if (nObj>0) 
	  _col->parameters().setValue( "DataDescription", objects[0]->getDataDescription() ) ;
	else
	  _col->parameters().setValue( "DataDescription", "UNKNOWN" ) ;
      }
    }


    // make sure that the transient flag is not set when collections are stored 
    flag.unsetBit( LCCollection::BITTransient ) ;

    unsigned flagWord = flag.getFlag() ;

    Xdr::tostream( s , &flagWord , 1 ) ;
    
    // write the collection parameters (defined in VCollectionStreamer)
    writeParameters( s ) ;
    
    Xdr::tostream( s , &nObj , 1 ) ;
    
    
    for( unsigned i=0 ; i< nObj ; i++ ){
      
      unsigned nInt, nFloat, nDouble ;
      
      if( ! isFixedSize || ( isFixedSize && (i==0) )   ) { 
	
	
	nInt = objects[i]->getNInt() ;
	nFloat = objects[i]->getNFloat() ;
	nDouble = objects[i]->getNDouble() ;
	
	Xdr::tostream( s , &nInt , 1 ) ;
	Xdr::tostream( s , &nFloat , 1 ) ;
	Xdr::tostream( s , &nDouble , 1 ) ;
      } 
      
      for(unsigned  int j=0 ; j< nInt ; j++){
	int iv =  objects[i]->getIntVal( j ) ;
	Xdr::tostream( s , &iv , 1 ) ;
      }
      
      for(unsigned int j=0 ; j< nFloat ; j++){
	float fv =  objects[i]->getFloatVal( j ) ;
	Xdr::tostream( s , &fv , 1 ) ;
      }
      for( unsigned int j=0 ; j< nDouble ; j++){
	double dv =  objects[i]->getDoubleVal( j ) ;
	Xdr::tostream( s , &dv , 1 ) ;
      }
      
    }
  }
  

  void LCGenericObjectStreamer::retrieveFromDB(const std::string &s) {

    _col = new LCCollectionVec( LCIO::LCGENERICOBJECT ) ;

    // current position within the string/stream
    unsigned pos(0) ;

    int flagWord ;
    pos = Xdr::fromstream( s , pos , &flagWord , 1 ) ;
    
    _col->setFlag( flagWord ) ;
    
    // read the paramers for this collection  (defined in VCollectionStreamer)
    pos = readParameters( s , pos ) ;

    bool isFixedSize =  LCFlagImpl(flagWord).bitSet( LCIO::GOBIT_FIXED ) ;
    
    
    unsigned int nObj ;

    pos = Xdr::fromstream( s , pos , &nObj , 1 ) ;
    

//     std::cout << " LCGenericObjectStreamer::retrieveFromDB() : " 
// 	      << this 
// 	      << " _col : " << _col 
//       	      << " flag : 0x" << std::hex << flagWord << std::dec 
// 	      << " nObj : " << nObj 
// 	      << std::endl ; 
    

    for( unsigned int i=0 ; i< nObj ; i++ ){
      
      LCGenericObjectImpl* obj ;

      static unsigned nInt, nFloat, nDouble ;

      // create temporary vectors to read the data from the stream
      std::vector<int> intVec ;
      std::vector<float> floatVec ;
      std::vector<double> doubleVec ;
      
      if( ! isFixedSize || ( isFixedSize && (i==0) )   ) { 
       	
	pos = Xdr::fromstream( s , pos , &nInt , 1 ) ;
	pos = Xdr::fromstream( s , pos , &nFloat , 1 ) ;
	pos = Xdr::fromstream( s , pos , &nDouble , 1 ) ;
      } 
      

      // FIXME: need to review LCGenericObjectImpl ...
      // for now the c'tor is the only way to set isFixedSize ...
      if( isFixedSize ) 
	obj = new LCGenericObjectImpl( nInt, nFloat, nDouble )   ;
      else
	obj = new LCGenericObjectImpl ;
	

      // read data into temp vecs
      intVec.resize( nInt ) ;
      floatVec.resize( nFloat ) ;
      doubleVec.resize( nDouble ) ;

      pos = Xdr::fromstream( s , pos , &intVec[0] , nInt ) ;
      pos = Xdr::fromstream( s , pos , &floatVec[0] , nFloat ) ;
      pos = Xdr::fromstream( s , pos , &doubleVec[0] , nDouble ) ;


      // copy to generic object 
      for(unsigned j = nInt-1 ; (j+1) > 0 ; j-- ){
	obj->setIntVal( j ,  intVec[j] ) ;
      }
      for(unsigned j = nFloat-1 ; (j+1) > 0 ; j-- ){
	obj->setFloatVal( j ,  floatVec[j] ) ;
      }
      for(unsigned j = nDouble-1 ; (j+1) > 0 ; j-- ){
	obj->setDoubleVal( j ,  doubleVec[j] ) ;
      }
      
      _col->addElement( obj ) ;
 
    }

    
  }
  
  
  
  LCGenericObjectStreamer* LCGenericObjectStreamer::create()  {
    return new LCGenericObjectStreamer ;
  }
  
  
}
#endif
