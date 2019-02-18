#ifndef CellIDDecoder_h
#define CellIDDecoder_h 1

//#include <vector>
//#include "EVENT/LCObject.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimTrackerHit.h"
#include "UTIL/BitField64.h"
#include "lcio.h"
#include <string>

// fixes problem in gcc 4.0.3
#include "EVENT/LCParameters.h"

//#include <sstream>
//#include <typeinfo>

using namespace lcio ;

namespace UTIL{


  /** Convenient class for decoding cellIDs from collection parameter LCIO::CellIDEncoding.
   *  See UTIL::BitField64 for a description of the encoding string. 
   * 
   *  @see BitField64
   *  @version $Id: CellIDDecoder.h,v 1.9 2006/12/08 10:18:59 gaede Exp $
   */
  template <class T> 
  class CellIDDecoder {
    
  public:  
    
    /** Constructor reads encoding string from collection parameter LCIO::CellIDEncoding.
     * optional BitField64* _defaultBitField can be given which stores the constructed Encoding.
     * This allows to construct a CellIDDecoder without constructing the complete encoding everytime.
     */
    CellIDDecoder( const LCCollection* col,  BitField64* _defaultBitField=nullptr ) : _oldHit(0) {

      std::string initString("") ; 

      if( col !=0 ) 
	initString = col->getParameters().getStringVal(  LCIO::CellIDEncoding ) ;
      
      if( initString.size() == 0 ) {
	
	initString = *_defaultEncoding ;

	std::cout << "    ----------------------------------------- " << std::endl  
		  << "       WARNING: CellIDDecoder - no CellIDEncoding parameter in collection ! " 
		  << std::endl 
		  << "         -> using default : \"" << initString << "\"" 
		  << std::endl 
		  << "    ------------------------------------------ "  
		  << std::endl ;
      }
      
      if(_defaultBitField==nullptr){
          _b = new BitField64(  initString ) ;
          _owningBitField=true;
      } else {
          _b=_defaultBitField;
          _b->reinitialize(initString);  //check if encoding is correct - else create encoding
          _owningBitField=false;
      }

    }
    
    ~CellIDDecoder(){ 
      if(_owningBitField) delete _b ;
    } 
    
    
    /** Provides access to the bit fields, e.g. <br>
     *   int layer =  myCellIDEncoding( hit )[ "layer" ] ;
     * 
     */
    inline const BitField64 & operator()( T* hit ){  
      
      if( hit != _oldHit && hit ) {
	

	long64 val = long64( hit->getCellID0() & 0xffffffff ) 
 	  |        ( long64( hit->getCellID1() ) << 32      ) ;
	
	_b->setValue( val ) ;

	_oldHit = hit ;
      }
      
      return  *_b ;
    }
    

    /** This can be used to set the default encoding that is used if no
     *  CellIDEncoding parameter is set in the collection, e.g. in older lcio files.
     */ 
    static void setDefaultEncoding(const std::string& defaultEncoding ) {
      
      delete _defaultEncoding ;
      
      _defaultEncoding = new std::string( defaultEncoding ) ;
    }
    
  protected:
    BitField64* _b ;

    T* _oldHit ;
    bool _owningBitField;
    
    static std::string*  _defaultEncoding ;
  } ; 
  
  template <class T>
  std::string* CellIDDecoder<T>::_defaultEncoding 
  = new std::string("byte0:8,byte1:8,byte2:8,byte3:8,byte4:8,byte5:8,byte6:8,byte7:8") ;

  

  /** Provides access to the bit fields, e.g. <br>
   *   int layer =  myCellIDEncoding( hit )[ "layer" ] ;
   * Specialization for SimTrackerHits that have only one cellID.
   */
  template<>
  inline const BitField64 & CellIDDecoder<SimTrackerHit>::operator()( SimTrackerHit* hit ){  
    
    if( hit != _oldHit && hit ) {
      
      long64 val = long64( hit->getCellID() & 0xffffffff )  ;
      
      _b->setValue( val ) ;
      
      _oldHit = hit ;
    }
    
    return  *_b ;
  }
  
  
  template <>
  std::string* CellIDDecoder<SimTrackerHit>::_defaultEncoding  ;

} // namespace
#endif


