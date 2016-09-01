#ifndef Xdr_h
#define Xdr_h

#ifdef USE_CONDDB

// 
#include "lccd.h"

#include <typeinfo>

namespace lccd {
  
  /** Utility functions for DB streamer objects. 
   *  These functions take care of the byte ordering and
   *  make sure the streamed data adheres to the xdr convention, e.g.
   *  every data item starts at a 32-bit word boundary.
   * 
   *  @author F.Gaede, DESY
   *  @version $Id: Xdr.hh,v 1.3 2005-09-09 08:09:37 gaede Exp $
   */
  
  class Xdr {
    
  public:
    
    /** Copies the the data in XDR, i.e. machine independent format 
     *  from one location to the another
     */
    static void copy( unsigned const char* from, unsigned char* dest, 
		      const int size, const int count);
    
    
    
    /** Writes count objects of (simple) type T to the string s
     *  in Xdr format.
     */
    template <class T>
    static void tostream( std::string& s , T* t, unsigned count ) {
      
      int bytlen, padlen;
      bytlen = count * sizeof( T ) ;
      padlen = (bytlen + 3) & 0xfffffffc ; 
      
      
      unsigned index = s.size() ;
      
      s.resize( s.size() + padlen ) ;
      
      unsigned char* dest = reinterpret_cast< unsigned char* > (  & s[ index ] ) ;
      unsigned char* from = reinterpret_cast< unsigned char* > ( t )  ;

//        std::cout << " tostream: 0x" << std::hex << *t  << std::dec << std::endl ; 
      
      Xdr::copy( from , dest , sizeof( T ) , count  ) ;
      
    }

    /** Read count objects of (simple) type t from the string/stream s at position.<br>
     *  Callers have to keep the current position within the string.
     * 
     *  @return the new position in the stream.
     */ 
    template <class T>
    static unsigned fromstream( const std::string& s , unsigned position, T* t, unsigned count ) {
      
      if( count == 0 )
	return position ;

      unsigned bytlen, padlen;
      bytlen = count * sizeof( T ) ;
      padlen = (bytlen + 3) & 0xfffffffc ; 
      

      unsigned const char* from = reinterpret_cast< unsigned const char* > (  & s[ position ] ) ;
      unsigned char* dest = reinterpret_cast< unsigned char* > ( t )  ;
      


      Xdr::copy( from , dest , sizeof( T ) , count  ) ;

     
//       std::cout << " fromstream: " ; 
//       for( unsigned int i=0; i < count; i++ )
//  	std::cout	<< std::hex << t[i]  << ", " << std::dec  ;
//       std::cout	<< "   - [" << position << "] - type:" 
//  		<< typeid( *t ).name()  
//  		<< std::endl ; 
      
      
      return position + padlen ;
    }
    
    
    
    /** Specialization for strings - count is ignored and taken to be 1.
     */
//     template<>
    static void tostream( std::string& s , std::string* strVal, unsigned count ) {
      
      int  strSize = (int) strVal->size() ;
      
      int bytlen, padlen;
      bytlen = strSize + sizeof( int )  ;
      
      padlen = (bytlen + 3) & 0xfffffffc ; 
      
      unsigned index = s.size() ;
      
      s.resize( index + padlen ) ;
      
      
      // write the string's size first 
      
      unsigned char* dest = reinterpret_cast< unsigned char* > (  & s[ index ] ) ;
      unsigned char* from = reinterpret_cast< unsigned char* > ( & strSize )  ;
      
      Xdr::copy( from , dest , sizeof( int ) , 1  ) ;
      
      
      dest = reinterpret_cast< unsigned char* > (  & s[ index + sizeof( int ) ] ) ;
      from = reinterpret_cast< unsigned char* > ( & strVal->operator[](0) )  ;
      
      //     std::cout << " tostream: - " << *strVal 
      // 	      << " - lenght: " << strSize
      // 	      << std::endl ;
      
      Xdr::copy( from , dest , 1 , strSize  ) ;
      
    }

  
    /** Specialization for strings - count is ignored and taken to be 1.
     *  @return the new position in the stream.
     */ 
//     template <>
    static unsigned  fromstream( const std::string& s, unsigned index, 
				 std::string* strVal, unsigned count ) {
      
      if( count == 0 )
	return index ;
      
      
      // read the string length first 
      int strSize ;
      unsigned const char* from = reinterpret_cast< unsigned const char* > (  & s[ index ] ) ;
      unsigned char* dest = reinterpret_cast< unsigned char* > ( & strSize )  ;
      
      Xdr::copy( from , dest , sizeof( int ) , 1  ) ;
      
      index += sizeof( int ) ;
      
      
      unsigned bytlen, padlen;
      bytlen =  strSize ;
      padlen = (bytlen + 3) & 0xfffffffc ; 
      
      
      strVal->resize( strSize ) ;
      
      from = reinterpret_cast< unsigned const char* > (  & s[ index ] ) ;
      dest = reinterpret_cast< unsigned char* > ( & strVal->operator[](0)  )  ;
      
      Xdr::copy( from , dest , 1 , strSize  ) ;
      
      //     std::cout << " fromstream: - " << *strVal 
      // 	      << " - lenght: " << strSize
      // 	      << " - [" << index + padlen << "] " 
      // 	      << std::endl ;
      
      
      return index + padlen ;
    }
    

  };  // class Xdr
}  //end namespace

#endif // USE_CONDDB
#endif // Xdr_h
