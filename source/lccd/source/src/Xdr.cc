#ifdef USE_CONDDB
#include"lccd/Xdr.hh"

#include <iostream>

namespace lccd{

  //------------------------------------------
  // Code copied from SIO_functions.
  //------------------------------------------

#if defined(__alpha__) ||   defined(__i386__)  ||   defined(_M_ALPHA)  ||   defined(_M_IX86) || defined(__x86_64__)
#define XDR_LITTLE_ENDIAN
#endif
  
#if defined(_AIX)      ||   defined(__sparc__) ||   defined(_M_PPC)
#define XDR_BIG_ENDIAN
#endif
  
#ifdef XDR_BIG_ENDIAN
  void Xdr::copy( unsigned const char*     from,
		  unsigned char*     dest,
		  const int          size,
		  const int          count ) {
    int totlen;

    totlen = size * count;

    // Keep the byte ordering during the copy.
    memcpy( dest, from, totlen );
  }
#endif

#ifdef XDR_LITTLE_ENDIAN
  void Xdr::copy( unsigned const char*     from,
		  unsigned char*     dest,
		  const int          size,
		  const int          count ) {

    int icnt, ibyt, jump;
    dest += size;
    jump  = size << 1;

    //  Reverse the byte ordering during the copy.
    for( icnt = 0; icnt < count; icnt++ ) {

      for( ibyt = 0; ibyt < size; ibyt++  ) {

        *--dest = *from++;
      }

      dest += jump;
    }
    
  }
#endif


}
#endif
