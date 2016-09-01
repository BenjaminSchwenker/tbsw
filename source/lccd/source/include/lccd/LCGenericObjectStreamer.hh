#ifndef LCGenericObjectStreamer_h
#define LCGenericObjectStreamer_h

#ifdef USE_CONDDB

#include "lccd.h"
#include "VCollectionStreamer.hh"

// -- lcio headers
#include "lcio.h"
#include "EVENT/LCCollection.h"


namespace lccd {
  

  /** Base class for LCCollection streamer classes.
   * 
   *  @author F.Gaede, DESY
   *  @version $Id: LCGenericObjectStreamer.hh,v 1.2 2007-04-27 13:21:31 gaede Exp $
   */

  class LCGenericObjectStreamer : public VCollectionStreamer {
    
  public:
    
    LCGenericObjectStreamer() ;

    virtual ~LCGenericObjectStreamer() {} ;


    /** Writes collection to stream ( string) using XDR. */
    virtual void storeToDB(std::string &s) const ;


    /** Reads collection from stream ( string) using XDR. */
    virtual void retrieveFromDB(const std::string &s) ;


    /** Covariant factory method. Subclasses need to return 
     *  a new instance of their type.
     */
    virtual LCGenericObjectStreamer* create()  ;


  protected:


  };

} //end namespace

#endif // USE_CONDDB
#endif // LCGenericObjectStreamer_h
