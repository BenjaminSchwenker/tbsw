#ifndef StreamerMgr_h
#define StreamerMgr_h

#ifdef USE_CONDDB

// 
#include "lccd.h"
#include "VCollectionStreamer.hh"

// -- lcio headers
#include "lcio.h"
#include "EVENT/LCCollection.h"

// -- CondDB headers
#include "ConditionsDB/ICondDBStreamable.h"

// -- C++ headers
#include <map>

namespace lccd {
  
  typedef  std::map< std::string, VCollectionStreamer* > StreamerMap ;


  /** Manager for streamer classes of LCCollections.<br>
   *  Singleton.
   * 
   *  @author F.Gaede, DESY
   *  @version $Id: StreamerMgr.hh,v 1.2 2007-04-27 13:21:31 gaede Exp $
   */

  class StreamerMgr {
    
  public:
    
    StreamerMgr() ;

    virtual ~StreamerMgr() {} ;

    /** The instance of this manager class.
     */
    static StreamerMgr* instance() ;


    /** Factory method that returns a new streamer object for the given LCIO type, 
     *  e.g. LCIO::LCGenericObject.<br>
     *  Null if type is not known/registered.
     */
    VCollectionStreamer* getStreamer( const std::string&  type ) ; 
    

    /** Register a streamer with the given type name. 
     *  Overwrites existing streamers for that type.
     */
    void registerStreamer( VCollectionStreamer* streamer , const std::string&  type ) ;


  protected:


  private:
    StreamerMap _map ;

    static StreamerMgr* _me ;

  };

}  //end namespace

#endif // USE_CONDDB
#endif // StreamerMgr_h
