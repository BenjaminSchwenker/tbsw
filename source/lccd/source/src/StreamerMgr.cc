#ifdef USE_CONDDB

#include "lcio.h"
#include "EVENT/LCIO.h"

#include "lccd/StreamerMgr.hh"
#include "lccd/LCGenericObjectStreamer.hh"

using namespace lcio ;

namespace lccd{

  StreamerMgr* StreamerMgr::_me  = 0 ;
  

  StreamerMgr* StreamerMgr::instance(){
    if( _me == 0) 
      _me = new StreamerMgr ;
    return _me ;
  }
  
  StreamerMgr::StreamerMgr(){

    // provide some default streamers
    _map[ LCIO::LCGENERICOBJECT ] = new lccd::LCGenericObjectStreamer  ;
  }
  

  VCollectionStreamer* StreamerMgr::getStreamer( const std::string&  type ) {

    VCollectionStreamer* st ;
    if ( ( st = _map[ type ] ) != 0 )
      return st->create()  ;
    else 
      return 0 ;  //FIXME: or should we throw an Exception ?
  }
  
  
  void StreamerMgr::registerStreamer( VCollectionStreamer* streamer , const std::string&  type ) {
    

    StreamerMap::iterator it =  _map.find( type ) ;

    if( it != _map.end() ) {
      delete it->second ;
      _map.erase( it )  ;
    }

    _map[ type ] = streamer ;
  }
  
  
}
#endif // USE_CONDDB

