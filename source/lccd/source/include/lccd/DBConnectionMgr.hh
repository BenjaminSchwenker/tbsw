#ifndef DBConnectionMgr_h
#define DBConnectionMgr_h

#ifdef USE_CONDDB

// 
#include "lccd.h"
// #include "lccd/DBInterface.hh"

// forward declarations
class ICondDBMgr ;

#include <iostream>
#include <set>
#include <map>


namespace lccd {
  
  class DBInterface ;
  
  
  /** Singleton that manages data base connection on the basis of access strings.
   *
   *  @author F.Gaede, DESY
   *  @version $Id: DBConnectionMgr.hh,v 1.2 2007-04-27 13:21:31 gaede Exp $
   */

  class DBConnectionMgr {
    
  private:
    typedef std::set< DBInterface* > DBSet ;
    
    struct DBRefs{
      DBRefs() : Mgr(0) {} 
      ICondDBMgr* Mgr ;
      DBSet Set ; 
    } ;
    
//     typedef std::map< std::string , std::pair< DBSet,ICondDBMgr*> > DBMap ;
    typedef std::map< std::string , DBRefs > DBMap ;
    
    DBConnectionMgr() ;

  public:

    /** The only way to access the manager object */
    static DBConnectionMgr* instance() ;

    /** Get a connection for the given dbInit - create if it does not exist yet */
    ICondDBMgr* getDBConnection( const std::string& dbInit , DBInterface* dbInterface ) ;
    
    /** Release the connection for the given interface */
    void releaseConnection( const std::string& dbInit, DBInterface* dbInterface ) ;
    
    /** D'tor*/
    ~DBConnectionMgr() ;
    
  private:
    DBMap _map ;

  };

}  //end namespace

#endif // USE_CONDDB
#endif // DBConnectionMgr_h
