#ifdef USE_CONDDB
#include "lccd/DBConnectionMgr.hh"
#include "ConditionsDB/ICondDBMgr.h"
#include "ConditionsDB/CondDBMySQLMgrFactory.h"
#include "ConditionsDB/CondDBObjFactory.h"

#include <iostream>


namespace lccd {

//     typedef std::set< ICondDBMgr* > DBSet ;
    
//     typedef std::map< std::string , DBSet > DBMap ;
//        typedef std::map< std::string , std::pair< DBSet,ICondDBMgr*> > DBMap ;

  DBConnectionMgr::DBConnectionMgr(){

  }

  DBConnectionMgr::~DBConnectionMgr(){

    // in case there are still open connections we close them
    for(  DBMap::iterator it = _map.begin() ;
	  it != _map.end() ; ++ it ){

      DBRefs& refs = it->second ;

      if( refs.Mgr != 0 ){
	
	CondDBMySQLMgrFactory::destroyCondDBMgr( refs.Mgr );

	// std::cout << " DBConnectionMgr::~DBConnectionMgr closing leftover db connection " 
	//	  << std::endl ;
	
      }
    }
  }
  
  DBConnectionMgr* DBConnectionMgr::instance(){
    // using a local static to ensure destructor is called at end of job:
    static DBConnectionMgr me ;
    return &me ;
    //     static DBConnectionMgr* me ;
    //     if(me==0) me = new DBConnectionMgr ;
    //     return me ;
  }
  
  ICondDBMgr* DBConnectionMgr::getDBConnection( const std::string& dbInit , DBInterface* dbInterface ){
  
    DBRefs& refs = _map[ dbInit ] ; 

    if( refs.Mgr == 0 ) {
      
      // create and initialize a new db connection:
      
      ICondDBMgr* condDBmgr = CondDBMySQLMgrFactory::createCondDBMgr();
      condDBmgr->init( dbInit );
      
      if( ! condDBmgr->isCondDBcreated() )
	condDBmgr->createCondDB();
      
      condDBmgr->startRead();
      condDBmgr->openDatabase();
      condDBmgr->commit();
      
      refs.Mgr = condDBmgr ;

      // std::cout << " DBConnectionMgr::getDBConnection( " << dbInit << ")" << std::endl ;
    }

    refs.Set.insert( dbInterface ) ; // keep track of users ...

    return refs.Mgr ;

  }
  
  void DBConnectionMgr::releaseConnection( const std::string& dbInit , DBInterface* dbInterface )  {
    
    DBRefs& refs = _map[ dbInit ] ; 
    
    refs.Set.erase( dbInterface ) ;

    if( refs.Set.size() == 0 ) { // no one's using this connection any more ...

      CondDBMySQLMgrFactory::destroyCondDBMgr( refs.Mgr );

      refs.Mgr = 0 ;

      // std::cout << " DBConnectionMgr::releaseConnection( " << dbInit << ")" << std::endl ;


    }
  }
  
  



}
#endif
