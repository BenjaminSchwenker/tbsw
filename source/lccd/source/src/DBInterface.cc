#include "lccd/DBInterface.hh"

#ifdef USE_CONDDB

// lccd headers
#include "lccd.h"
#include "lccd/StreamerMgr.hh"
#include "lccd/DBConnectionMgr.hh"

// LCIO headers
#include "lcio.h"
#include "EVENT/LCIO.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "IO/LCWriter.h"

#include "IMPL/LCRunHeaderImpl.h"
#include "IMPL/LCEventImpl.h"
#include "UTIL/LCTime.h"
#include "UTIL/LCTOOLS.h"

// -- CondDB headers
#include "ConditionsDB/ICondDBMgr.h"
#include "ConditionsDB/CondDBMySQLMgrFactory.h"
#include "ConditionsDB/CondDBObjFactory.h"
#include "CondDBInterface.h"



#include <memory>
#include <algorithm>
#include <functional>
#include <stdlib.h>
#include <iomanip>

#include <iostream>
#include <sstream>

namespace lccd {

  /** Helper class to sort collections of conditions data w.r.t. to their validity time intervall */
  struct less_wrt_validity : public binary_function<lcio::LCCollection*,lcio::LCCollection*,bool>{
    bool operator() (lcio::LCCollection*  c0, lcio::LCCollection* c1) const {
      return ( std::atoll( c0->parameters().getStringVal(lccd::DBSINCE).c_str() ) <
	       std::atoll( c1->parameters().getStringVal(lccd::DBSINCE).c_str() ) ) ;
    }
  } ;

  DBInterface::DBInterface( const std::string& folder, bool update ) : 
    _dbInit( lccd::getDBInitString() ) , // as defined in $LCCD_DB_INIT 
    _folder( folder ) ,
    _update( update ) ,
    _dbName("") {
    
    init() ;
  }
  
  DBInterface::DBInterface(const std::string& dbInit, const std::string& folder, bool update ) : 
    _dbInit( dbInit ) ,
    _folder( folder ),
    _update( update ) ,
    _dbName("") {
    
    init() ;
  }
  
  DBInterface::~DBInterface() {

    DBConnectionMgr::instance()->releaseConnection(  _dbInit , this ) ;   
//     CondDBMySQLMgrFactory::destroyCondDBMgr( _condDBmgr );
 } 
  
  void DBInterface::init() {
    
    //  Database initialization
//     _condDBmgr = CondDBMySQLMgrFactory::createCondDBMgr();
//     _condDBmgr->init( _dbInit );
    
//     if( ! _condDBmgr->isCondDBcreated() )
      
//       _condDBmgr->createCondDB();
	
//     // FIXME: the following calls are useless as ConDBMySQL doesn't use transactions so far...
//     // is this a problem ?
//     // void CondDBInterface::startUpdate() throw(CondDBException) {}
//     // void CondDBInterface::startRead() throw(CondDBException) {}
//     // void CondDBInterface::abort() throw(CondDBException) {}
//     // void CondDBInterface::commit() throw(CondDBException) {}
    
//     _condDBmgr->startRead();
//     _condDBmgr->openDatabase();
//     _condDBmgr->commit();

    //fg: use the DBConnectionMgr to limit the database connections for this job
     _condDBmgr = DBConnectionMgr::instance()->getDBConnection(  _dbInit , this ) ;   
    
    _condDataAccess = _condDBmgr->getCondDBDataAccess() ;
    _condFolderMgr  = _condDBmgr->getCondDBFolderMgr() ; 
    _condTagMgr     = _condDBmgr->getCondDBTagMgr() ;
    

    // we can't get the database name through the abstract interface as this is 
    // a feature typical for MySQL, i.e. CondDBInterface* is used
    CondDBInterface* interface = dynamic_cast<CondDBInterface*>(_condDBmgr) ;
    if( interface != 0 )
      _dbName = interface->getGeneralDBname() ;

    /*    
     * std::cout << "DBInterface::init: connected to database " 
     *	      << _dbName 
     *	      << " using folder: " << _folder
     *        << std::endl;
     */

    if( _update ) {
       _condDBmgr->startUpdate();

      if (! _condFolderMgr->exist( _folder) ) {
	
 	std::cout << "DBInterface::init: folder does not exist, creating folder " << _folder << std::endl;
	
 	_condFolderMgr->createCondDBFolder( _folder, "", "", true); 
	
  	_condDBmgr->commit();
      }
//      else
//        std::cout << "CondDBFolder already exists...proceding!" << std::endl;
    }    
  }


  void DBInterface::storeCollection( LCCDTimeStamp since, LCCDTimeStamp till,
				     lcio::LCCollection* col ,
				     const std::string& description ){

    if( ! _update ) 
      throw lcio::Exception(" DBInterface::storeCollection: not in update mode !" ) ;

    lccd::VCollectionStreamer* colStr = 
      lccd::StreamerMgr::instance()->getStreamer( col->getTypeName() )  ;
    
    if( colStr == 0  )
      throw lcio::Exception( " DBInterface::storeCollection: no streamer found for collection of type " 
			     + col->getTypeName() ) ;
    
    // create auto pointer to prevent memory leaks if exceptions are thrown
    std::auto_ptr<lccd::VCollectionStreamer> colStreamer( colStr ) ;
    
    colStreamer->setCollection( col ) ;
    

    // prepend the colleciton type to the description 
    std::string dbDescription( col->getTypeName() ) ;
    dbDescription += ": " ;
    dbDescription += description ;
    
    try{
      
      ICondDBObject* condObject 
	= CondDBObjFactory::createCondDBObject(since, till, *colStreamer,
					       dbDescription  );
      
      condDataAccess()->storeCondDBObject( _folder, condObject );
      
      CondDBObjFactory::destroyCondDBObject(condObject);
      
      std::cout << " DBInterface::storeCollection: Sucessfully stored collection in database " 
		<< std::endl; 
      
      condDBMgr()->commit();
      
    } catch (CondDBException &e) {
      
      // this exception seems to never get thrown - 
      // instead it is 'caught' beforehand by some CondDBMySQL code ...
      // FIXME: should understand this
      throw lcio::Exception("DBInterface::storeCollection" +  std::string(e.getMessage()) ) ;
    }
    
  }

  lcio::LCCollection* DBInterface::findCollection( LCCDTimeStamp timeStamp, 
						   const std::string& tag ) {
    LCCDTimeStamp ignore1, ignore2 ;

    return findCollection( timeStamp, ignore1 ,ignore2 , tag ) ; 
  }


  lcio::LCCollection*  DBInterface::findCollection( LCCDTimeStamp timeStamp, 
						    LCCDTimeStamp& since, LCCDTimeStamp& till,
						    const std::string& tag ) { 
    
    ICondDBObject* condObject = 0;
    
    _condDBmgr->startRead(); // dummy method for now 
    
    CondDBKey point =  timeStamp ;
    
    condDataAccess()->findCondDBObject( condObject, _folder, point , tag );
    
    if( condObject == 0 ){
      
      std::cout << "DBInterface::findCollection: Warning: No object found in database folder " << _folder << " for timestamp : " 
		<< point
		<< std::endl ;
      return 0 ;
    }
    
    since = condObject->validSince() ;
    till = condObject->validTill() ;
    
    lcio::LCCollection*  col = collectionFromCondDBObject( condObject, tag )  ;
    
    CondDBObjFactory::destroyCondDBObject(condObject);
    
    return  col ;
  }


  lcio::LCCollection*  DBInterface::findNextValidCollection( LCCDTimeStamp timeStamp, 
							     LCCDTimeStamp& since, LCCDTimeStamp& till,
							     const std::string& tag ) { 
    
    ICondDBObject* condObject = 0;

    _condDBmgr->startRead(); // dummy method for now 

    CondDBKey searchPoint =  timeStamp ;

    condDataAccess()->findNextValidCondDBObject( condObject, _folder, searchPoint , tag );

    if( condObject == 0 ){
      
      std::cout << "DBInterface::findNextValidCollection: Warning: No further object found in database folder " << _folder << " from timestamp : " 
		<< searchPoint
		<< std::endl ;
      return 0 ;
    }

    since = condObject->validSince() ;
    till = condObject->validTill() ;
    
    lcio::LCCollection*  col = collectionFromCondDBObject( condObject, tag )  ;
    
    CondDBObjFactory::destroyCondDBObject(condObject);
    
    return  col ;

  }


  void DBInterface::findCollections( ColVec& colVec, const std::string& tag ){ 

    _condDBmgr->startRead(); // dummy method for now 

    ICondDBDataIterator* objIter ;

    condDataAccess()->browseObjectsInTag( objIter, _folder, tag ) ;

    ICondDBObject* condObject ; 


    // Note: due to a bug/feature in ConditionsDBMysQL implementation
    // the following runs into an endless loop:
    //  while( (condObject = objIter->current() ) != 0 ) {  
    //    ...
    //     objIter->next() ;
    //  } 

    // This way of looping through the iterator produces a memory leak,
    // as current() and next() instantiate condObjects ...
    //     do {
    //       condObject = objIter->current() ;
    //       ...
    //     } while ( objIter->next() );


    // this seems to be the only proper way of looping through the iterator:
    for( condObject = objIter->current() ; condObject != 0  ; condObject = objIter->next() ){
      
      //       cout << " ---- adding collection:  [" <<  condObject->validSince() 
      // 	   << "," << condObject->validTill() << "]" << endl ;
      
      lcio::LCCollection*  col = collectionFromCondDBObject( condObject, tag )  ;
      
      colVec.push_back( col ) ;
      
      CondDBObjFactory::destroyCondDBObject( condObject );
      
    }
    
    delete objIter ;

    // now we have to sort the array w.r.t. to the since time stamp
    // FIXME: this shouldn't be necessary according to the ConditionsDB API description ... ?
    std::sort(  colVec.begin() , colVec.end() ,   less_wrt_validity()  ) ;
    
  }



  void DBInterface::findCollections( ColVec& colVec, LCCDTimeStamp timeStamp ){

    _condDBmgr->startRead(); // dummy method for now 

    ICondDBDataIterator* objIter ;

    condDataAccess()->browseObjectsAtPoint( objIter, _folder, timeStamp ) ;

    ICondDBObject* condObject ; 

    for( condObject = objIter->current() ; condObject != 0  ; condObject = objIter->next() ){
      
      // FIXME: need to determine proper tag - so far DBTag=UNKNOWN for all collections ...
      lcio::LCCollection*  col = collectionFromCondDBObject( condObject, "UNKNOWN" )  ;
      
      colVec.push_back( col ) ;
      
      CondDBObjFactory::destroyCondDBObject( condObject );
      
    }
    delete objIter ;

  }

  void 	DBInterface::findCollections (ColVec &colVec, LCCDTimeStamp since, LCCDTimeStamp till, const std::string &tag) {
    
    _condDBmgr->startRead(); // dummy method for now 
    
    ICondDBDataIterator* objIter = 0;
    
    condDataAccess()->browseHistory(objIter, _folder, since, till, tag ) ;
    
    if (!objIter)  return; // return if no elements found
    
    ICondDBObject* condObject ; 
    
    
    // Note: due to a bug/feature in ConditionsDBMysQL implementation
    // the following runs into an endless loop:
    //  while( (condObject = objIter->current() ) != 0 ) {  
    //    ...
    //     objIter->next() ;
    //  } 
    
    // This way of looping through the iterator produces a memory leak,
    // as current() and next() instantiate condObjects ...
    //     do {
    //       condObject = objIter->current() ;
    //       ...
    //     } while ( objIter->next() );
    
    
    // this seems to be the only proper way of looping through the iterator:
    for( condObject = objIter->current() ; condObject != 0  ; condObject = objIter->next() ){
      
      //       cout << " ---- adding collection:  [" <<  condObject->validSince() 
      // 	   << "," << condObject->validTill() << "]" << endl ;
      
      lcio::LCCollection*  col = collectionFromCondDBObject( condObject, tag )  ;
      
      colVec.push_back( col ) ;
      
      CondDBObjFactory::destroyCondDBObject( condObject );
      
    }
    
    delete objIter ;
    
    // now we have to sort the array w.r.t. to the since time stamp
    // FIXME: this shouldn't be necessary according to the ConditionsDB API description ... ?
    std::sort(  colVec.begin() , colVec.end() ,   less_wrt_validity()  ) ;
    
  }

  void DBInterface::findAllCollections( ColVec& colVec ){ 
      
    _condDBmgr->startRead(); // dummy method for now 
    
    ICondDBDataIterator* objIter ;
    
    condDataAccess()->browseAllObjects( objIter, _folder) ;
    
    ICondDBObject* condObject ; 
    
    // this seems to be the only proper way of looping through the iterator:
    for( condObject = objIter->current() ; condObject != 0  ; condObject = objIter->next() ){
      
      //       cout << " ---- adding collection:  [" <<  condObject->validSince() 
      // 	   << "," << condObject->validTill() << "]" << endl ;
      
      lcio::LCCollection*  col = collectionFromCondDBObject( condObject ,"UNKNWON")  ;
      
      colVec.push_back( col ) ;
      
      CondDBObjFactory::destroyCondDBObject( condObject );
      
    }
    
    delete objIter ;
      
    // now we have to sort the array w.r.t. to the since time stamp
    // FIXME: this shouldn't be necessary according to the ConditionsDB API description ... ?
    std::sort(  colVec.begin() , colVec.end() ,   less_wrt_validity()  ) ;
    
  }



  lcio::LCCollection* DBInterface::collectionFromCondDBObject( ICondDBObject* condObject, 
							       const std::string& tag ) {
    
    
    // need to decode the collection type from object description: "LCIOTYPE: some description"
    std::string desc ;
    condObject->description( desc ) ;
    
    std::string colType( desc , 0 ,  desc.find_first_of(':') ) ;
    
    
    lccd::VCollectionStreamer* colStr = 
      lccd::StreamerMgr::instance()->getStreamer( colType )  ;
    
    if( colStr == 0  ){
      throw lcio::Exception( " DBInterface::findCollection: no streamer found for collection of type " 
 			     + colType ) ;
    }
    
    // create auto pointer to prevent memory leaks if exceptions are thrown
    std::auto_ptr<lccd::VCollectionStreamer> colStreamer( colStr ) ;
    
    condObject->data( *colStreamer ) ;
    
    LCCDTimeStamp since = condObject->validSince() ;
    LCCDTimeStamp till = condObject->validTill() ;
    
    //---- add some parameters to the collection --------
    lcio::LCCollection* col =  colStreamer->getCollection() ;
    
    
    lcio::StringVec strVec ;
    
    std::stringstream sinceStr ;
    sinceStr << since ;
    strVec.push_back(  sinceStr.str() ) ;
    strVec.push_back( lcio::LCTime( since ).getDateString() ) ;
    col->parameters().setValues( lccd::DBSINCE ,  strVec ) ;
    strVec.clear() ;

    std::stringstream tillStr ;
    tillStr << till ;
    strVec.push_back(  tillStr.str() ) ;
    strVec.push_back( lcio::LCTime( till ).getDateString() ) ;
    col->parameters().setValues( lccd::DBTILL ,  strVec ) ;
    strVec.clear() ;

    
    std::stringstream nowStr ;
    lcio::LCTime now ;
    nowStr << now.timeStamp() ;
    strVec.push_back(  nowStr.str() ) ;
    strVec.push_back( now.getDateString() ) ;
    col->parameters().setValues( lccd::DBQUERYTIME ,  strVec ) ;
    strVec.clear() ;
    
    
    SimpleTime ins ;
    condObject->insertionTime( ins ) ; 

    std::stringstream insertStr ;
    lcio::LCTime insert(  fromSimpleTime( ins )   ) ;
    insertStr << insert.timeStamp() ;
    strVec.push_back(  insertStr.str() ) ;
    strVec.push_back( insert.getDateString() ) ;
    col->parameters().setValues( lccd::DBINSERTIONTIME ,  strVec ) ;
    strVec.clear() ;
    
    
    
    std::string dbTag( tag ) ;
    if( dbTag.size() == 0 ) 
      dbTag = "HEAD" ;
    col->parameters().setValue( lccd::DBTAG ,  dbTag ) ;
    
    
    col->parameters().setValue( lccd::DBFOLDER ,  _folder ) ;
    
    col->parameters().setValue( lccd::DBNAME ,  _dbName ) ;

    col->parameters().setValue( lccd::DBLAYER,   (int) condObject->layer() ) ;
    
    
    //---------------------------------------------------
    
    return col ;
  }
  

  void DBInterface::createDBFile( const std::string& tag ) {


    lcio::LCWriter* wrt = lcio::LCFactory::getInstance()->createLCWriter() ;
    

    // create oputput file name :  conddb_COLNAME_TAG_YYYYMMDD_HHMMSS.slcio
    std::stringstream file ;       
    
    std::string dbTag( tag ) ;
    if( dbTag.size() == 0 )  dbTag = "HEAD" ;
    
    std::string colName( _folder, _folder.find_last_of( '/' ) + 1 ,  _folder.size() ) ;
    
    lcio::LCTime now ;
    
    file << "condDB_" 
	 << colName  << "_"
	 << dbTag<< "_"
	 << std::setfill('0') 
	 << std::setw(4) <<  now.year() 
	 << std::setw(2) <<  now.month() 
	 << std::setw(2) <<  now.day() << "_"
	 << std::setw(2) <<  now.hour() //<< "-"
	 << std::setw(2) <<  now.min() //<< "-"
	 << std::setw(2) <<  now.sec() 
	 << ".slcio" ;

//     wrt->open( file.str()  , lcio::LCIO::WRITE_NEW )  ; 
    wrt->open( file.str() ) ;
    
        
    ColVec colVec ;
    findCollections( colVec, tag ) ; 


    // --- write a run header with some information related to the query
    lcio::LCRunHeaderImpl* rHdr = new lcio::LCRunHeaderImpl ;

    lcio::StringVec sinceVec ;
    lcio::StringVec tillVec ;
    
    // add map with validity time intervalls and events - needed by DBFileHandler
    for( ColVec::iterator it = colVec.begin() ; it != colVec.end() ; it++) {
      
      sinceVec.push_back ( (*it)->parameters().getStringVal( lccd::DBSINCE ) ) ;
      tillVec.push_back ( (*it)->parameters().getStringVal( lccd::DBTILL ) ) ;
    }
    
    rHdr->parameters().setValues(  lccd::DBSINCE , sinceVec ) ;
    rHdr->parameters().setValues(  lccd::DBTILL , tillVec ) ;
    
    // write some additional information to the run header
    rHdr->parameters().setValue( lccd::DBTAG ,  dbTag ) ;
    rHdr->parameters().setValue( lccd::DBFOLDER ,  _folder ) ;
    rHdr->parameters().setValue( lccd::DBNAME ,  _dbName ) ;
    
    
    std::stringstream desc ;
    desc << " Conditions data extracted from database on " << now.getDateString() << " GMT \n" ;
    rHdr->setDescription( desc.str() ) ;
    
    //---------------------------------------------------
    
    wrt->writeRunHeader( rHdr ) ;
    

    // loop over collections....
//     std::string colName("ConditionsData") ; // FIXME

    int evtNum(0) ;
    for( ColVec::iterator it = colVec.begin() ; it != colVec.end() ; it++) {
      
      lcio::LCEventImpl* evt = new lcio::LCEventImpl ;
      evt->setEventNumber( evtNum++ ) ;
      
      evt->addCollection(  *it , colName  ) ;  
      
      lcio::LCTOOLS::dumpEventDetailed( evt  ) ;

      wrt->writeEvent( evt ) ;

      delete evt ; // this deletes the collection as well
  
    }

    wrt->close() ;
    // clean up
    delete wrt ;
    delete rHdr ;
    
  }

  void DBInterface::dump2File( ) {


    lcio::LCWriter* wrt = lcio::LCFactory::getInstance()->createLCWriter() ;
    
    
    // create oputput file name :  conddb_COLNAME_YYYYMMDD_HHMMSS.slcio
    std::stringstream file ;       
      
      
    std::string colName( _folder, _folder.find_last_of( '/' ) + 1 ,  _folder.size() ) ;
      
    lcio::LCTime now ;
    
    file << "condDB_" 
	 << colName  << "_"
	 << std::setfill('0') 
	 << std::setw(4) <<  now.year() 
	 << std::setw(2) <<  now.month() 
	 << std::setw(2) <<  now.day() << "_"
	 << std::setw(2) <<  now.hour() //<< "-"
	 << std::setw(2) <<  now.min() //<< "-"
	 << std::setw(2) <<  now.sec() 
	 << ".slcio" ;
    
    //     wrt->open( file.str()  , lcio::LCIO::WRITE_NEW )  ; 
    wrt->open( file.str() ) ;
    
    
    ColVec colVec ;
    findAllCollections( colVec ) ; 
      

    // --- write a run header with some information related to the query
    lcio::LCRunHeaderImpl* rHdr = new lcio::LCRunHeaderImpl ;
    
    lcio::StringVec sinceVec ;
    lcio::StringVec tillVec ;
    
    // add map with validity time intervalls and events - needed by DBFileHandler
    for( ColVec::iterator it = colVec.begin() ; it != colVec.end() ; it++) {
      
      sinceVec.push_back ( (*it)->parameters().getStringVal( lccd::DBSINCE ) ) ;
      tillVec.push_back ( (*it)->parameters().getStringVal( lccd::DBTILL ) ) ;
    }
    
    rHdr->parameters().setValues(  lccd::DBSINCE , sinceVec ) ;
    rHdr->parameters().setValues(  lccd::DBTILL , tillVec ) ;
    
    // write some additional information to the run header
    //      rHdr->parameters().setValue( lccd::DBTAG ,  dbTag ) ;
    rHdr->parameters().setValue( lccd::DBFOLDER ,  _folder ) ;
    rHdr->parameters().setValue( lccd::DBNAME ,  _dbName ) ;
    
    
    std::stringstream desc ;
    desc << " Conditions data extracted from database on " << now.getDateString() << " GMT \n" ;
    rHdr->setDescription( desc.str() ) ;
    
    //---------------------------------------------------
    
    wrt->writeRunHeader( rHdr ) ;
    
      
    // loop over collections....
    //     std::string colName("ConditionsData") ; // FIXME
    
    int evtNum(0) ;
    for( ColVec::iterator it = colVec.begin() ; it != colVec.end() ; it++) {
      
      lcio::LCEventImpl* evt = new lcio::LCEventImpl ;
      evt->setEventNumber( evtNum++ ) ;
      
      evt->addCollection(  *it , colName  ) ;  
      
      lcio::LCTOOLS::dumpEventDetailed( evt  ) ;
      
      wrt->writeEvent( evt ) ;
      
      delete evt ; // this deletes the collection as well
  
    }

    wrt->close() ;
    // clean up
    delete wrt ;
    delete rHdr ;
      
  }



  void DBInterface::createSimpleFile( LCCDTimeStamp timeStamp, const std::string& tag, bool allLayers ) {
    
    lcio::LCWriter* wrt = lcio::LCFactory::getInstance()->createLCWriter() ;
    
    
    // create oputput file name :  conddb_COLNAME_TAG_TIMESTAMP.slcio
    std::stringstream file ;       
    
    std::string dbTag( tag ) ;
    if( dbTag.size() == 0 )  dbTag = "HEAD" ;
    
    std::string colName( _folder, _folder.find_last_of( '/' ) + 1 ,  _folder.size() ) ;
    
    file << "condDB_" 
	 << colName  << "_"
	 << dbTag<< "_"
	 << std::setfill('0') 
	 << std::setw(19) << timeStamp
	 << ".slcio" ;
    
    //     wrt->open( file.str()  , lcio::LCIO::WRITE_NEW )  ; 
    wrt->open( file.str() ) ;
    
    // --- write a run header with some information related to the query
    lcio::LCRunHeaderImpl* rHdr = new lcio::LCRunHeaderImpl ;
    
    rHdr->parameters().setValue( lccd::DBTAG ,  dbTag ) ;
    rHdr->parameters().setValue( lccd::DBFOLDER ,  _folder ) ;
    rHdr->parameters().setValue( lccd::DBNAME ,  _dbName ) ;
    rHdr->parameters().setValue( lccd::DBNAME ,  _dbName ) ;
    std::stringstream ts ;     
    ts << timeStamp ;
    rHdr->parameters().setValue( lccd::DBTIMESTAMP ,  ts.str() ) ;

    lcio::LCTime now ;
    std::stringstream desc ;
    desc << " Conditions data extracted from database on " << now.getDateString() << " GMT \n" ;

    lcio::LCTime tsDate( timeStamp ) ;
    desc << "   for timestamp: " << std::setfill('0') << std::setw(19) << timeStamp   << "\n"
	 << "              ( " <<  tsDate.getDateString() << ") \n" ;
    
    if( allLayers ) {
      desc << "  Includes all layers for the given time stamp \n" ;
    }	

    rHdr->setDescription( desc.str() ) ;
      
    wrt->writeRunHeader( rHdr ) ;
    //---------------------------------------------------
      
    

    lcio::LCEventImpl* evt = new lcio::LCEventImpl ;

    lcio::LCCollection* col = findCollection( timeStamp, tag ) ; 

    evt->addCollection(  col , colName  ) ;  


    if( allLayers ) {

      ColVec colVec ;
      findCollections( colVec, timeStamp ) ; 

      for( ColVec::iterator it = colVec.begin() ; it != colVec.end() ; it++) {
	
	std::stringstream cn ;
	cn <<  colName <<  "_layer_"  
	   << std::setfill('0') << std::setw(3) 
	   << (*it)->parameters().getIntVal( lccd::DBLAYER )  ;
	
	evt->addCollection(  *it , cn.str()  ) ;  
      }
    }

    //---------------------------------------------------

    wrt->writeEvent( evt ) ;

    wrt->close() ;

    delete evt ; // this deletes the collections as well
    delete wrt ;
    delete rHdr ;


  }


  void DBInterface::tagFolder( const std::string& tag,  const std::string& description, std::string usingTagName) {
   
    if( ! _update ) 
      throw lcio::Exception(" DBInterface::tagFolder: not in update mode !" ) ;

    condTagMgr()->createCondDBTag( tag, description ) ;
    condTagMgr()->tag( _folder , tag , usingTagName ) ;

  }

}
#endif
