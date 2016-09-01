#ifndef DBCondHandler_h
#define DBCondHandler_h

#ifdef USE_CONDDB

#include "lccd.h"
#include "ConditionsHandlerBase.hh"

#include "lcio.h"
#include "Exceptions.h"


namespace lccd {

  class DBInterface ;
  
  /** Implementation of ConditionsHandlerBase that reads the conditions data 
   *  from the database. On every call to update() it checks whether the provided
   *  time stamp exceeds the validity of the current conditions data. 
   *  If so, new constants are retrieved from the data base.
   * 
   *  @author F.Gaede, DESY
   *  @version $Id: DBCondHandler.hh,v 1.2 2007-04-27 13:21:31 gaede Exp $
   */
  
  class DBCondHandler : public ConditionsHandlerBase {
    
  public:
    
    /** Default constructor. You have to specify the following parameters:
     * 
     * @param dbInit: initialization string for the database, e.g. "localhost:lccd_test:calvin:hobbes" 
     * @param folder:  name of the database folder
     * @param name:  of this conditions data set - used, e.g. as collection name 
     * in update(LCEvent* evt)
     * @param tag: tag to be used for queries, HEAD if omitted 
     */
    DBCondHandler(const std::string& dbInit,  const std::string& folder,
		  const std::string& name,  const std::string& tag="" ) ;
    
    
    virtual ~DBCondHandler() ;
    
    
    /** Retrieves the conditions data from the data base if time stamp 
     *  exceeds current validity time range.<br>
     *  Notifies all registered listeners if conditions have changed.
     */
    virtual void update( LCCDTimeStamp timeStamp ) ;
//     virtual void update( LCCDTimeStamp timeStamp )  throw (lcio::Exception)  ;
    
    
  protected :
    
    DBCondHandler() {} 
    
    lccd::DBInterface* _db ;

    // ----  data members ------
    std::string _dbInit ;
    std::string _folder ;
    std::string _tag ;
    

  };
}  //end namespace

#endif  // USE_CONDDB

#endif // DBCondHandler_h
