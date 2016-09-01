#ifndef DBFileHandler_h
#define DBFileHandler_h

//#include "lccd.h"
#include "ConditionsHandlerBase.hh"

using namespace lcio ;

namespace lccd {
  

  // some convenient typedefs
  typedef std::pair<LCCDTimeStamp,LCCDTimeStamp> ValidityInterval ;
  typedef std::vector< ValidityInterval  > ValidityVector ;
  

  /** Implementation of ConditionsHandlerBase that reads the conditions data 
   *  from a special LCIO file. The file has to be extracted from a data base 
   *  beforehand. e.g. with lccd::DBInterface::createDBFile() 
   *  or created otherwise by fullfilling the following requirements:<br>
   *  <ul>
   *  <li>the run header has a pair of string vector parameters DBSince and DBTill with time stamps
   *      for every event in the file</li>
   *  <li>the collections in the events have to have a common name</li>
   *  <li>consecutive events hold condition collections for consecutive validity time intervalls</li>
   *  </lu> <br>
   *  This can be used for example to run a calibration job with a well defined tag of calibration
   *  constants from a database without actually having access to the database, e.g. in a grid environment.
   *   
   *  @author F.Gaede, DESY
   *  @version $Id: DBFileHandler.hh,v 1.3 2007-04-27 13:21:31 gaede Exp $
   */
  
  class DBFileHandler : public ConditionsHandlerBase {
    
  public:
    
    /** Default constructor. You have to specify the following parameters:
     * 
     * @param fileName: the LCIO file that holds the conditions data.
     * @param name:  of this conditions data set - used, e.g. as collection name 
     * in update(LCEvent* evt)
     * @param inputCollection: name of collection in input file - if not given the name
     * of this conditions data set is used
     */
    DBFileHandler(const std::string& fileName, 
		      const std::string& name,
		      const std::string& inputCollection="" ) ;

    
    virtual ~DBFileHandler() ;
    
    /** Reads the event with the conditions data valid for the given time stamp, if:<br>
     *  <ul>
     *  <li>the time stamp lies outside the validity range of the current collection </li>
     *  <li>and there is an event in the file that has a corresponding validity range as defined
     *      by the string vector parameters DBSince and DBTill in the run header</li>
     *  </lu>
     *  
     */
    virtual void update( LCCDTimeStamp timestamp ) ;
    

  protected :
    
    DBFileHandler() {} 

    /** Initializes the validity intervall map from the run header */
    void init() ;

    /** The conditions event number for the given time stamp */
    int findEventNumber( LCCDTimeStamp timeStamp ) ;

    // ----  data members ------

    /** Name of LCIO input file with conditions data */
    std::string _fileName ;

    /** Name of LCCollection in first event */
    std::string _inputCollection ;

    /** Vector that holds the validity time intervalls for all events */
    ValidityVector _valVec ;

    /** The LCIO input file reader */
    LCReader* _lcReader ;

  };
}  //end namespace

#endif // DBFileHandler_h
