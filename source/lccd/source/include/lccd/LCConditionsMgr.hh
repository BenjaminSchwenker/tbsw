#ifndef LCConditionsMgr_h
#define LCConditionsMgr_h

// 
#include "lccd.h"
#include "IConditionsHandler.hh"

// -- lcio headers
// #include "lcio.h"
// #include "EVENT/LCCollection.h"

// -- C++ headers
#include <map>

namespace lccd {
  
  typedef  std::map< std::string, IConditionsHandler* > HandlerMap ;


  /** Manager for IConditionsHandler objects.<br>
   *  Can be used to keep all conditions handler objects in an aplication
   *  in one place and treat them through a common interface, e.g. 
   *  call IConditionsHandler::updateEvent().
   *  
   *  Singleton.
   * 
   *  @author F.Gaede, DESY
   *  @version $Id: LCConditionsMgr.hh,v 1.2 2007-04-27 13:21:31 gaede Exp $
   */

  class LCConditionsMgr {

  protected:

    LCConditionsMgr() ;
    
  public:

    ~LCConditionsMgr() {} ;

    /** The instance of this manager class.
     */
    static LCConditionsMgr* instance() ;


    /** Register a change listener with the conditions handler 'name'.
     *  Throws an exception if no handler with the given name has been
     *  registered before.
     */
    void registerChangeListener( IConditionsChangeListener* cl, const std::string&  name) ;


    /** Register a handler with the given type name. 
     */
    void registerHandler(const std::string&  name, IConditionsHandler* handler ) ;


    /** Remove the handler with the given name.
     */
     void removeHandler(const std::string&  name ) ;


    /** Returns the handler for the given name - Null
     *  if not registered.
     */
    IConditionsHandler * getHandler( const std::string&  name) ;

    /** Call IConditionsHandler::updateEvent() for all registered handlers. 
     */
    void updateEvent( lcio::LCEvent* evt ) ;

    /** Call IConditionsHandler::update() for all registered handlers. 
     */
    void update( LCCDTimeStamp timestamp ) ;

    /** Clear the map and delete all handlers */
    void clear() ;

  private:
    HandlerMap _map ;

    static LCConditionsMgr* _me ;

  };

} //end namespace

#endif // LCConditionsMgr_h
