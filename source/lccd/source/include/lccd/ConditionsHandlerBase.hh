#ifndef ConditionsHandlerBase_h
#define ConditionsHandlerBase_h

#include "lccd.h"

#include "IConditionsHandler.hh"

#include <vector>


//using namespace lcio ;

namespace lccd {
  
  /** Base implementation of IConditionsHandler. Provides common code, e.g.
   *  registration and notfication of IConditionsChangeListeners.<br>
   *  Subclasses need to implement IConditionsHandler::update( LCCDTimeStamp timestamp ).
   * 
   *  @author F.Gaede, DESY
   *  @version $Id: ConditionsHandlerBase.hh,v 1.3 2007-04-27 13:21:31 gaede Exp $
   */
  
  class ConditionsHandlerBase : public IConditionsHandler {
    
  public:
    
    /** Default c'tor - need a name for the conditions data set.
     */
    ConditionsHandlerBase(const std::string& name) : 
      _name( name),
      _validSince( LCCDPlusInf ) ,   // initially we don't have a valid time range
      _validTill(  LCCDMinusInf ) {
    } 
    
    virtual ~ConditionsHandlerBase() {}
    

    /** Retrieves the new conditions data if required by evt->getTimestamp() 
     * and adds a collection to event with its name.
     * Notifies all registered listeners
     */
    virtual void updateEvent( lcio::LCEvent* evt ) ;
   

    /** Retrieves the new conditions data if required by timestamp.
     *  Notifies all registered listeners<br>
     *  See implementation classes for details.
     */
    virtual void update( LCCDTimeStamp timestamp ) =0 ;
    

    /** Returns the current collection of conditions data
     */
    virtual lcio::LCCollection* currentCollection() { return _col ; } 
     
    
    /** Every IConditionsChangeListener will be notified if the conditions data
     *  of this instance has changed.
     */
    virtual void registerChangeListener( IConditionsChangeListener* cl ) ;
    
    /** Remove the specified listener from list of registered listeners ;
     */
    virtual void removeChangeListener( IConditionsChangeListener* cl) ;

    /** The name of this conditions handler - used, e.g. as name of the collection
     *  in update( lcio::LCEvent* evt ).
     *
     *  @see update( lcio::LCEvent* evt )
     */
    virtual const std::string name() { return _name ; }  

    /** First timestamp for which the current conditions data is valid.
     */
    virtual LCCDTimeStamp validSince() { return _validSince ; }  
    
    /** Last timestamp for which the current conditions data is valid.
     */
    virtual LCCDTimeStamp validTill() { return _validTill ; }
    

  protected :
    
    ConditionsHandlerBase() {} 

    /** Tell all registered listeners that the conditions data has changed.
     */
    virtual void notifyListeners() ;
    

    // ----  data members ------

    /** The current collection of conditions data.*/
    lcio::LCCollection* _col ;

    /** Name assigned with conditions data*/
    std::string _name ;

    /** First timestamp for which the current data is valid.*/
    LCCDTimeStamp _validSince  ;

    /** Last timestamp for which the current data is valid.*/
    LCCDTimeStamp _validTill ;

    /** The registered change listeners */
    std::vector< IConditionsChangeListener* > _changeListeners ;
    

  };
}  //end namespace

#endif // ConditionsHandlerBase_h
