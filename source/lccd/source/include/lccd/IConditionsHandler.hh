#ifndef IConditionsHandler_h
#define IConditionsHandler_h

// 
#include "lccd.h"
#include "IConditionsChangeListener.hh"


// -- lcio headers
#include "lcio.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"


namespace lccd {
  
  /** Abstract handler for conditions data. Need one implementation instance 
   *  for every data set. Check documentation of concrete implemetation for
   *  behaviour on update, e.g. a simple file handler might ignore calls to 
   *  update and allways hold the same data set.
   * 
   *  @author F.Gaede, DESY
   *  @version $Id: IConditionsHandler.hh,v 1.3 2007-04-27 13:21:31 gaede Exp $
   */

  class IConditionsHandler {
    
  public:
    
    virtual ~IConditionsHandler() {} ;


    /** Retrieves the new conditions data if required by evt->getTimestamp() 
     * and adds a collection to event with its name.
     * Notifies all registered listeners
     */
    virtual void updateEvent( lcio::LCEvent* evt ) = 0;
    
    /** Retrieves the new conditions data if required by timestamp.
     *  Notifies all registered listeners<br>
     *  See implementation classes for details.
     */
    virtual void update( LCCDTimeStamp timestamp ) = 0;
    
    /** Returns the current collection of conditions data
     */
    virtual lcio::LCCollection* currentCollection() = 0;
    
    
    /** Every IConditionsChangeListener will be notified if the conditions data
     *  of this instance has changed.
     */
    virtual void registerChangeListener( IConditionsChangeListener* cl) = 0;
    
    
    /** Remove the specified listener from list of registered listeners ;
     */
    virtual void removeChangeListener( IConditionsChangeListener* cl) = 0;
    
    
    /** The name of this conditions handler - used, e.g. as name of the collection
     *  in update( lcio::LCEvent* evt ).
     *
     *  @see update( lcio::LCEvent* evt )
     */
    virtual const std::string name() = 0; 
    
    /** First timestamp for which the current conditions data is valid.
     */
    virtual LCCDTimeStamp validSince() =0 ;
    
    /** Last timestamp for which the current conditions data is valid.
     */
    virtual LCCDTimeStamp validTill() =0;
    
  };

} //end namespace

#endif // IConditionsHandler_h
