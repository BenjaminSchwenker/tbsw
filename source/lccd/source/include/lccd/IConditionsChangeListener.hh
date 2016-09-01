#ifndef IConditionsChangeListener_h
#define IConditionsChangeListener_h

#include "lcio.h"
#include "EVENT/LCCollection.h"

namespace lccd {
  
/** Simple interface that allows notification of implementation classes if
 *  a conditions data set has changed.
 * 
 *  @author F.Gaede, DESY
 *  @version $Id: IConditionsChangeListener.hh,v 1.2 2007-04-27 13:21:31 gaede Exp $
 */

class IConditionsChangeListener {
  
public:

  virtual ~IConditionsChangeListener() {} 

  /** The callback - use updated conditions data collection as needed, e.g.
   *  repopulate a calibration map. Implementors can access the collection name 
   *  through the collection parameter "CollectionName": <br>
   *  <pre>  
   *    std::string colName = col->getParameters().getStringVal("CollectionName") ;
   *  </pre>
   *  if needed. 
   */ 
  virtual void conditionsChanged( lcio::LCCollection* col ) = 0 ;
  
};

}  //end namespace

#endif // IConditionsChangeListener_h
