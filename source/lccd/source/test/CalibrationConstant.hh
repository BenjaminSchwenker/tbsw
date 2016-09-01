#ifndef CalibrationConstant_h
#define CalibrationConstant_h 1

#include "lcio.h"
#include "UTIL/LCFixedObject.h"

#define NINT 1 
#define NFLOAT 2
#define NDOUBLE 0

#define ID_INDEX 0
#define OFFSET_INDEX 0
#define GAIN_INDEX 1

using namespace lcio ;

class CalibrationConstant ;

/** Example for a simple calibration class based on the LCFixedObject template.
 *  <p>
 *  LCFixedObject uses an instance of LCGenericObjectImpl that holds the data, thus 
 *  there is no overhead when the data is read from a database or file
 *  for copying it  to some local structure (Decorator pattern).<br>
 *  
 */
class CalibrationConstant : public LCFixedObject<NINT,NFLOAT,NDOUBLE> {
  
public: 
  
  /** Convenient c'tor.
   */
  CalibrationConstant(int cellID, float offset, float gain) {

    obj()->setIntVal( ID_INDEX , cellID  ) ;
    obj()->setFloatVal( OFFSET_INDEX ,  offset ) ;  
    obj()->setFloatVal( GAIN_INDEX ,  gain ) ;  
  }

  /** 'Copy constructor' needed to interpret LCCollection read from file/database.
   */
  CalibrationConstant(LCObject* obj) : LCFixedObject<NINT,NFLOAT,NDOUBLE>(obj) { } 

  /** Important for memory handling*/
  virtual ~CalibrationConstant() { /* no op*/  }
  

  // the class interface:
  int getCellID()   { return obj()->getIntVal( ID_INDEX ) ;  } 
  float getOffset() { return obj()->getFloatVal( OFFSET_INDEX )  ; } 
  float getGain()   { return obj()->getFloatVal( GAIN_INDEX )  ; } 
  
  void print(  std::ostream& os ) ;
  

  // -------- need to implement abstract methods from LCGenericObject

  const std::string getTypeName() const { 
    return"CalibrationConstant" ;
  } 
  
  const std::string getDataDescription() const {
    return "i:cellID,f:offset,f:gain" ;
  }

}; // class


#endif 
//=============================================================================
