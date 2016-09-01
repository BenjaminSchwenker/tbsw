#ifndef GEAR_GearParametersImpl_H
#define GEAR_GearParametersImpl_H 1

#include <map>

#include "gear/GearParameters.h"
#include "gear/GEAR.h"

namespace gear {
    
  /** Implementation of  GearParameters - a set off parameters that can be used to describe the 
   *  geometrical properties of a specific subdetector type for reconstruction.
   *  It basically allows to retrieve named parameters of type int, string and double
   *  and arrays/vectors of these types and set these. 
   *
   * @author F. Gaede, DESY
   * @version $Id: GearParametersImpl.h,v 1.1.1.1 2005-09-13 14:41:03 gaede Exp $
   */
  class GearParametersImpl : virtual public GearParameters {
	
  public: 

    typedef std::map< std::string, int > IntMap ;
    typedef std::map< std::string, double > DoubleMap ;
    typedef std::map< std::string, std::string > StringMap ;
    typedef std::map< std::string, IntVec > IntVecMap ;
    typedef std::map< std::string, DoubleVec > DoubleVecMap ;
    typedef std::map< std::string, StringVec > StringVecMap ;
	
	
    /// Destructor.
    virtual ~GearParametersImpl() ;
	
    /** Integer value for key.     
     *  @throws UnknownParameterException
     */
    virtual int getIntVal(const std::string & key) const throw (UnknownParameterException, std::exception ) ;
	
    /** Double value for key      
     *  @throws UnknownParameterException
     */
    virtual double getDoubleVal(const std::string & key) const throw (UnknownParameterException, std::exception ) ;
	
    /** String value for key      
     *  @throws UnknownParameterException
     */
    virtual const std::string & getStringVal(const std::string & key) const 
      throw (UnknownParameterException, std::exception ) ;
	
    /** Integer values for key      
     *  @throws UnknownParameterException
     */
    virtual const std::vector<int>  & getIntVals(const std::string & key) const 
      throw (UnknownParameterException, std::exception ) ;
	
    /** Double values for key      
     *  @throws UnknownParameterException
     */
    virtual const std::vector<double>  & getDoubleVals(const std::string & key) const 
      throw (UnknownParameterException, std::exception ) ;
	
    /** String values for key      
     *  @throws UnknownParameterException
     */
    virtual const std::vector<std::string>  & getStringVals(const std::string & key) const 
      throw (UnknownParameterException, std::exception ) ;
	
    // setters

    /** Set Integer value for key.     
     */
    virtual void  setIntVal(const std::string & key , int val ) ;
	
    /** Double value for key      
     */
    virtual void  setDoubleVal(const std::string & key, double val ) ;
	
    /** String value for key      
     */
    virtual void  setStringVal(const std::string & key , const std::string & val) ;
	
    /** Integer values for key      
     */
    virtual void  setIntVals(const std::string & key, const std::vector<int>& vals) ;
	
    /** Double values for key      
     */
    virtual void  setDoubleVals(const std::string & key, const std::vector<double>& vals) ;
	
    /** String values for key      
     */
    virtual void  setStringVals(const std::string & key, const std::vector<std::string>& vals) ;



    /** All keys of int variables. */ 
    virtual const std::vector<std::string>  & getIntKeys() const ;

    /** All keys of double variables. */ 
    virtual const std::vector<std::string>  & getDoubleKeys() const ;

    /** All keys of string variables. */ 
    virtual const std::vector<std::string>  & getStringKeys() const ;

    /** All keys of IntVec variables. */ 
    virtual const std::vector<std::string>  & getIntVecKeys() const ;

    /** All keys of DoubleVec variables. */ 
    virtual const std::vector<std::string>  & getDoubleVecKeys() const ;

    /** All keys of StringVec variables. */ 
    virtual const std::vector<std::string>  & getStringVecKeys() const ;

  protected:
    
    // maps holding named parameters
    std::map< std::string, int > _intMap ;
    std::map< std::string, double > _doubleMap ;
    std::map< std::string, std::string > _stringMap ;
    std::map< std::string, IntVec > _intVecMap ;
    std::map< std::string, DoubleVec > _doubleVecMap ;
    std::map< std::string, StringVec > _stringVecMap ;

    mutable StringVec _intKeys ;
    mutable StringVec _doubleKeys ;
    mutable StringVec _stringKeys ;
    mutable StringVec _intVecKeys ;
    mutable StringVec _doubleVecKeys ;
    mutable StringVec _stringVecKeys ;
              

  }; // class
} // namespace gear

#endif /* ifndef GEAR_GEARPARAMETERS_H */
