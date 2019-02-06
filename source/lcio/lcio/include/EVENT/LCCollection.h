// -*- C++ -*-
// AID-GENERATED
// =========================================================================
// This class was generated by AID - Abstract Interface Definition          
// DO NOT MODIFY, but use the org.freehep.aid.Aid utility to regenerate it. 
// =========================================================================
#ifndef EVENT_LCCOLLECTION_H
#define EVENT_LCCOLLECTION_H 1

#include <string>

#include "EVENT/LCParameters.h"
#include "Exceptions.h"

namespace EVENT {

class LCObject;

/** The generic collection used in LCIO. Holds arbitrary objects of type LCObject.
 * 
 * @author gaede 
 * @version $Id: LCCollection.aid,v 1.11 2005/03/02 16:22:58 gaede Exp $
 * @see LCObject
 * @see LCIO
 */
           
                                                    
     
class LCCollection {

public: 
    /// Destructor.
    virtual ~LCCollection() { /* nop */; }

    /**Returns the number of elements in the collection.
     * Same as size().
     */
    virtual int getNumberOfElements() const = 0;

    /** Returns the type name of the collection - valid names are defined in LCIO.
     */
    virtual const std::string & getTypeName() const = 0;

    /** Returns pointer to element at index - no range check, use getNumberOfEntries().
     */
    virtual LCObject * getElementAt(int index) const = 0;

    /** Returns flag word for collection. Bits 16-31 are reserved for LCIO
     *  Depending on the object type stored they have a special meaning, e.g. 
     *  for SimCalorimeterHits: <br>
     *  CHBIT_LONG = 31   -  store position <br>
     *  CHBIT_BARREL = 30 -  endcap or barrel <br>
     *  CHBIT_ID1 = 29 -   cellid1 is sored <br>
     *  CHBIT_PDG = 28 - store pdg of secondaries <br>
     *  &nbsp;<br>
     *  Bit 16 is used to flag collection as transient <br>
     *  Bit 17 is used to flag collection as default <br>
     *  Bit 18 is used to flag collection as subset <br>
     *  Bits 0-15 are subdetector/user specific.
     * @see isTransient()
     */
    virtual int getFlag() const = 0;

    /**Transient bit in flag word.
     */
          

    static const int BITTransient = 16 ;
    static const int BITDefault   = 17 ;
    static const int BITSubset    = 18 ;
    /** True if collection is transient, i.e. will not be written to any LCIO file.
     *  Convenient method that checks bit BITTransient of the flag word.
     */
    virtual bool isTransient() const = 0;

    /** True if collection is the default collection for the given type.
     *  This implies that the collection is complete and unambigous.
     *  Convenient method that checks bit BITDefault of the flag word.
     */
    virtual bool isDefault() const = 0;

    /** True if the collection holds a subset of objects from other collections. 
     *  If the collection is not transient only the pointers/references to the original
     *  objects will be stored.
     *  Convenient method that checks bit BITSubset of the flag word.
     */
    virtual bool isSubset() const = 0;

    /** Adds the given element to (end of) the collection. Throws an exception 
     * if the collection (event) is 'read only'.
     *
     * @throws ReadOnlyException
     */
    virtual void addElement(LCObject * obj) = 0;

    /** Removes the i-th element from the collection. Throws an exception 
     * if the collection (event) is 'read only'.
     *
     * @throws ReadOnlyException
     */
    virtual void removeElementAt(int i) = 0;

    /** Set the flag word. This is allowed in 'read only' mode.
     */
    virtual void setFlag(int flag) = 0;

    /** Parameters defined for this collection.
     */
    virtual const LCParameters & getParameters() const = 0;

   /** Parameters defined for this collection.
    * Can be used to modify the collection paramters. 
    */
    virtual LCParameters & parameters() = 0;
}; // class
} // namespace EVENT
#endif /* ifndef EVENT_LCCOLLECTION_H */
