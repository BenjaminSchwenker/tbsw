// AID-GENERATED
// =========================================================================
// This class was generated by AID - Abstract Interface Definition          
// DO NOT MODIFY, but use the org.freehep.aid.Aid utility to regenerate it. 
// =========================================================================

package hep.lcio.event;


/** Generic class for raw tracker data. It can be used to store the full FADC spectrum as it comes out of the
 *  tracker DAQ or just one ore more single FADC readout values.
 *  @see TrackerData
 *  @see TrackerPulse
 * 
 * @author gaede
 * @version $Id: TrackerRawData.aid,v 1.2 2006/03/24 13:25:53 gaede Exp $
 */

public interface TrackerRawData extends LCObject {

    /** Returns the first detector specific (geometrical) cell id.
     */
    public int getCellID0();

    /** Returns the second detector specific (geometrical) cell id. Optional, check/set 
     *  flag(LCIO::TRAWBIT_ID1)==1.
     */
    public int getCellID1();

    /** Returns a time measurement associated with the adc values, e.g. the 
     *  t0 of the spectrum for the TPC. Subdetector dependent.
     */
    public int getTime();

    /** The actual FADC spectrum.
     */
    public short[] getADCValues();
} // class or interface

