// AID-GENERATED
// =========================================================================
// This class was generated by AID - Abstract Interface Definition          
// DO NOT MODIFY, but use the org.freehep.aid.Aid utility to regenerate it. 
// =========================================================================

package hep.lcio.event;


/** A generic simulated tracker hit. 
 * 
 * @author gaede
 * @version $Id: SimTrackerHit.aid,v 1.9 2006/03/24 13:25:52 gaede Exp $
 */

public interface SimTrackerHit extends LCObject {

    /**Returns the detector specific (geometrical) cell id.
     */
    public int getCellID();

    /** Returns the hit  position in [mm].	
     */
    public double[] getPosition();

    /** Returns  the dE/dx of the hit in [GeV].
     */ 	
    public float getdEdx();

    /** Returns the  time of the hit in [ns]. TO DO needs definition.
     */
    public float getTime();

    /** Returns the MC particle that caused the hit.
     *
     * @see MCParticle
     */
    public MCParticle getMCParticle();

    /** Returns the 3-momentum of the particle at the hits position in [GeV] - 
     * optional, only if bit LCIO::THBIT_MOMENTUM is set.	
     */ 
    public float[] getMomentum();

    /** The path length of the particle in the sensitive material that resulted in this hit.
     *  This is only stored together with momentum, i.e. if  LCIO::THBIT_MOMENTUM is set.
     */
    public float getPathLength();
} // class or interface

