// AID-GENERATED
// =========================================================================
// This class was generated by AID - Abstract Interface Definition          
// DO NOT MODIFY, but use the org.freehep.aid.Aid utility to regenerate it. 
// =========================================================================

package hep.lcio.event;

import java.util.List;

/** The LCIO cluster.
 * 
 * @author gaede
 * @version $Id: Cluster.aid,v 1.14 2006/08/03 16:53:34 gaede Exp $
 */


public interface Cluster extends LCObject {

    /** Flagword that defines the type of cluster. Bits 0-15 can be used to denote the subdetectors
     *  that have contributed hits to the cluster. For the definition of the bits 
     *  check/Set the collection variables ClusterTypeBitNames and ClusterTypeBitIndices.
     *  </br>Bits 16-31 are used internally.
     */
    public int getType();

    /** Energy of the cluster.
     */
    public float getEnergy();

    /** Position of the cluster.
     */
    public float[] getPosition();

    /** Covariance matrix of the position (6 Parameters)
    */
    public float[] getPositionError();

    /** Intrinsic direction of cluster at position: Theta.
     * Not to be confused with direction cluster is seen from IP.
     */
    public float getITheta();

    /** Intrinsic direction of cluster at position: Phi.
     * Not to be confused with direction cluster is seen from IP.
     */
    public float getIPhi();

    /** Covariance matrix of the direction (3 Parameters)
     */
    public float[] getDirectionError();

    /** Shape parameters - check/set  collection parameter
     *  ClusterShapeParameters for size and names of parameters.
     */
    public float[] getShape();

//     /** Type hypotheses: 3 Parameters: compatible with EM, HAD, muon cluster
//      */
//     public const FloatVec& getParticleType() const ;
    /** The particle Id's sorted by their likelihood.
     * @see ParticleID
     */
    public List getParticleIDs();

    /** The clusters that have been combined to this cluster.
     */
    public List getClusters();

    /** The hits that have been combined to this cluster.
     *  Only available if collection flag bit LCIO::CLBIT_HITS==1 and if 
     *  the CalorimeterHit objects have not been saved with LCIO::RCHBIT_NO_PTR==1.
     *  @see CalorimeterHit
     */
    public List getCalorimeterHits();

    /** Returns the energy contribution of the hits 
     *  Runs parallel to the CalorimeterHitVec from getCalorimeterHits()
     */
    public float[] getHitContributions();

    /** A vector that holds the energy observed in a particular subdetectors.
     *  The mapping of indices to subdetectors is implementation dependent.
     *  To be used as convenient information or if hits are not stored in 
     *  the data set, e.g. DST or FastMC. 
     *  Check/set collection parameter ClusterSubdetectorNames for decoding the
     *  indices of the array.
     */
    public float[] getSubdetectorEnergies();
} // class or interface

