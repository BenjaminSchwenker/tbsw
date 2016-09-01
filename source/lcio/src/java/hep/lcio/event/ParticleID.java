// AID-GENERATED
// =========================================================================
// This class was generated by AID - Abstract Interface Definition          
// DO NOT MODIFY, but use the org.freehep.aid.Aid utility to regenerate it. 
// =========================================================================

package hep.lcio.event;


/** Persistent interface for LCIO ParticleIDs.
 *  Used by ReconstructedParticle and Cluster  
 *  for different hypotheses on the particle type.
 * 
 * @author gaede
 * @version $Id: ParticleID.aid,v 1.11 2008/05/30 13:23:27 gaede Exp $
 * @see ReconstructedParticle.getParticleIDs()
 * @see Cluster.getParticleIDs()
 */
public interface ParticleID extends LCObject {

    /** Type - userdefined.
     */
    public int getType();

    /** The PDG code of this id - UnknownPDG ( 999999 ) if unknown.
     */
    public int getPDG();

    /**The likelihood  of this hypothesis - in a user defined normalization.
     */
    public float getLikelihood();

    /** Type of the algorithm/module that created this hypothesis - NOTE: must be unique within one 
     *  collection.  
     *  Check/set collection parameters PIDAlgorithmTypeName and PIDAlgorithmTypeID.
     */
    public int getAlgorithmType();

    /** Parameters associated with this hypothesis.
     * Check/set collection parameters ParameterNames_PIDAlgorithmTypeName for decoding the indices.
     */
    public float[] getParameters();

    /** Constant to be used if the PDG code is not known or undefined.
     */
          
     
                                          
      
    public final static int UnknownPDG = 999999;
} // class or interface

