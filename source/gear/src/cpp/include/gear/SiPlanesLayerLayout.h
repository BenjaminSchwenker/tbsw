#ifndef GEAR_SIPLANESLAYERLAYOUT_H
#define GEAR_SIPLANESLAYERLAYOUT_H 1

#include <vector>
#include <tuple>

namespace gear {

/**Abstract description of layers in a pixel beam telescope.
 * @author T. Klimkovich, DESY
 * @author B. Schwenker, Uni Göttingen
 */
class SiPlanesLayerLayout {

public: 
    /// Destructor.
    virtual ~SiPlanesLayerLayout() { /* nop */; }

    /** The total number of layers. */
    
    virtual int getNLayers() const = 0;
    
    /** Thickness of nonsensitive volume of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */	
    virtual double getLayerThickness(int layerIndex) const = 0;
    
    /** The radiation length of nonsensitive volume of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source.
     */
    virtual double getLayerRadLength(int layerIndex) const = 0;
    
    /** Size in x direction of nonsensitive volume of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */	
    virtual double getLayerSizeX(int layerIndex) const = 0;
     
    /** Size in y direction of nonsensitive volume of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */	
    virtual double getLayerSizeY(int layerIndex) const = 0;
    
    /** Atomic number of nonsensitive volume of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */	
    virtual double getLayerAtomicNumber(int layerIndex) const = 0;
    
    /** Atomic mass of nonsensitive volume of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */	
    virtual double getLayerAtomicMass(int layerIndex) const = 0;
    
    /** ID of sensitive volume of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source.*/
    
    virtual int getSensitiveID(int layerIndex) const = 0;
     
    /** x position of the center of sensitive volume of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source.
     */
    virtual double getSensitivePositionX(int layerIndex) const = 0;

    /** y position of the center of sensitive volume of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source.
     */
    virtual double getSensitivePositionY(int layerIndex) const = 0;

    /** z position of the center of sensitive volume of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source.
     */
    virtual double getSensitivePositionZ(int layerIndex) const = 0;
    
    /** Thickness of sensitive volume of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */	
    virtual double getSensitiveThickness(int layerIndex) const = 0;

    /** The radiation length of sensitive volume of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source.
     */
    virtual double getSensitiveRadLength(int layerIndex) const = 0;
    
    /** Atomic number of sensitive volume of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */	
    virtual double getSensitiveAtomicNumber(int layerIndex) const = 0;
    
    /** Atomic mass of sensitive volume of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */	
    virtual double getSensitiveAtomicMass(int layerIndex) const = 0;
     
    /** 3D Euler rotation angle alpha  
     *  of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */
    virtual double getSensitiveRotationAlpha(int layerIndex) const = 0;

    /** 3D Euler rotation angle beta  
     *  of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */
    virtual double getSensitiveRotationBeta(int layerIndex) const = 0;

    /** 3D Euler rotation angle gamma  
     *  of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */
    virtual double getSensitiveRotationGamma(int layerIndex) const = 0;

    /** First element (cos(theta)) of rotation matrix of sensitive volume 
     *  of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */
    virtual double getSensitiveRotation1(int layerIndex) const = 0;

    /** Second element (-sin(theta)) of rotation matrix of sensitive volume 
     *  of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */
    virtual double getSensitiveRotation2(int layerIndex) const = 0;

    /** Third element (sin(theta)) of rotation matrix of sensitive volume 
     *  of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */
    virtual double getSensitiveRotation3(int layerIndex) const = 0;

    /** Fourth element (cos(theta)) of rotation matrix of sensitive volume 
     *  of layer layerIndex - layer indexing starts at 0
     *  for the layer closest to the beam source. 
     */
    virtual double getSensitiveRotation4(int layerIndex) const = 0;

    /** Vector with definition of readout cells along sensor u axis 
     */
    virtual std::vector< std::tuple<int,int,double> > getSensitiveUCells(int layerIndex) const = 0; 

    /** Vector with definition of readout cells along sensor u axis 
     */
    virtual std::vector< std::tuple<int,int,double> > getSensitiveVCells(int layerIndex) const = 0; 

}; // class
} // namespace gear
#endif /* ifndef GEAR_SIPLANESLAYERLAYOUT_H */
