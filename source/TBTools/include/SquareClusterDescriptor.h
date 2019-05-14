#ifndef SquareClusterDescriptor_H
#define SquareClusterDescriptor_H

#include "Det.h"
#include "PixelCluster.h"

#include <string>


namespace depfet {
  
 
  //! SquareClusterDescriptor
  /*! This class represents a sorted square pixel cluster 
   * 
   *  @Author Benjamin Schwenker, Universitaet Goettingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */ 
  
  class SquareClusterDescriptor  {
  
  public:

    /** Default constructor */
    SquareClusterDescriptor():
      m_originU(0), m_originV(0)
    {}
    
    /** Constructor */
    SquareClusterDescriptor(const PixelCluster& Cluster, const Det& Sensor) ;

    /** Get origin coordinate
     */
    float getOriginU() const { return  m_originU; }

    /** Get origin coordinate
     */
    float getOriginV() const { return  m_originV; }
    
    /** Get cluster shape string 
     * @return shape string containing all features of cluster used for position reconstruction 
     */
    std::string getShape(int vCellPeriod = 1, int uCellPeriod = 1, int etaBin=0) const; 
     
    /** Get cluster type string. 
     * @return type string, same as shape but without eta information
     */
    std::string getType(int vCellPeriod = 1, int uCellPeriod = 1) const;

    /** Get cluster shape string
     * @return shape string containing all features of cluster used for position reconstruction
     */
    std::string getEtaBinString(int etaBin=0) const;
    
    /** Get cluster eta
    */
    double computeEta(double thetaU, double thetaV) const; 
    
    /** Get cluster eta bin  
    */
    static int computeEtaBin(double eta, const std::vector<double>& etaBinEdges); 
    
    /** Get index of head pixel in cluster
    */
    int getHeadPixelIndex(double thetaU, double thetaV) const; 
    
    /** Get index of tail pixel in cluster
    */
    int getTailPixelIndex(double thetaU, double thetaV) const;
    
    /** Get indes of last pixel of local row vOffset  
    */
    int getLastPixelWithVOffset(int vOffset) const; 
     
    /** Get index of first pixel of local row vOffset 
    */
    int getFirstPixelWithVOffset(int vOffset) const; 
    
  protected:
    
    // Origin of cluster coordinates
    float m_originU;   
    // Origin of cluster coordinates
    float m_originV; 
    // Cluster size in vcells     
    unsigned short m_vSize;  
    // Start ucell of the cluster      
    unsigned short m_uStart; 
    // Start vcell of the cluster      
    unsigned short m_vStart;   
    // Pixeltype 
    int m_pixelType;
   
     
    // Sorted vector of raw digits 
    std::vector<RawDigit> m_sortedDigits;
  };
}
#endif
