#ifndef DEPFET_NAMESPACE_H
#define DEPFET_NAMESPACE_H

//! The depfet namespace.
/*! This namespace is used in order not to pollute neither the lcio,
 *  nor the Marlin, not the standard namespaces. It contains all
 *  classes defined by the DEPFET collaboration in order to
 *  develop both the analysis/reconstruction software.
 *
 *  @author Benjamin Schwenker, Göttingen University <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

namespace depfet {}
#endif

#ifndef DEPFET_H
#define DEPFET_H

// system includes <>
#include <iostream>
#include <string>
#include <vector>

namespace depfet
{

  //! Global constants used in the Eutelescope package
  /*!
   * This class has only static data members used only to define global
   * constant to be used within the DEPFET package. Please add here
   * whatever constant you want to use.  A typical useful of this class
   * is to define name of collection to be retrieved/saved from/to
   * files.
   *
   */
  
  class DEPFET
  {
  
  public:
    //! Default destructor.
    /*! This is the default destructor, but it is actually a NO-OP
     *  since there is nothing to be destroyed.
     */
    virtual ~ DEPFET ()  {;  }
    
  public:
    
    
    // Encoding strings
    
    //! Default encoding for full data 
    /*! 
     *  "sensorID:17,uMax:12,vMax:12"
     */
    static const char * MATRIXDEFAULTENCODING;
         
    //! Default encoding for zero suppress data
    /*!
     *  "sensorID:5,sparsePixelType:5"
     */
    static const char * ZSDATADEFAULTENCODING;
    
    //! Cluster default encoding
    /*! 
     *  "sensorID:5,clusterID:8,sparsePixelType:5,quality:5"
     */
    static const char * ZSCLUSTERDEFAULTENCODING;
    
  }; // End class
  
  //! Sparse pixel type enum
  /*! This enumerator is used to define the sparsified pixel type.
   *
   */
  enum SparsePixelType {
    kBaseSparsePixel   = 0,
    kSimpleSparsePixel = 1,
    // add here your implementation
    kUnknownPixelType       = 31
  };

}

#endif
