#ifndef _AlignableDet_h
#define _AlignableDet_h



#include <Eigen/Core>

typedef Eigen::Matrix<double,6,1> SensorAlignmentParameters;
typedef Eigen::Matrix<double,6,6> SensorAlignmentCovariance;
typedef Eigen::Matrix<double,2,6> SensorAlignmentJacobian; 	
#include<Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(SensorAlignmentParameters)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(SensorAlignmentCovariance)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(SensorAlignmentJacobian)

#include <vector>
namespace depfet {
 	
//! Class AlignableDet
/* 
 * Class AlignableDet manages the results of track based 
 * alignment algorithms. 
 *  
 * @Author B. Schwenker, University of Göttingen
 * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class AlignableDet
{
 public:
  
  /** Constructor  
   */
  AlignableDet(int nAlignables);
  
  /** Getter
   */
  SensorAlignmentParameters GetAlignState(int i);
  SensorAlignmentCovariance GetAlignCovariance(int i); 
  
  /** Setter
   */
  void SetAlignState(int i, const SensorAlignmentParameters& a);
  void SetAlignCovariance(int i, const SensorAlignmentCovariance& E);
   
  /* Alignment parameters (all alignables) 
   */
  std::vector< SensorAlignmentParameters > alignmentParametersVec;    

  /* Alignment covariance (all alignables)   
   */  
  std::vector<SensorAlignmentCovariance> alignmentCovarianceVec;      
};

} // Namespace
#endif

