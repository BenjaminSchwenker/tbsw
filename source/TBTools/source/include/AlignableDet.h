#ifndef _AlignableDet_h
#define _AlignableDet_h

#include <vector>

#include <Eigen/Core>

typedef Eigen::Matrix<double,6,1> SensorAlignmentParameters;
typedef Eigen::Matrix<double,6,6> SensorAlignmentCovariance;
typedef Eigen::Matrix<double,2,6> SensorAlignmentJacobian; 	

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
  void SetAlignState(int i, SensorAlignmentParameters a);
  void SetAlignCovariance(int i, SensorAlignmentCovariance E);  
   
  /* Alignment parameters (all alignables) 
   */
  std::vector< SensorAlignmentParameters > alignmentParametersVec;    

  /* Alignment covariance (all alignables)   
   */  
  std::vector<SensorAlignmentCovariance> alignmentCovarianceVec;      
};

} // Namespace
#endif

