#ifndef _ThreeDModel_h
#define _ThreeDModel_h
	
#include <Eigen/Core>
using Eigen::Matrix3d;
/**  ThreeDModel.h
 *   
 *   Declaration of functions and procedures for 3D manipulation with matrices.
 */

namespace depfet {

// Definition of angles and rotations from Veikko Karimaki
void FillRotMatrixKarimaki(Matrix3d & m, double alpha, double beta, double gamma);
void GetAnglesKarimaki(const Matrix3d & m, double & alpha, double & beta, double & gamma);
void testAnglesKarimaki();

} // Namespace

#endif
