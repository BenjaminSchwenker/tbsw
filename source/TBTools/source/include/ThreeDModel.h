#ifndef _ThreeDModel_h
#define _ThreeDModel_h
	
// CLHEP includes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
 	
/**  ThreeDModel.h
 *   
 *   Declaration of functions and procedures for 3D manipulation with matrices.
 */
	
namespace depfet {

// Definition of angles and rotations from Veikko Karimaki
void FillRotMatrixKarimaki(CLHEP::HepMatrix & m, double alpha, double beta, double gamma);
void GetAnglesKarimaki(const CLHEP::HepMatrix & m, double & alpha, double & beta, double & gamma);
void testAnglesKarimaki();

} // Namespace

#endif
