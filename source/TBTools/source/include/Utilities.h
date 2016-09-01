#ifndef _Utilities_h
#define _Utilities_h
 	
// C++ / STL includes
#include <string>
#include <iostream>
#include <sstream>
#include <assert.h>
 	
// ROOT includes
#include "TString.h" // for char *Form(...)
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TObject.h"
#include "TTree.h"
 	
// CLHEP includes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>

namespace depfet {

/** Utilities.h
 *
 * Various preprocessor macros, variables, functions, and procedures
 * that misfit elsewhere.
 **/

// Used to generate key strings for decoding/encoding
std::string to_string( int key );  	
 	 
bool equal(double val1, double val2, double precision);
 	
//////////////////////////////////////////////////////////////////////
// Conversion of matrices
 	
void CLHEPtoROOT( CLHEP::HepMatrix & oldM, TMatrixD * newM);
void CLHEPtoROOT( CLHEP::HepSymMatrix & oldM, TMatrixDSym * newM);
void CLHEPtoROOT( CLHEP::HepVector & oldV, TVectorD * newV);
 	
//////////////////////////////////////////////////////////////////////
// Find objects from ROOT files
 	
TObject * get_object( const char * objectname, const char * filename);
TTree * get_tree( const char * treename, const char * filename);
TTree * init_align_tree( const char * filename);

} // Namespace 
	
#endif
