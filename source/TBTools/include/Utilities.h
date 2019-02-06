#ifndef _Utilities_h
#define _Utilities_h
 	
// C++ / STL includes
#include <string>
#include <iostream>
#include <sstream>
#include <assert.h>
#include <vector> 	

// ROOT includes
#include "TString.h" // for char *Form(...)
#include "TObject.h"
#include "TTree.h"
 	
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
// Find objects from ROOT files
 	
TObject * get_object( const char * objectname, const char * filename);
TTree * get_tree( const char * treename, const char * filename);
TTree * init_align_tree( const char * filename);

//////////////////////////////////////////////////////////////////////
// Split a string into using char delimiter

std::vector<std::string> splitpath( const std::string& str, const std::set<char> delimiters);


} // Namespace 
	
#endif
