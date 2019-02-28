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

// to_string(intType) used do convert any integer type  (up to unsigned long long) to string
// waaaaay faster than normal stringstream method. From stackoverflow.

namespace _utils_toStr_impl {
    template <typename T>
    T reduce2(T v) {
        T k = ((v * 410) >> 12) & 0x000F000F000F000Full;
        return (((v - k * 10) << 8) + k);
    }
    template <typename T>
    T reduce4(T v) {
        T k = ((v * 10486) >> 20) & 0xFF000000FFull;
        return reduce2(((v - k * 100) << 16) + (k));
    }
    typedef unsigned long long ull;
    inline ull reduce8(ull v) {
        ull k = ((v * 3518437209u) >> 45);
        return reduce4(((v - k * 10000) << 32) + (k));
    }
}

template <typename T>
std::string to_string(T o) {
    union {
        char str[16];
        unsigned short u2[8];
        unsigned u4[4];
        unsigned long long u8[2];
    };

    unsigned v = o < 0 ? ~o + 1 : o;

    u8[0] = (_utils_toStr_impl::ull(v) * 3518437209u) >> 45;
    u8[0] = (u8[0] * 28147497672ull);
    u8[1] = v - u2[3] * 100000000;

    u8[1] = _utils_toStr_impl::reduce8(u8[1]);
    char* f;
    if (u2[3]) {
        u2[3] = _utils_toStr_impl::reduce2(u2[3]);
        f = str + 6;
    } else {
        unsigned short* k = u4[2] ? u2 + 4 : u2 + 6;
        f = *k ? (char*)k : (char*)(k + 1);
    }
    if (!*f) f++;

    u4[1] |= 0x30303030;
    u4[2] |= 0x30303030;
    u4[3] |= 0x30303030;
    if (o < 0) *--f = '-';
    return std::string(f, (str + 16) - f);
}


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
