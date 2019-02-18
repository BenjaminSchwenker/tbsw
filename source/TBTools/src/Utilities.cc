/**  Utilities.cc Contains often used functions for vectors,  matrices and other. */

#include "Utilities.h"

 	
// ROOT includes
#include <TMath.h>
#include <TFile.h>
#include <TROOT.h>
 	
using namespace std;


namespace depfet {

 	
std::vector<std::string> splitpath( const std::string& str, const std::set<char> delimiters)
{
  std::vector<std::string> result;
  
  char const* pch = str.c_str();
  char const* start = pch;
  for(; *pch; ++pch)
  {
    if (delimiters.find(*pch) != delimiters.end())
    {
      if (start != pch)
      {
        std::string str(start, pch);
        result.push_back(str);
      }
      else
      {
        result.push_back("");
      }
      start = pch + 1;
    }
  }
  result.push_back(start);
  
  return result;
}
 	
/** Test if two values are equal within a certain precision. */
bool equal( double val1, double val2, double precision)
{
  if (val1 == 0)
    return val2 < precision;
  if (val2 == 0)
    return val1 < precision;
  
  //   cout << "val1= " << val1 << "   val2= " << val2 << "  equality= " << TMath::Abs(2*(val1-val2)/(val1+val2)) << endl;
  return TMath::Abs(2*(val1-val2)/(val1+val2)) < precision;
}

	
// Utility routines
//==================
 	
TObject * get_object( const char * objectname, const char * filename)
{
  // try to find out if file is already opened
  TFile * f = (TFile *) gROOT->GetListOfFiles()->FindObject(filename);
  if (f != 0) {
    //     cout << "File " << fname << " was already open" << endl;
  }
  else 
  {
    f = new TFile(filename, "READ");
    if (f == 0)
    {
      cout << "Error creating TFile object" << endl;
      return 0;
    }
    if (!f->IsOpen()) 
    {
      cout << "Could not open file " << filename << endl;
      return 0;
    }
  }
  TObject * obj = (TTree *) f->Get(objectname);
  if (obj == 0) {
    cout << "Error reading object " << objectname << " from file " << filename << endl;
  }
  return obj;
}

TTree * get_tree( const char * treename, const char * filename)
{
  return (TTree *) get_object(treename, filename);
}
 	
TTree * init_align_tree( const char * filename)
{
  TTree * t = get_tree("AlignTree", filename);
  return t;
}

} // Namespace;
 	
