/**  Utilities.cc Contains often used functions for vectors,  matrices and other. */

#include "Utilities.h"

// C++ includes
#include <assert.h>
#include <iostream>
 	
// ROOT includes
#include <TMath.h>
#include <TFile.h>
#include <TROOT.h>
 	
using namespace std;
using namespace CLHEP;

namespace depfet {

/** Used to generate string keys.*/
std::string to_string( int key ) {

  std::ostringstream s;
  s  <<  key;
  return s.str();

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

 
//////////////////////////////////////////////////////////////////////
// Conversion of matrices
 	
/** \brief Convert a CLHEP HepMatrix to ROOT TMatrixD. */
void CLHEPtoROOT( HepMatrix & oldM, TMatrixD * newM)
{
  newM->ResizeTo(oldM.num_row(), oldM.num_col());
  for (int i = 0; i < oldM.num_row(); i++) {
    for (int j = 0; j < oldM.num_col(); j++) {
      (*newM)[i][j] = oldM[i][j];
    }
  }
}
 	
/** \brief Convert a CLHEP HepSymMatrix to ROOT TMatrixDSym. */
void CLHEPtoROOT( HepSymMatrix & oldM, TMatrixDSym * newM)
{
  newM->ResizeTo(oldM.num_row(), oldM.num_col());
  for (int i = 0; i < oldM.num_row(); i++) {
    for (int j = i; j < oldM.num_col(); j++) {
      (*newM)[i][j] = oldM[i][j];
    }
  }
}
 	
/** \brief Convert a CLHEP HepVector to ROOT TVectorD. */
void CLHEPtoROOT( HepVector & oldV, TVectorD * newV)
{
  newV->ResizeTo(oldV.num_row());
  for (int i = 0; i < oldV.num_row(); i++) {
    (*newV)[i] = oldV[i];
  }
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
 	
