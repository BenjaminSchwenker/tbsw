#include "AlignableDet.h"


using namespace CLHEP;

namespace depfet {	

/** Constructor 
 */
AlignableDet::AlignableDet(int nAlignables, int nParameters)
{ 
  // Initially, no correction is known
  alignmentParameters = HepVector(nAlignables*nParameters,0);      
  alignmentCovariance = HepSymMatrix(nAlignables*nParameters,0);  
}


/** Getter
 */
HepVector AlignableDet::GetAlignState(int ipl) 
{
  HepVector a(6,0);     
  for (int k=0; k<6; ++k) { 
    a[k] = alignmentParameters[ipl*6+k];
  }
  return a; 
}


HepSymMatrix AlignableDet::GetAlignCovariance(int ipl)
{
  HepSymMatrix E(6,0); 
  for (int k=0; k<6; ++k) {
    for (int l=0; l<=k; ++l) {    
      E[k][l] = alignmentCovariance[ipl*6+k][ipl*6+l];
    }
  }
  return E; 
} 
 

 
/** Setter
 */
void AlignableDet::SetAlignState(int ipl, HepVector a)
{
  for (int k=0; k<6; ++k) { 
    alignmentParameters[ipl*6+k] = a[k];
  }
}


void AlignableDet::SetAlignCovariance(int ipl, HepSymMatrix E)
{
   for (int k=0; k<6; ++k) {
    for (int l=0; l<=k; ++l) {    
      alignmentCovariance[ipl*6+k][ipl*6+l] = E[k][l];
    }
  }
}  


} // Namespace;
