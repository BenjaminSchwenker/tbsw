#ifndef _AlignableDet_h
#define _AlignableDet_h

// CLHEP includes 
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/SymMatrix.h>
 	
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
  AlignableDet(int nAlignables, int nParameters);
  
  /** Getter
   */
  CLHEP::HepVector GetAlignState(int i);
  CLHEP::HepSymMatrix GetAlignCovariance(int i); 
  
  /** Setter
   */
  void SetAlignState(int i, CLHEP::HepVector a);
  void SetAlignCovariance(int i, CLHEP::HepSymMatrix E);  
   
  /* Alignment parameters (all alignables) 
   */
  CLHEP::HepVector alignmentParameters;    

  /* Alignment covariance (all alignables)   
   */  
  CLHEP::HepSymMatrix alignmentCovariance;      
};

} // Namespace
#endif

