// TBTools includes
#include "TBVertex.h"


namespace depfet {

/**
 * File TBVertex.cc
 * 
 * Definition of class TBVertex.
 * 
 *  @Author E. Heinrich, University of Göttingen
 *  <mailto:erik.heinrich@stud.uni-goettingen.de>
 */

// Constructors

TBVertex::TBVertex() 
{

 Pos = VertexParameter::Zero();
 GlobalPos = VertexParameter::Zero();
  
 Cov = VertexCovariance::Zero();
 GlobalCov = VertexCovariance::Zero();
  
 chi2 = 0;
 ndf = -3;
  
 Res = VertexResidual::Zero();
}

TBVertex::TBVertex(VertexParameter aPos, VertexParameter aGlobalPos, VertexCovariance aCov, VertexCovariance aGlobalCov, double achi2) 
{ 
  // Set vertex position 
  Pos = aPos;
  // Set global vertex position 
  GlobalPos = aGlobalPos;
  // Set vertex covariance
  Cov = aCov;
  // Set global vertex covariance
  GlobalCov = aGlobalCov;
  // Set fit chi2
  chi2 = achi2;

  ndf = -3;
  Res = VertexResidual::Zero();
}

} // Namespace
