#include "AlignableDet.h"

namespace depfet {	

/** Constructor 
 */
AlignableDet::AlignableDet(int nAlignables)
{ 
  
  // We assume six parameters per alignable for a planar sensor.
  // And we do not store the cross correlations between the alignment of different sensors. 
  
  for (int i=0;i<nAlignables;i++) {
    // Initially, no correction is known
    alignmentParametersVec.push_back(SensorAlignmentParameters::Zero());
    alignmentCovarianceVec.push_back(SensorAlignmentCovariance::Zero()); 
  }
}


/** Getter
 */
SensorAlignmentParameters AlignableDet::GetAlignState(int ipl) const
{
  return alignmentParametersVec.at(ipl); 
}

SensorAlignmentCovariance AlignableDet::GetAlignCovariance(int ipl) const
{
  return alignmentCovarianceVec.at(ipl); 
} 
 

 
/** Setter
 */
void AlignableDet::SetAlignState(int ipl, const SensorAlignmentParameters& a)
{
  alignmentParametersVec[ipl] = a; 
}

void AlignableDet::SetAlignCovariance(int ipl, const SensorAlignmentCovariance& E)
{
  alignmentCovarianceVec[ipl] = E; 
}  


} // Namespace;
