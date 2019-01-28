#ifndef TBVERTEX_H
#define TBVERTEX_H 1

// TBTools includes
#include "TBTrack.h"


#include <Eigen/Core>


typedef Eigen::Matrix<double,3,1> VertexParameter;
typedef Eigen::Matrix<double,3,3> VertexCovariance;
typedef Eigen::Matrix<double,5,3> VertexJacobian;
typedef Eigen::Matrix<double,4,1> VertexResidual;



namespace depfet { 


/** Class TBVertex
 *  
 *  The class TBVertex represents a vertex position
 *  relative to the scattering plane used for the 
 *  calculation of said vertex.
 *  It consists of the position r = (x,y,z), the
 *  corresponding covariance matrix C and 
 *  the Chi2-Value of the vertex fit.
 *  
 *  @Author E. Heinrich, University of Göttingen
 *  <mailto:erik.heinrich@stud.uni-goettingen.de>
 */

class TBVertex {

 public: //Members

  //Vertex position vector (in local and global coordinates)
  VertexParameter Pos;
  VertexParameter GlobalPos;
  
  //Corresponding covariance (in local and global coordinates)
  VertexCovariance Cov;
  VertexCovariance GlobalCov;
  //Vertex fit chi2
  double chi2;
  //Vertex fit number degrees of freedom
  int ndf;
  //Filter residual
  VertexResidual Res;

  //Vector containing all states for given measurement
  std::vector<TBTrackState> States;

 public: //Functions

  // Constructors
  TBVertex(); 
  TBVertex(VertexParameter aPos, VertexParameter aGlobalPos, VertexCovariance aCov, VertexCovariance aGlobalCov, double achi2);  
  
  // Get/Set local vertex position  
  void SetPos(const VertexParameter& aPos) { Pos= aPos; }; 
  VertexParameter&  GetPos() { return Pos; };

  // Get/Set global vertex position  
  void SetGlobalPos(const VertexParameter& aPos) { GlobalPos= aPos; }; 
  VertexParameter&  GetGlobalPos() { return GlobalPos; };
   
  // Get/Set local position covariance 
  void SetCov(const VertexCovariance& aCov ) { Cov = aCov; }; 
  VertexCovariance&  GetCov() { return Cov; };

  // Get/Set global position covariance 
  void SetGlobalCov(const VertexCovariance& aCov ) { GlobalCov = aCov; }; 
  VertexCovariance&  GetGlobalCov() { return GlobalCov; };

  // Get/Set chi²-value
  void SetChi2(double achi2 ) { chi2 = achi2; }; 
  double GetChi2() { return chi2; };
  double GetChi2Ndof() { return chi2/double(ndf); };

  // Get/Set filter residual
  void SetRes(const VertexResidual& aRes ) { Res = aRes; }; 
  VertexResidual& GetRes() { return Res; };

  // Get/Set ndf-value
  void SetNdf(int andf ) { ndf = andf; }; 
  int GetNdf() { return ndf; };

  // add state
  void AddTrackState(TBTrackState& state) { States.push_back(state); };

  // remove last state
  void RemoveLastTrackState() { States.pop_back(); };

  // get States vector
  std::vector<TBTrackState>& GetStates() { return States; };

};

} //Namespace

#endif
