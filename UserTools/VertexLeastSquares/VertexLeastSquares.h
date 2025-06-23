#ifndef VertexLeastSquares_H
#define VertexLeastSquares_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MatrixUtil.h"
#include "Position.h"
#include "Hit.h"

#include "TFile.h"
#include "TTree.h"



/**
 * \class VertexLeastSquares
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/
class VertexLeastSquares: public Tool {


 public:

  VertexLeastSquares(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  // Generate a cylinder with XYZ nodes in each dimension
  std::vector<Position> GenerateVetices();

  // Evaluate the Jacobian and the solution vector at the guess vertex
  //   we'll be solving A*dx = b where A = Jacobian(guess), b = -f(guess),
  //   f = vector of the hit time functions: 0 = 1/c * sqrt( (X-xi)^2 + (Y-yi)^2 + (Z-zi)^2 ) + T - ti
  //   X,Y,Z are the coords of the vertex and T is the emission time of the light
  //   xi,yi,zi are the PMT coords and ti is the PMT hit time 
  void EvalAtGuessVertex(util::Matrix &A, util::Vector &b, const Position &guess, const std::vector<Hit> &hits);
  void EvalAtGuessVertexMC(util::Matrix &A, util::Vector &b, const Position &guess, const std::vector<MCHit> &hits);

  // Actually run the thing until convergence
  void RunLoop();
  void RunLoopMC();

  // Filter hits based on time separation and causality
  std::vector<Hit> FilterHits(std::vector<Hit> hits);
  std::vector<MCHit> FilterHitsMC(std::vector<MCHit> hits);


 private:

  // Configuration parameters
  //std::string fClusterMapName;
  bool fUseMCHits;
  double fBreakDist;
  bool fDebugTree;
  double fYSpacing;
  int fNPlanarPoints;
  double fRegularizer;

  // The ANNIE geometry service
  Geometry *fGeom = nullptr;

  // The clusters we'll load from the CStore
  std::map<double, std::vector<Hit>>   *fClusterMap   = nullptr;
  std::map<double, std::vector<MCHit>> *fClusterMapMC = nullptr;

  // Vertices we'll save to the ANNIEEvent
  std::map<double, Position> *fVertexMap = nullptr;

  // The number of boundary points for the XZ plane verticies
  int fNBoundary;

  // Golden ratio and pi
  const double phi = (1. + sqrt(5.))/2.;
  const double phisq = pow(phi, 2);
  const double pi = 3.141592;

  // Speed of light in water, m/ns
  const double fSoL = 0.299792458 * 3/4; 

  /// \brief verbosity levels: if 'verbosity' < this level, the message type will be logged.
  int verbosity;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;
};


#endif
