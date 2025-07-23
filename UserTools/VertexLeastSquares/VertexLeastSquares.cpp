#include "VertexLeastSquares.h"

VertexLeastSquares::VertexLeastSquares():Tool(){}

//------------------------------------------------------------------------------
bool VertexLeastSquares::Initialise(std::string configfile, DataModel &data)
{

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // Load my config parameters
  bool gotVerbosity = m_variables.Get("verbosity",verbosity);
  if (!gotVerbosity) {
    verbosity = 0;
    Log("VertexLeastSquares: \"verbosity\" not set in the config, defaulting to 0", v_error, verbosity);
  }

  //bool gotClusterMapName = m_variables.Get("ClusterMapName", fClusterMapName);
 // if (!gotClusterMapName) {
    //Log("VertexLeastSquares: \"ClusterMapName\" not set in the config! Aborting!", v_error, verbosity);
    //return false;
  //}

  bool gotUseMCHits = m_variables.Get("UseMCHits", fUseMCHits);
  if (!gotUseMCHits) {
    Log("VertexLeastSquares: \"UseMCHits\" not set in the config! Aborting!", v_error, verbosity);
    return false;
  }

  bool gotBreakDist = m_variables.Get("BreakDist", fBreakDist);
  if (!gotBreakDist) {
    fBreakDist = 0.005;
    Log("VertexLeastSquares: \"BreakDist\" not set in the config! Using default value of 0.5 cm.", v_error, verbosity);
  }

  bool gotYSpacing = m_variables.Get("YSpacing", fYSpacing);
  if (!gotYSpacing) {
    fYSpacing = 0.5;
    Log("VertexLeastSquares: \"YSpacing\" not set in the config! Using default value of 50 cm.", v_error, verbosity);
  }

  bool gotNPlanarPoints = m_variables.Get("NPlanarPoints", fNPlanarPoints);
  if (!gotNPlanarPoints) {
    fNPlanarPoints = 20;
    Log("VertexLeastSquares: \"NPlanarPoints\" not set in the config! Using default value of 20.", v_error, verbosity);
  }

  bool gotRegularizer = m_variables.Get("Regularizer", fRegularizer);
  if (!gotRegularizer) {
    fRegularizer = 1.;
    Log("VertexLeastSquares: \"Regularizer\" not set in the config! Using default value of 1.", v_error, verbosity);
  }

  fNBoundary = int(sqrt(fRegularizer));
  
  // Load the geometry service
  bool gotGeometry = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",fGeom);
  if(!gotGeometry){
    Log("VertexLeastSquares: Error retrieving Geometry from ANNIEEvent! Aborting!", v_error, verbosity); 
    return false; 
  }

  // Set up the pointer we're going to save. No need to 
  // delete it at Finalize, the store will handle it
  fVertexMap = new std::map<double, Position>;

  return true;
}

//------------------------------------------------------------------------------
bool VertexLeastSquares::Execute()
{

  // An upstream tool may opt to skip this execution stage
  // For example the PMTWaveformSim tool will skip events with no MCHits or if
  // no waveforms are produced.
  bool skip = false;
  bool got_skip_status = m_data->Stores["ANNIEEvent"]->Get("SkipExecute", skip);
  if (got_skip_status && skip) {
    Log("VertexLeastSquares: An upstream tool told me to skip this event.",v_warning,verbosity);
    return true;
  } 

  fVertexMap->clear();
  if (fUseMCHits) {
    bool gotClusters = m_data->CStore.Get("ClusterMapMC", fClusterMapMC);
    if (!gotClusters) {
      logmessage = "VertexLeastSquares: no ClusterMapMC in the ANNIEEvent!!! Are you sure you ran ClusterFinder?";
      Log(logmessage, v_error, verbosity);
      return false;
    }
    RunLoopMC();
  } else {
    bool gotClusters = m_data->CStore.Get("ClusterMap", fClusterMap);
    if (!gotClusters) {
      logmessage = "VertexLeastSquares: no ClusterMap in the ANNIEEvent!!! Are you sure you ran ClusterFinder?";
      Log(logmessage, v_error, verbosity);
      return false;
    }

    RunLoop();
  }

  m_data->Stores.at("ANNIEEvent")->Set("VertexLeastSquaresMap",  fVertexMap);
  
  return true;
}

//------------------------------------------------------------------------------
bool VertexLeastSquares::Finalise()
{

  return true;
}

//-----------------------------------------------------------------------------
std::vector<Position> VertexLeastSquares::GenerateVetices()
{
  // generate the y values
  std::vector<double> ys;
  double yMin = fGeom->GetTankCentre().Y() - fGeom->GetTankHalfheight();
  double yMax = fGeom->GetTankCentre().Y() + fGeom->GetTankHalfheight();
  for (double yValue = yMin; yValue <= yMax; yValue += fYSpacing) 
    ys.push_back(yValue);
  
  // generate the x and z values using the sunflower pattern
  std::vector<double> xs;
  std::vector<double> zs;
  for (int n = 1; n < fNPlanarPoints+1; ++n) {
    double rad = ( (n > fNPlanarPoints + fNBoundary) ? 1.0 :
		   sqrt((n+0.5)/(fNPlanarPoints - (fNBoundary+1.)/2.)) );

    // Scale the radius to the tank and make sure it's inside
    rad *= fGeom->GetTankRadius();
    rad = rad > fGeom->GetTankRadius() ? fGeom->GetTankRadius() : rad;

    double angle = 2. * pi * n / phisq;

    xs.push_back(rad*cos(angle));
    zs.push_back(rad*sin(angle));
  }

  std::vector<Position> vertices;
  for (uint yIdx = 0; yIdx < ys.size(); ++yIdx) {

    // Apply a random rotation to each y level to cover more space
    double rotAng = (double(rand())/RAND_MAX)*2*pi;
    for (uint xzIdx = 0; xzIdx < fNPlanarPoints; ++xzIdx) {
      double x = xs[xzIdx];
      double z = zs[xzIdx];

      vertices.push_back(Position(x*cos(rotAng) - z*sin(rotAng) + fGeom->GetTankCentre().X(),
				  ys[yIdx],
				  z*cos(rotAng) + x*sin(rotAng) + fGeom->GetTankCentre().Z()));
    }
  }

  return vertices;
}

//------------------------------------------------------------------------------
void VertexLeastSquares::EvalAtGuessVertex(util::Matrix &A, util::Vector &b,
					   const Position &guess,
					   const std::vector<Hit> &hits)
{
  // The function f(x) = 0
  // 0 = 1/c * sqrt( (X-xi)^2 + (Y-yi)^2 + (Z-zi)^2 ) + T - ti
  // A = Jacobian(guess), b = -f(guess)
  double mean_T = 0;
  for (uint hitIdx = 0; hitIdx < hits.size(); ++hitIdx) {
    Detector *det = fGeom->ChannelToDetector(hits[hitIdx].GetTubeId());
    Position detpos = det->GetDetectorPosition();
    double dist = (detpos - guess).Mag();    
    double ti = hits[hitIdx].GetTime();
    
    A(hitIdx, 0) = (guess.X() - detpos.X())/(dist * fSoL);
    A(hitIdx, 1) = (guess.Y() - detpos.Y())/(dist * fSoL);
    A(hitIdx, 2) = (guess.Z() - detpos.Z())/(dist * fSoL);
    A(hitIdx, 3) = 1.;

    b(hitIdx) = -dist/fSoL + ti;
    
    mean_T += ti - dist/fSoL;
  }// end loop over hits
  
  mean_T = mean_T / hits.size();

  // Subtract off the mean emission time (T)
  for (uint hitIdx = 0; hitIdx < hits.size(); ++hitIdx) 
    b(hitIdx) -= mean_T;  

  // Tack on the regularization term
  A(hits.size(),     0) = fRegularizer;
  A(hits.size() + 1, 1) = fRegularizer;
  A(hits.size() + 2, 2) = fRegularizer;
  A(hits.size() + 3, 3) = fRegularizer;
  b(hits.size()    ) = 0;
  b(hits.size() + 1) = 0;
  b(hits.size() + 2) = 0;
  b(hits.size() + 3) = 0;


}


//------------------------------------------------------------------------------
void VertexLeastSquares::EvalAtGuessVertexMC(util::Matrix &A, util::Vector &b,
					   const Position &guess,
					   const std::vector<MCHit> &hits)
{
  // The function f(x) = 0
  // 0 = 1/c * sqrt( (X-xi)^2 + (Y-yi)^2 + (Z-zi)^2 ) + T - ti
  // A = Jacobian(guess), b = -f(guess)
  double mean_T = 0;
  for (uint hitIdx = 0; hitIdx < hits.size(); ++hitIdx) {
    Detector *det = fGeom->ChannelToDetector(hits[hitIdx].GetTubeId());
    Position detpos = det->GetDetectorPosition();
    double dist = (detpos - guess).Mag();    
    double ti = hits[hitIdx].GetTime();
    
    A(hitIdx, 0) = (guess.X() - detpos.X())/(dist * fSoL);
    A(hitIdx, 1) = (guess.Y() - detpos.Y())/(dist * fSoL);
    A(hitIdx, 2) = (guess.Z() - detpos.Z())/(dist * fSoL);
    A(hitIdx, 3) = 1.;

    b(hitIdx) = -dist/fSoL + ti;
    
    mean_T += ti - dist/fSoL;
  }// end loop over hits
  
  mean_T = mean_T / hits.size();

  // Subtract off the mean emission time (T)
  for (uint hitIdx = 0; hitIdx < hits.size(); ++hitIdx) 
    b(hitIdx) -= mean_T;

  // Tack on the regularization term
  A(hits.size(),     0) = fRegularizer;
  A(hits.size() + 1, 1) = fRegularizer;
  A(hits.size() + 2, 2) = fRegularizer;
  A(hits.size() + 3, 3) = fRegularizer;
  b(hits.size()    ) = 0;
  b(hits.size() + 1) = 0;
  b(hits.size() + 2) = 0;
  b(hits.size() + 3) = 0;


}

//------------------------------------------------------------------------------
void VertexLeastSquares::RunLoop()
{

  for (auto clusterpair : *fClusterMap) {    
    // Filter hits
    auto filt_hits = FilterHits(clusterpair.second);

    // We need at least 4 hits
    if (filt_hits.size() < 4) {
      std::cout << "Not enough hits. We only have: " << filt_hits.size() << std::endl;
      fVertexMap->emplace(clusterpair.first, Position(-9, -9, -9));
      continue;
    }
    
    // Loop over vertex seeds
    std::vector<Position> seeds = GenerateVetices();

    std::vector<Position> bestVtxs;
    double bestStDev = 999999;
    int bestIdx = 0;
    int seednum = 0;
    for (auto guess : seeds) {
      Position lastGuess;
      int loop_num = 0;
      bool inTank = true;
      auto hits = filt_hits;
      while (true) {
	lastGuess = guess;
  inTank = true;
	
	// Set up Matrix and solution vector
	util::Matrix A(filt_hits.size() + 4, 4);
	util::Vector b(filt_hits.size() + 4);
	EvalAtGuessVertex(A, b, guess, hits);
	
	// Solve it and update the guess vertex
	util::Vector solution(4);
	util::least_squares(A, b, solution);

	// If the solution has nans it's because all rows of A evaluate to the same thing
	// revert to the last vertex and break
	if (std::isnan(solution(0)) || std::isnan(solution(1)) || std::isnan(solution(2))) {
	  inTank = false;
	  break;

	}

  // Otherwise, update the guess and carry on
	guess = lastGuess + Position(solution(0), solution(1), solution(2));
	
	// Exit if the new guess is close enough to the last one
	if ((guess - lastGuess).Mag() < fBreakDist) break;
	
	// Exit if we're outside of the tank
	if(!fGeom->GetTankContained(guess)) {
	  inTank = false;
	  break;
	}

	// Don't let the loop last forever
	++loop_num;
	if (loop_num > 100) break;
      }// end while loop

      // Skip if no tank contained vertex was found for this seed
      if (!inTank) continue;
      
      bestVtxs.push_back(guess);
      
      // Check the std deviation of b for the best guess from this seed
      util::Matrix A(filt_hits.size() + 4, 4);
      util::Vector b(filt_hits.size() + 4);
      EvalAtGuessVertex(A, b, bestVtxs.back(), filt_hits);
      double sum = 0;
      double sum_sq = 0;
      for (uint idx = 0; idx < b.size - 4; ++idx) {
	sum += b(idx);
	sum_sq += b(idx) * b(idx) ;
      }
      double mean = sum / b.size;
      double stdev = sum_sq / b.size - mean * mean;
      if (stdev < bestStDev) {
	bestStDev = stdev;
	bestIdx = bestVtxs.size()-1;
      }
      ++seednum;
    }// end loop over seed vertices

    // Save the vertex if we have one, otherwise default to 9s
    if (bestVtxs.size())
      fVertexMap->emplace(clusterpair.first, bestVtxs[bestIdx]);
    else
      fVertexMap->emplace(clusterpair.first, Position(-9, -9, -9));
  }// end loop over the clusters
}

//------------------------------------------------------------------------------
void VertexLeastSquares::RunLoopMC()
{
  for (auto clusterpair : *fClusterMapMC) {
    // Filter hits
    auto filt_hits = FilterHitsMC(clusterpair.second);

    // We need at least 4 hits
    if (filt_hits.size() < 4) {
      std::cout << "Not enough hits. We only have: " << filt_hits.size() << std::endl;
      fVertexMap->emplace(clusterpair.first, Position(-9, -9, -9));
      continue;
    }

    // Loop over vertex seeds
    std::vector<Position> seeds = GenerateVetices();

    std::vector<Position> bestVtxs;
    double bestStDev = 999999;
    int bestIdx = 0;
    int seednum = 0;
    for (auto guess : seeds) {
      Position lastGuess;
      int loop_num = 0;
      bool inTank = true;
      auto hits = filt_hits;
      while (true) {
	lastGuess = guess;
	inTank = true;
	
	// Set up Matrix and solution vector
	util::Matrix A(hits.size() + 4, 4);
	util::Vector b(hits.size() + 4);
	EvalAtGuessVertexMC(A, b, guess, hits);
	
	// Solve it and update the guess vertex
	util::Vector solution(4);
	util::least_squares(A, b, solution);

	// If the solution has nans it's because all rows of A evaluate to the same thing
	// revert to the last vertex and break
	if (std::isnan(solution(0)) || std::isnan(solution(1)) || std::isnan(solution(2))) {
	  inTank = false;
	  break;

	}

	// Otherwise, update the guess and carry on
	guess = lastGuess + Position(solution(0), solution(1), solution(2));
	
	// Exit if the new guess is close enough to the last one
	if ((guess - lastGuess).Mag() < fBreakDist) break;

	// If we're outside of the tank then remove the last hit and try again
	// but if there are only 4 hits left then just give up
	if(!fGeom->GetTankContained(guess)) {
	    inTank = false;
	    break;
	}

	// Don't let the loop last forever
	++loop_num;
	if (loop_num > 100) break;
      }// end while loop

      // Skip if no tank contained vertex was found for this seed
      if (!inTank) continue;
      
      // Check the std deviation of b for the best guess from this seed
      bestVtxs.push_back(guess);
      util::Matrix A(filt_hits.size() + 4, 4);
      util::Vector b(filt_hits.size() + 4);
      EvalAtGuessVertexMC(A, b, bestVtxs.back(), filt_hits);
      double sum = 0;
      double sum_sq = 0;
      for (uint idx = 0; idx < b.size - 4; ++idx) {
	sum += b(idx);
	sum_sq += b(idx) * b(idx) ;
      }
      double mean = sum / b.size;
      double stdev = sum_sq / b.size - mean * mean;
      if (stdev < bestStDev) {
	bestStDev = stdev;
	bestIdx = bestVtxs.size()-1;
      }
      ++seednum;
    }// end loop over seed vertices

    // Save the vertex if we have one, otherwise default to -9s
    if (bestVtxs.size()) {      
      fVertexMap->emplace(clusterpair.first, bestVtxs[bestIdx]);
    }
    else {
      std::cout << "No good vertex found!!" << std::endl;
      fVertexMap->emplace(clusterpair.first, Position(-9, -9, -9));
    }
  }// end loop over the clusters
}

//------------------------------------------------------------------------------
std::vector<Hit> VertexLeastSquares::FilterHits(std::vector<Hit> hits)
{ 
  // Sort the hits by time
  std::sort(hits.begin(), hits.end(),
	    [](const Hit &a, const Hit &b) {return a.GetTime() < b.GetTime();});

  // Check for causally independent hits
  // time difference between the hits must be less than the relative time of flight between the PMTs
  std::map<int, std::set<int>> compatibleIds;
  int mostCompatibleIds = 0;
  int bestIdx = 0;
  for (uint idx = 0; idx < hits.size(); ++idx) {
    Detector *det_i = fGeom->ChannelToDetector(hits[idx].GetTubeId());
    Position pos_i = det_i->GetDetectorPosition();

    for (uint jdx = 0; jdx < hits.size(); ++jdx) {
      Detector *det_j = fGeom->ChannelToDetector(hits[jdx].GetTubeId());
      Position pos_j = det_j->GetDetectorPosition();
      
      double deltaT = abs(hits[jdx].GetTime() - hits[idx].GetTime());
      double tof = (pos_i - pos_j).Mag() / fSoL;
      
      if (deltaT <= tof) 
	compatibleIds[idx].insert(jdx);
    }

    if (compatibleIds[idx].size() > mostCompatibleIds) {
      mostCompatibleIds = compatibleIds[idx].size();
      bestIdx = idx;
    }
  }

  // Keep the set of most compatible ids if they're within 10 ns of the first hit
  std::vector<Hit> filtered;
  int firstHitIdx = *compatibleIds[bestIdx].begin();
  double t0 = hits[firstHitIdx].GetTime();
  for (auto idx : compatibleIds[bestIdx]) {
    if (abs(hits[idx].GetTime() - t0) > 13) continue;
    filtered.push_back(hits[idx]);
  }

  return filtered;
}


//------------------------------------------------------------------------------
std::vector<MCHit> VertexLeastSquares::FilterHitsMC(std::vector<MCHit> hits)
{ 
  
  std::vector<MCHit> filtered;

  // Sort the hits by time
  std::sort(hits.begin(), hits.end(),
	    [](const Hit &a, const Hit &b) {return a.GetTime() < b.GetTime();});

  // Keep the first four and find the mean time
  double mean = 0;
  for (uint idx = 0; idx < 4; ++idx) {
    mean += hits[idx].GetTime();
    filtered.push_back(hits[idx]);
  }
  mean = mean/4.;

  // Look through the rest of the hits,
  // check that they are within 10 ns (for reflections),
  // and check that they could come from the same point
  for (uint idx = 4; idx < hits.size(); ++idx) {
    if ((hits[idx].GetTime() - mean) > 10)
      continue;

    bool keep = true;
    for (uint jdx = 0; jdx < 4; ++jdx) {
      Detector *det_i = fGeom->ChannelToDetector(hits[idx].GetTubeId());
      Position pos_i = det_i->GetDetectorPosition();

      Detector *det_j = fGeom->ChannelToDetector(filtered[jdx].GetTubeId());
      Position pos_j = det_j->GetDetectorPosition();

      double deltaT = abs(filtered[jdx].GetTime() - hits[idx].GetTime());
      double tof = (pos_i - pos_j).Mag() / fSoL;

      if (deltaT > tof) {
	keep = false;
	break;
      }
    }

    if (!keep) continue;
    filtered.push_back(hits[idx]);
  }

  return filtered;
}
