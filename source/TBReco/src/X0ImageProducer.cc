// X0ImageProduce
//
// See X0ImageProducer.h for full documentation of processor.
//
// Author: Benjamin Schwenker, Göttingen University
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "X0ImageProducer.h"

// TBTools includes
#include "MaterialEffect.h"
#include "TBDetector.h"
#include "TBHit.h"
#include "TBKalmanMSC.h"
#include "TBTrack.h"
#include "TBVertex.h"
#include "TBVertexFitter.h"
#include "TrackInputProvider.h"
#include "Utilities.h"

// Include basic C
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

// Include LCIO classes
#include <IMPL/LCFlagImpl.h>
#include <UTIL/CellIDDecoder.h>

// ROOT includes
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>

// Used namespaces
using namespace std;
using namespace lcio;
using namespace marlin;

namespace depfet {

//
// Instantiate this object
//
X0ImageProducer aX0ImageProducer;

//
// Constructor
//
X0ImageProducer::X0ImageProducer() : Processor("X0ImageProducer") {

  // Processor description
  _description =
      "X0ImageProducer: X/X0 measurement for EUDET/AIDA telescope data";

  //
  // Input collections

  registerInputCollection(LCIO::TRACK, "DownStreamTrackCollection",
                          "Name of downstream track collection",
                          _downStreamTrackColName, std::string("down_tracks"));

  registerInputCollection(LCIO::TRACK, "UpStreamTrackCollection",
                          "Name of upstream track collection",
                          _upStreamTrackColName, std::string("up_tracks"));

  // Processor parameters:

  registerProcessorParameter("DUTPlane",
                             "Plane number of DUT along the beam line", _idut,
                             static_cast<int>(3));

  registerProcessorParameter("VertexFitSwitch",
                             "Choose upstream-downstream track matching "
                             "approach. true: vertexfit, false: distance cut",
                             _vertexfitswitch, static_cast<bool>(false));

  registerProcessorParameter("ToyScatteringSwitch",
                             "Switch to fast toy simulation mode: Replace "
                             "reconstructed scatter angles from track fitting "
                             "by angles sampled from a Gauss distribution",
                             _m_toy, static_cast<bool>(false));

  registerProcessorParameter(
      "ToyScatterModel", "Choose model for multiple scattering used in toy "
                         "simulations: Highland(0), SumOfSingleScatterings(1)",
      _m_toyScatterModel, static_cast<int>(0));

  registerProcessorParameter(
      "ToyRecoError", "Angle reconstruction error used in toy simulations, if "
                      "its smaller than 0, use real angle reco error instead "
                      "(only used in toy simulation mode)",
      _m_reco_error, static_cast<double>(-1.0));

  registerProcessorParameter("ToyBetheHeitlerSwitch",
                             "Flag (true/false) for simulating fractional "
                             "Bethe Heitler energy loss (only used in toy "
                             "simulation mode)",
                             _m_ToyBetheHeitler, static_cast<bool>(true));

  registerProcessorParameter("RootFileName", "Output root file name",
                             _rootFileName, std::string("X0.root"));

  registerProcessorParameter("MaxVertexChi2", "Maximal Chi2/Ndof value for "
                                              "vertex fit of matched up and "
                                              "downstream tracks, only used "
                                              "when VertexFitSwitch is true",
                             _maxVertexChi2, static_cast<double>(10.0));

  registerProcessorParameter(
      "MaxDist", "Maximum distance between up and downstream tracks at dut "
                 "plane [mm], only used when VertexFitSwitch is false",
      _maxDist, static_cast<double>(0.1));
}

//
// Method called at the beginning of data processing
//
void X0ImageProducer::init() {

  // Initialize variables
  _nRun = 0;
  _nEvt = 0;
  _timeCPU = clock() / 1000;

  // Print set parameters
  printProcessorParams();

  // Load DUT module
  const Det &dut = TBDetector::GetInstance().GetDet(_idut);

  // Print out geometry information
  streamlog_out(MESSAGE3) << "Scatter DUT plane  ID = " << dut.GetSensorID()
                          << "  at position = " << _idut << endl
                          << endl;

  // Initialize new random generator with unique seed
  gRandom = new TRandom3();
  gRandom->SetSeed(0);

  // Print random generator info
  streamlog_out(MESSAGE3) << "Random Generator setup with seed: "
                          << (gRandom->GetSeed()) << std::endl
                          << std::endl;

  bookHistos();
}

//
// Method called for each run
//
void X0ImageProducer::processRunHeader(LCRunHeader *run) {
  // Print run number
  streamlog_out(MESSAGE3) << "Processing run: " << (run->getRunNumber())
                          << std::endl
                          << std::endl;

  _nRun++;
}

//
// Method called for each event
//
void X0ImageProducer::processEvent(LCEvent *evt) {

  // Print event number
  if ((evt->getEventNumber()) % 100 == 0)
    streamlog_out(MESSAGE3) << "Events processed: " << (evt->getEventNumber())
                            << std::endl
                            << std::endl;

  streamlog_out(MESSAGE2) << "Events processed: " << (evt->getEventNumber())
                          << std::endl
                          << std::endl;

  _nEvt++;

  //
  // Get telescope track collection
  //

  LCCollection *downtrackcol = 0;
  try {
    downtrackcol = evt->getCollection(_downStreamTrackColName);
  } catch (lcio::DataNotAvailableException &e) {
    streamlog_out(MESSAGE2) << "Not able to get collection "
                            << _downStreamTrackColName << " from event "
                            << evt->getEventNumber() << " in run "
                            << evt->getRunNumber() << endl
                            << endl;

    throw SkipEventException(this);
  }

  LCCollection *uptrackcol = 0;
  try {
    uptrackcol = evt->getCollection(_upStreamTrackColName);
  } catch (lcio::DataNotAvailableException &e) {
    streamlog_out(MESSAGE2) << "Not able to get collection "
                            << _upStreamTrackColName << " from event "
                            << evt->getEventNumber() << " in run "
                            << evt->getRunNumber() << endl
                            << endl;

    throw SkipEventException(this);
  }

  // Configure Kalman track fitter
  TBKalmanMSC TrackFitterMSC(TBDetector::GetInstance());
  TrackInputProvider TrackLCIOReader;

  // Store tracks
  std::vector<TBTrack> downTrackStore;

  // Loop over tracks in input track collection

  int nDownTracks = downtrackcol->getNumberOfElements();
  for (int itrk = 0; itrk < nDownTracks; itrk++) {

    // Retrieve track from LCIO
    Track *lciotrk = dynamic_cast<Track *>(downtrackcol->getElementAt(itrk));

    // Convert LCIO -> TB track
    TBTrack trk =
        TrackLCIOReader.MakeTBTrack(lciotrk, TBDetector::GetInstance());

    // Refit track in nominal alignment
    bool trkerr = TrackFitterMSC.ProcessTrack(trk, -1, 0);
    if (trkerr) {
      streamlog_out(MESSAGE1) << "Fit failed. Skipping track!" << endl;
      continue;
    }

    downTrackStore.push_back(trk);
  }

  // Store tracks
  std::vector<TBTrack> upTrackStore;

  // Loop over tracks in input track collection

  int nUpTracks = uptrackcol->getNumberOfElements();
  for (int itrk = 0; itrk < nUpTracks; itrk++) {

    // Retrieve track from LCIO
    Track *lciotrk = dynamic_cast<Track *>(uptrackcol->getElementAt(itrk));

    // Convert LCIO -> TB track
    TBTrack trk =
        TrackLCIOReader.MakeTBTrack(lciotrk, TBDetector::GetInstance());

    // Refit track in nominal alignment
    bool trkerr = TrackFitterMSC.ProcessTrack(trk, 1, 0);
    if (trkerr) {
      streamlog_out(MESSAGE1) << "Fit failed. Skipping track!" << endl;
      continue;
    }

    upTrackStore.push_back(trk);
  }

  //
  // Match upstream and downstream tracks at DUT plane
  int nMatch = 0;
  vector<vector<int>> up2down(upTrackStore.size());
  vector<vector<int>> down2up(downTrackStore.size());

  Det &dut = TBDetector::GetInstance().GetDet(_idut);

  // Initialize Vertex Fitter
  TBVertexFitter VertexFitter(_idut);

  // Vertex fit track matching
  if (_vertexfitswitch) {
    // Continue matching tracks and hits until all tracks are matched
    // or no hit is close enough to a track!!
    double chi2min = numeric_limits<double>::max();

    do {
      int bestup = -1;
      int bestdown = -1;

      chi2min = numeric_limits<double>::max();

      for (size_t iup = 0; iup < upTrackStore.size(); iup++) {
        // Get upstream track
        TBTrack &uptrack = upTrackStore[iup];

        for (size_t idown = 0; idown < downTrackStore.size(); idown++) {

          // If downtrack matched to some uptrack, skip downtrack
          // Downtracks should only be used once.
          if (down2up[idown].size() > 0)
            continue;

          // Get upstream track
          TBTrack &downtrack = downTrackStore[idown];

          // Initialize Vertex
          TBVertex Vertex;

          // Add upstream track state to vertex
          Vertex.AddTrackState(uptrack.GetTE(_idut).GetState());

          // Additionally add any downstream tracks, which have already been
          // sucessfully matched with this upstream track
          // This ensures that all matched downstream tracks really come from
          // the same vertex
          for (size_t n = 0; n < up2down[iup].size(); n++)
            Vertex.AddTrackState(
                downTrackStore[up2down[iup][n]].GetTE(_idut).GetState());

          // Add current downstream track state to vertex
          Vertex.AddTrackState(downtrack.GetTE(_idut).GetState());

          // Perform vertex fit
          // Check whether the combination of all
          bool vfiterr = VertexFitter.FitVertex(Vertex);
          if (vfiterr)
            continue;
          double vertex_chi2 = Vertex.GetChi2Ndof();

          if (vertex_chi2 < chi2min) {
            chi2min = vertex_chi2;
            bestup = iup;
            bestdown = idown;
          }
        }
      }

      streamlog_out(MESSAGE2) << "In matching loop: best up " << bestup
                              << " to best down " << bestdown << endl;
      streamlog_out(MESSAGE2) << "  chi2min: " << chi2min << endl;

      // Check if best match is good enough
      if (chi2min < _maxVertexChi2) {

        streamlog_out(MESSAGE2) << "  match found!!!" << endl;
        nMatch++;
        up2down[bestup].push_back(bestdown);
        down2up[bestdown].push_back(bestup);
      }

    } // End matching loop
    while (chi2min < _maxVertexChi2);
  } else {
    // Continue matching tracks and hits until all tracks are matched
    // or no hit is close enough to a track!!
    double distmin = numeric_limits<double>::max();

    do {
      int bestup = -1;
      int bestdown = -1;

      distmin = numeric_limits<double>::max();

      for (size_t iup = 0; iup < upTrackStore.size(); iup++) {

        // if ( up2down[iup].size() > 0 ) continue;

        for (int idown = 0; idown < (int)downTrackStore.size(); idown++) {

          // If matched, skip track
          if (down2up[idown].size() > 0)
            continue;

          TBTrack &uptrack = upTrackStore[iup];
          TBTrack &downtrack = downTrackStore[idown];

          // In and OutStates of the reconstructed Track at the current detector
          TBTrackState &InState = uptrack.GetTE(_idut).GetState();
          TBTrackState &OutState = downtrack.GetTE(_idut).GetState();

          double u_in = InState.GetPars()[2];
          double v_in = InState.GetPars()[3];
          double u_out = OutState.GetPars()[2];
          double v_out = OutState.GetPars()[3];

          double hitdist = std::abs(u_in - u_out) + std::abs(v_in - v_out);

          if (hitdist < distmin) {
            distmin = hitdist;
            bestup = iup;
            bestdown = idown;
          }
        }
      }

      streamlog_out(MESSAGE2) << "In matching loop: best up " << bestup
                              << " to best down " << bestdown << endl;
      streamlog_out(MESSAGE2) << "  distmin: " << distmin << endl;

      // Check if best match is good enough
      if (distmin < _maxDist) {

        streamlog_out(MESSAGE2) << "  match found!!!" << endl;
        nMatch++;
        up2down[bestup].push_back(bestdown);
        down2up[bestdown].push_back(bestup);
      }

    } // End matching loop
    while (distmin < _maxDist);
  }

  // Fill event tree
  _rootRunNumber = evt->getRunNumber();
  _rootEventNumber = evt->getEventNumber();
  _rootnUpTracks = upTrackStore.size();
  _rootnDownTracks = downTrackStore.size();
  _rootNMatched = nMatch;
  _rootFile->cd("");
  _rootEventTree->Fill();

  // Average momentum before and after energy loss.
  double weighted_mean_mom = 0;

  for (size_t iup = 0; iup < upTrackStore.size(); iup++) {

    // Check upstream track is matched
    if (up2down[iup].size() < 1)
      continue;

    TBTrack &uptrack = upTrackStore[iup];

    // Scattering Vertex fitting
    // The In and every Out State given for one vertex is added to a vertex
    // class
    TBVertex Vertex;
    Vertex.AddTrackState(uptrack.GetTE(_idut).GetState());

    // Add downstream trackstates to vertex
    for (size_t idown = 0; idown < up2down[iup].size(); idown++) {
      Vertex.AddTrackState(
          downTrackStore[up2down[iup][idown]].GetTE(_idut).GetState());
    }

    // Calculate vertex parameters
    // In case the vertex multiplicity is larger than 1, these values will be
    // set for every upstream downstream track combination
    bool vfiterr = VertexFitter.FitVertex(Vertex);
    if (vfiterr) {
      streamlog_out(MESSAGE3) << "Vertex fit failed." << endl;
    }
    auto vertexpos = Vertex.GetPos();
    auto vertexcov = Vertex.GetCov();
    auto vertexglobalpos = Vertex.GetGlobalPos();
    auto vertexglobalcov = Vertex.GetGlobalCov();
    auto vertexres = Vertex.GetRes();

    _root_vertex_u = vertexpos[0];
    _root_vertex_v = vertexpos[1];
    _root_vertex_w = vertexpos[2];
    _root_vertex_x = vertexglobalpos[0];
    _root_vertex_y = vertexglobalpos[1];
    _root_vertex_z = vertexglobalpos[2];

    _root_vertex_u_var = vertexcov(0, 0);
    _root_vertex_v_var = vertexcov(1, 1);
    _root_vertex_w_var = vertexcov(2, 2);
    _root_vertex_x_var = vertexglobalcov(0, 0);
    _root_vertex_y_var = vertexglobalcov(1, 1);
    _root_vertex_z_var = vertexglobalcov(2, 2);

    _root_vertex_chi2ndf = Vertex.GetChi2Ndof();
    _root_vertex_prob = TMath::Prob(Vertex.GetChi2(), Vertex.GetNdf());
    _root_vertex_u_res = vertexres[2];
    _root_vertex_v_res = vertexres[3];

    for (size_t idown = 0; idown < up2down[iup].size(); idown++) {
      TBTrack &downtrack = downTrackStore[up2down[iup][idown]];

      // comboChi2 is chi2 combination of track in upstream and downstream
      // telescope arm
      double comboChi2 = uptrack.GetChiSqu() + downtrack.GetChiSqu();

      // In and OutStates of the reconstructed Track at the current detector
      TBTrackState &InState = uptrack.GetTE(_idut).GetState();
      TBTrackState &OutState = downtrack.GetTE(_idut).GetState();

      // MSC Analysis for the reconstructed angles
      // Here we use the In and Out State and the GetScatterKinks function of
      // the TBKalmanMSC Class

      // Angles and angle errors
      auto theta = TrackFitterMSC.GetScatterKinks(InState, OutState);
      auto Cov = TrackFitterMSC.GetScatterKinkCov(InState, OutState);

      // Get the track parameters of the fitted track on the current sensor
      // The u and v positions are needed for a position-resolved measurement
      auto p_in = InState.GetPars();
      auto p_out = OutState.GetPars();

      // Get the covariance entries of the intersection coordinates
      auto instate_covs = InState.GetCov();
      auto outstate_covs = OutState.GetCov();

      _root_u_var = 1.0 / (1.0 / instate_covs(2, 2) +
                           1.0 / outstate_covs(2, 2)); // weighted mean
      _root_v_var = 1.0 / (1.0 / instate_covs(3, 3) +
                           1.0 / outstate_covs(3, 3)); // weighted mean

      // Fill root variables
      _root_vertex_multiplicity = up2down[iup].size();
      _root_vertex_id = iup;
      _root_momentum = uptrack.GetMomentum();
      _rootTrackProbUp = TMath::Prob(uptrack.GetChiSqu(), uptrack.GetNDF());
      _rootTrackProbDown =
          TMath::Prob(downtrack.GetChiSqu(), downtrack.GetNDF());
      _rootTrackProbCombo =
          TMath::Prob(comboChi2, uptrack.GetNDF() + downtrack.GetNDF());

      _root_u_in = p_in[2];
      _root_v_in = p_in[3];
      _root_u_out = p_out[2];
      _root_v_out = p_out[3];
      _root_u =
          (instate_covs(2, 2) * p_in[2] + outstate_covs(2, 2) * p_out[2]) /
          (instate_covs(2, 2) + outstate_covs(2, 2)); // weighted mean
      _root_v =
          (instate_covs(3, 3) * p_in[3] + outstate_covs(3, 3) * p_out[3]) /
          (instate_covs(3, 3) + outstate_covs(3, 3)); // weightes mean

      _root_angle1 = theta[0];
      _root_angle2 = theta[1];
      _root_angle1_var = Cov(0, 0);
      _root_angle2_var = Cov(1, 1);

      _root_u_in = p_in[2];
      _root_v_in = p_in[3];
      _root_u_out = p_out[2];
      _root_v_out = p_out[3];
      _root_u =
          (instate_covs(2, 2) * p_in[2] + outstate_covs(2, 2) * p_out[2]) /
          (instate_covs(2, 2) + outstate_covs(2, 2)); // weighted mean
      _root_v =
          (instate_covs(3, 3) * p_in[3] + outstate_covs(3, 3) * p_out[3]) /
          (instate_covs(3, 3) + outstate_covs(3, 3)); // weightes mean

      _root_angle1_var = Cov(0, 0);
      _root_angle2_var = Cov(1, 1);

      // Construct the u and v residuals and calculate a chi2 value from them
      auto res = p_in - p_out;
      auto res_covs = instate_covs + outstate_covs;

      // bool invertible;
      // Use only the sub matrices, which describe the spatial coordinates of
      // the trackstate
      auto jchisq = res.block<2, 1>(2, 0).transpose() *
                    res_covs.block<2, 2>(3, 3).inverse() *
                    res.block<2, 1>(2, 0);

      streamlog_out(MESSAGE1) << "Complete Covariance matrix: " << res_covs
                              << endl;
      streamlog_out(MESSAGE1) << "Part of Covariance matrix used here: "
                              << res_covs.block<2, 2>(3, 3) << endl;

      _root_chi2 = jchisq[0];
      _root_prob = TMath::Prob(jchisq[0], 2);

      if (_m_toy) {
        double dudw = p_in[0];
        double dvdw = p_in[1];
        double u = p_in[2];
        double v = p_in[3];
        // double mom = uptrack.GetMomentum();
        double l0 =
            dut.GetThickness(u, v) * std::sqrt(1 + dudw * dudw + dvdw * dvdw);

        // Simulate energy loss by bremsstrahlung (Bethe Heitler theory)
        weighted_mean_mom = (1.0-materialeffect::Epsilon_Weightfactor) * uptrack.GetCharge() / p_in[4];

   		streamlog_out ( MESSAGE1 ) << endl;
   		streamlog_out ( MESSAGE1 ) << "Materialeffects of target plane " << endl;
   		streamlog_out ( MESSAGE1 ) << "Momentum before transition: " << uptrack.GetCharge() / p_in[4] << " GeV" << endl;
   		streamlog_out ( MESSAGE1 ) << "First term of weighted mean: " << weighted_mean_mom << " GeV" << endl;

        if (_m_ToyBetheHeitler) {
		  streamlog_out ( MESSAGE3 ) << "Inside! " << endl;
          double t = l0 / dut.GetRadLength(u, v);
          double rndm = gRandom->Rndm(1);
          materialeffect::SimulateBetherHeitlerEnergyLoss(
              p_in, t, uptrack.GetMass(), uptrack.GetCharge(), rndm);
        }
        // Take average of momentum before and after scattering
        weighted_mean_mom += materialeffect::Epsilon_Weightfactor * uptrack.GetCharge() / p_in[4];

   		streamlog_out ( MESSAGE1 ) << "Momentum after transition: " << uptrack.GetCharge() / p_in[4] << " GeV" << endl;
   		streamlog_out ( MESSAGE1 ) << "Second term of weighted mean: " << materialeffect::Epsilon_Weightfactor * uptrack.GetCharge() / p_in[4] << " GeV" << endl;
   		streamlog_out ( MESSAGE1 ) << "Overall weighted mean: " << weighted_mean_mom << " GeV" << endl;

        double kink_u = 0;
        double kink_v = 0;

        if (_m_toyScatterModel == 0) {
          // Highland model scattering
          double theta2 = materialeffect::GetScatterTheta2(
              weighted_mean_mom, l0, dut.GetRadLength(u, v), uptrack.GetMass(),
              uptrack.GetCharge());
          kink_u = gRandom->Gaus(0, TMath::Sqrt(theta2));
          kink_v = gRandom->Gaus(0, TMath::Sqrt(theta2));
        } else if (_m_toyScatterModel == 1) {
          kink_u = materialeffect::GetScatterKink_SC(
              l0, dut.GetRadLength(u, v), dut.GetAtomicNumber(u, v),
              dut.GetAtomicMass(u, v), uptrack.GetMass(), uptrack.GetCharge(),
              weighted_mean_mom);
          kink_v = materialeffect::GetScatterKink_SC(
              l0, dut.GetRadLength(u, v), dut.GetAtomicNumber(u, v),
              dut.GetAtomicMass(u, v), uptrack.GetMass(), uptrack.GetCharge(),
              weighted_mean_mom);
        }

        // Scatter track ('in' state -> 'out' state)
        auto toystate = p_in;
        materialeffect::ScatterTrack(toystate, kink_u, kink_v);

        TBTrackState outstate_toy;
        outstate_toy.Pars = toystate;

        // Calculate scattering angles from
        theta = TrackFitterMSC.GetScatterKinks(InState, outstate_toy);

        double reco_error1;
        double reco_error2;

        //  In case of a selected angle reco error smaller than 0,
        //  use the real reco error
        if (_m_reco_error < 0) {
          reco_error1 = TMath::Sqrt(Cov(0, 0));
          reco_error2 = TMath::Sqrt(Cov(1, 1));
        }

        // Else use the selected reco error
        else {
          reco_error1 = _m_reco_error;
          reco_error2 = _m_reco_error;
        }

        _root_angle1 = theta[0] + gRandom->Gaus(0, reco_error1);
        _root_angle2 = theta[1] + gRandom->Gaus(0, reco_error2);

        _root_angle1_var = reco_error1 * reco_error1;
        _root_angle2_var = reco_error2 * reco_error2;

        _root_momentum = weighted_mean_mom;
      }

      _rootMscTree->Fill();

    } // end down track loop

  } // end up track loop

  return;
}

//
// Method called after each event to check the data processed
//
void X0ImageProducer::check(LCEvent *) {}

//
// Method called after all data processing
//
void X0ImageProducer::end() {

  // Print the summer
  streamlog_out(MESSAGE3) << endl
                          << endl
                          << "Total number of processed events:     "
                          << setiosflags(ios::right) << _nEvt
                          << resetiosflags(ios::right) << endl
                          << endl
                          << endl;

  // CPU time end
  _timeCPU = clock() / 1000 - _timeCPU;

  // Print message
  streamlog_out(MESSAGE3) << std::endl
                          << " "
                          << "Time per event: "
                          << std::setiosflags(std::ios::fixed |
                                              std::ios::internal)
                          << std::setprecision(3) << _timeCPU / _nEvt << " ms"
                          << std::endl
                          << std::setprecision(3) << std::endl
                          << " "
                          << "Processor succesfully finished!" << std::endl;

  // Close root  files
  _rootFile->Write();
}

//
// Method printing processor parameters
//
void X0ImageProducer::printProcessorParams() const {

  streamlog_out(MESSAGE3)
      << std::endl
      << " "
      << "X0ImageProducer Development Version, be carefull!!"
      << " " << std::endl
      << std::endl;
}

// ROOT_OUTPUT
void X0ImageProducer::bookHistos() {

  _rootFile = new TFile(_rootFileName.c_str(), "recreate");

  //
  // Event Summary Tree
  _rootEventTree = new TTree("Event", "Event info");
  _rootEventTree->Branch("iRun", &_rootRunNumber, "iRun/I");
  _rootEventTree->Branch("iEvt", &_rootEventNumber, "iEvt/I");
  _rootEventTree->Branch("nDownTracks", &_rootnDownTracks, "nDownTracks/I");
  _rootEventTree->Branch("nUpTracks", &_rootnUpTracks, "nUpTracks/I");
  _rootEventTree->Branch("nMatched", &_rootNMatched, "nMatched/I");

  //
  // Scattering angle Tree
  _rootMscTree = new TTree("MSCTree", "Multiple scattering data");
  _rootMscTree->Branch("iRun", &_rootRunNumber, "iRun/I");
  _rootMscTree->Branch("iEvt", &_rootEventNumber, "iEvt/I");
  _rootMscTree->Branch("prob_up", &_rootTrackProbUp, "prob_up/D");
  _rootMscTree->Branch("prob_down", &_rootTrackProbDown, "prob_down/D");
  _rootMscTree->Branch("prob_combo", &_rootTrackProbCombo, "prob_combo/D");

  _rootMscTree->Branch("u_in", &_root_u_in, "u_in/D");
  _rootMscTree->Branch("v_in", &_root_v_in, "v_in/D");
  _rootMscTree->Branch("u_out", &_root_u_out, "u_out/D");
  _rootMscTree->Branch("v_out", &_root_v_out, "v_out/D");
  _rootMscTree->Branch("u", &_root_u, "u/D");
  _rootMscTree->Branch("v", &_root_v, "v/D");
  _rootMscTree->Branch("u_var", &_root_u_var, "u_var/D");
  _rootMscTree->Branch("v_var", &_root_v_var, "v_var/D");

  _rootMscTree->Branch("theta1", &_root_angle1, "theta1/D");
  _rootMscTree->Branch("theta2", &_root_angle2, "theta2/D");
  _rootMscTree->Branch("theta1_var", &_root_angle1_var, "theta1_var/D");
  _rootMscTree->Branch("theta2_var", &_root_angle2_var, "theta2_var/D");
  _rootMscTree->Branch("momentum", &_root_momentum, "momentum/D");

  _rootMscTree->Branch("chi2", &_root_chi2, "chi2/D");
  _rootMscTree->Branch("prob", &_root_prob, "prob/D");

  //
  // Vertexing Tree
  _rootMscTree->Branch("vertex_u", &_root_vertex_u, "vertex_u/D");
  _rootMscTree->Branch("vertex_v", &_root_vertex_v, "vertex_v/D");
  _rootMscTree->Branch("vertex_w", &_root_vertex_w, "vertex_w/D");
  _rootMscTree->Branch("vertex_u_var", &_root_vertex_u_var, "vertex_u_var/D");
  _rootMscTree->Branch("vertex_v_var", &_root_vertex_v_var, "vertex_v_var/D");
  _rootMscTree->Branch("vertex_w_var", &_root_vertex_w_var, "vertex_w_var/D");

  _rootMscTree->Branch("vertex_x", &_root_vertex_x, "vertex_x/D");
  _rootMscTree->Branch("vertex_y", &_root_vertex_y, "vertex_y/D");
  _rootMscTree->Branch("vertex_z", &_root_vertex_z, "vertex_z/D");
  _rootMscTree->Branch("vertex_x_var", &_root_vertex_x_var, "vertex_x_var/D");
  _rootMscTree->Branch("vertex_y_var", &_root_vertex_y_var, "vertex_y_var/D");
  _rootMscTree->Branch("vertex_z_var", &_root_vertex_z_var, "vertex_z_var/D");

  _rootMscTree->Branch("vertex_chi2ndf", &_root_vertex_chi2ndf,
                       "vertex_chi2ndf/D");
  _rootMscTree->Branch("vertex_prob", &_root_vertex_prob, "vertex_prob/D");
  _rootMscTree->Branch("vertex_u_res", &_root_vertex_u_res, "vertex_u_res/D");
  _rootMscTree->Branch("vertex_v_res", &_root_vertex_v_res, "vertex_v_res/D");

  _rootMscTree->Branch("vertex_multiplicity", &_root_vertex_multiplicity,
                       "_root_vertex_multiplicity/I");
  _rootMscTree->Branch("vertex_id", &_root_vertex_id, "_root_vertex_id/I");
}

} // Namespace
