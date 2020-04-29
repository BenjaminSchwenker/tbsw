// GoeClusterCalibrator Processor  
// 
// See GoeClusterCalibrator.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "GoeClusterCalibrator.h"

// TBTools includes
#include "TBDetector.h"
#include "Utilities.h"
#include "TBTrack.h"
#include "TrackInputProvider.h"
#include "GenericTrackFitter.h"
#include "PixelCluster.h"
#include "PolyClusterDescriptor.h"

// Include basic C
#include <iostream>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <sstream>

// Include LCIO classes
#include "lcio.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <Exceptions.h>

// Include ROOT classes
#include <TMath.h>
#include <TFile.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>


// Used namespaces
using namespace std; 
using namespace lcio ;
using namespace marlin ;

namespace depfet {

  //
  // Instantiate this object
  //
  GoeClusterCalibrator aGoeClusterCalibrator ;
  
  //
  // Constructor
  //
  GoeClusterCalibrator::GoeClusterCalibrator() : Processor("GoeClusterCalibrator")
  {
    
    // Processor description
    _description = "GoeClusterCalibrator: Create clusterDB for clusters using reco tracks" ;
    
    //
    // Input collections 
    
    registerInputCollection(LCIO::TRACK,"InputTrackCollectionName",
                            "Track input collection",
                            _inputTrackCollectionName,std::string("tracks"));
    
    registerProcessorParameter( "ClusterDBFileName",
                                "Output clusterDB file name",
                                _clusterDBFileName, std::string("clusterDB.root"));  
     
    registerProcessorParameter ("MinClusters",
                                "Minimum number of clusters needed to calibrate a cluster type",
                                _minClusters,  static_cast < int > (200));
    
    registerProcessorParameter ("MaxEtaBins",
                                "Maximum number of eta bins per cluster type",
                                _maxEtaBins,  static_cast < int > (1));
    
    registerProcessorParameter ("vCellPeriod",
                                "Periodicity for vCells used for clusterDB",
                                _vCellPeriod,  static_cast < int > (1));
     
    registerProcessorParameter ("uCellPeriod",
                                "Periodicity for uCells used for clusterDB",
                                _uCellPeriod,  static_cast < int > (1));
    
    registerProcessorParameter ("MinSigmaU", "Minimum value of sigmaU for 2x2 cluster covariance matrix [mm]",
                                _minSigmaU, static_cast < float > (0.5E-3) );
    
    registerProcessorParameter ("MinSigmaV", "Minimum value of sigmaV for 2x2 cluster covariance matrix [mm]",
                                _minSigmaV, static_cast < float > (0.5E-3) );
     
    
    std::vector<int> initSelectIDVec;
    registerProcessorParameter ("SelectPlanes",
                                "Select clusters from list of planes",
                                _selectIDVec, initSelectIDVec);
    
  }
  
  //
  // Method called at the beginning of data processing
  //
  void GoeClusterCalibrator::init() {
    
    // Initialize variables
    _nRun = 0 ;
    _nEvt = 0 ;
    _timeCPU = clock()/1000 ;
    
    if ( _minSigmaU <= 0 ) {  
      _minSigmaU = 1.0E-3;  
    } 
    
    if ( _minSigmaV <= 0 ) {
      _minSigmaV = 1.0E-3;  
    } 
    
    if ( _maxEtaBins<=0) _maxEtaBins=1;

    // Create a useful name for a file that should be deleted after Marlin is finished
    // There should be one collector output for each clusterDB
    // Add prefix "tmp" to indicate that the file is temporary 
    
    std::set<char> delims{'/'};
    _collectorOutputFileName = "tmpCollectorOutputFor_" + splitpath(_clusterDBFileName, delims).back();
    
    _rootCollectorOutputFile = new TFile( _collectorOutputFileName.c_str(),"recreate");
    _rootCollectorOutputFile->cd("");
    
    _trackVarUHisto = new TH1F("trackVarUHisto","trackVarUHisto",1,0,1);
    _trackVarUHisto->StatOverflows(); 	            
    
    _trackVarVHisto = new TH1F("trackVarVHisto","trackVarVHisto",1,0,1); 
    _trackVarVHisto->StatOverflows();         
    
    _trackCovUVHisto = new TH1F("trackCovUVHisto","trackCovUVHisto",1,0,1);
    _trackCovUVHisto->StatOverflows(); 	

    _trackDuDwHisto = new TH1F("trackDuDwHisto","trackDuDwHisto",1,0,1);
    _trackDuDwHisto->StatOverflows(); 
     
    _trackDvDwHisto = new TH1F("trackDvDwHisto","trackDvDwHisto",1,0,1);
    _trackDuDwHisto->StatOverflows();      
    
    m_rootTree = new TTree("tree","Cluster info");
    m_rootTree->Branch<string>("TypeName", &m_typeName);
    m_rootTree->Branch< std::vector<float> >("ClusterEtaSS", &m_clusterEtaSS);
    m_rootTree->Branch< std::vector<float> >("ClusterEtaOS", &m_clusterEtaOS);
    m_rootTree->Branch<float>("OffsetU", &m_positionOffsetU);
    m_rootTree->Branch<float>("OffsetV", &m_positionOffsetV);
    m_rootTree->Branch<float>("CogResidualU", &m_cogResidualU);
    m_rootTree->Branch<float>("CogResidualV", &m_cogResidualV);
    m_rootTree->Branch<int>("RunNumber", &m_runNumber);
    m_rootTree->Branch<int>("PlaneNumber", &m_planeNumber);
    
    // Print set parameters
    printProcessorParams();
    
   
  }
  
  //
  // Method called for each run
  //
  void GoeClusterCalibrator::processRunHeader(LCRunHeader * run)
  {
   
    // Print run number
    streamlog_out(MESSAGE3) << "Processing run: "
                            << (run->getRunNumber())
                            << std::endl << std::endl;
    
    _nRun++ ;

  }
  
  //
  // Method called for each event
  //
  void GoeClusterCalibrator::processEvent(LCEvent * evt)
  {
    
    // Set run number 
    m_runNumber = evt->getRunNumber();   
    _nEvt ++ ;
    
    TrackInputProvider TrackIO; 
    
    GenericTrackFitter TrackFitter(TBDetector::GetInstance());
    TrackFitter.SetNumIterations(1); 
    
    LCCollection* inputCollection;
    try {
      inputCollection = evt->getCollection(_inputTrackCollectionName);
    } catch (DataNotAvailableException& e) {
      throw SkipEventException(this);
    }
     
    // Main loop over all tracks
    int nTracks = inputCollection->getNumberOfElements(); 
    for (int itrk = 0; itrk < nTracks; itrk++) {   
      
      // Retrieve track from LCIO 
      Track * inputtrack = dynamic_cast<Track*> (inputCollection->getElementAt(itrk));
      
      // Convert LCIO -> TB track  
      TBTrack track = TrackIO.MakeTBTrack( inputtrack, TBDetector::GetInstance() );  
      
      // Refit track 
      bool trkerr = TrackFitter.Fit(track);
      if ( trkerr ) {
        continue;
      } 
      
      //
      // Loop over all clusters in the track
      for (int ipl= 0; ipl< TBDetector::GetInstance().GetNSensors(); ++ipl) {   
        
        
        // Get sensor data 
        //------------------------
        TBTrackElement& TE = track.GetTE(ipl);  
        const Det & Sensor = TBDetector::Get(ipl);  
        
        bool selectID = false;
        for (auto id :  _selectIDVec)  {
          if  (id == ipl) selectID = true; 
        }
         
        // Ignore track elements w/o measurment
        if ( TE.HasHit() && selectID ) { 
          
          // This is the plane number of one plane to which 
          // the clusterDB would be applied
          _setOfPlaneNumbers.insert(ipl);
          
          // Get local track parameters 
          double trk_tu = TE.GetState().GetPars()[0];  // rad
          double trk_tv = TE.GetState().GetPars()[1];  // rad
          double trk_u = TE.GetState().GetPars()[2];   // mm
          double trk_v = TE.GetState().GetPars()[3];   // mm
          
           
          double sigma2_u = TE.GetState().GetCov()(2,2); 
          double sigma2_v = TE.GetState().GetCov()(3,3); 
          double cov_uv   = TE.GetState().GetCov()(2,3);  
           
          // We use temporary histograms to compute an averaged 2x2 
          // covariance matrix for all reco tracks at given sensor.
             
          _trackVarUHisto->Fill(sigma2_u);
          _trackVarVHisto->Fill(sigma2_v); 
          _trackCovUVHisto->Fill(cov_uv);  
          _trackDuDwHisto->Fill(trk_tu);  
          _trackDvDwHisto->Fill(trk_tv);       
             
          PixelCluster Cluster = TE.GetHit().GetCluster();  
          PolyClusterDescriptor Descriptor(Cluster, Sensor);
          
          // Compute center of gravity hit coordinate
          double u_cog{0.0}, v_cog{0.0}, sig2_u_cog{0.0}, sig2_v_cog{0.0}, cov_uv_cog{0.0};
          Cluster.getCenterOfGravity(Sensor, u_cog, v_cog, sig2_u_cog, sig2_v_cog, cov_uv_cog); 
          
          // Fill collector tree with all data needed for this calibration
          
          // A string to identify the cluster type, it quantifies the configuration of firing pixels 
          // but does not use the measured pixel signals. Details depend on the implementation of 
          // the cluster descriptor. 
          m_typeName = Descriptor.getType(_vCellPeriod, _uCellPeriod);
          
          // The eta=s_tail/(s_head+s_tail) is computed from the charge deposited where the particle 
          // enters (tail charge) and leaves (head charge) the sensor. There are differnt ways to 
          // compute head/tail charges depending on cluster type and sign of incidence angles. Here, 
          // we compute a number of relevant choices and leave the best selection to the final 
          // calibration. 
          m_clusterEtaSS.clear();
          m_clusterEtaOS.clear();
          for (int etaIndex=0; etaIndex<9; etaIndex++) {
            m_clusterEtaSS.push_back( Descriptor.computeEta(+1, +1, etaIndex) );
            m_clusterEtaOS.push_back( Descriptor.computeEta(+1, -1, etaIndex) ); 
          }
          
          m_positionOffsetU = trk_u - Descriptor.getOriginU();  
          m_positionOffsetV = trk_v - Descriptor.getOriginV();  
          m_cogResidualU = trk_u - u_cog;
          m_cogResidualV = trk_v - v_cog;
          m_planeNumber = ipl;
          m_rootTree->Fill(); 
           
          // Create the residual histograms for bias correction if needed  
          auto key = std::make_pair(m_runNumber, m_planeNumber);
          if ( _biasMapU.find(key) == _biasMapU.end() ) {
            string tmpName = "hBiasU_runno_" + to_string(m_runNumber) + "_plane_" + to_string(m_planeNumber); 
            _biasMapU[key] = new TH1F(tmpName.c_str(),tmpName.c_str(),1,0,1);
            _biasMapU[key]->StatOverflows();
            _biasMapU[key]->SetDirectory(nullptr);
          }
          if ( _biasMapV.find(key) == _biasMapV.end() ) {
            string tmpName = "hBiasV_runno_" + to_string(m_runNumber) + "_plane_" + to_string(m_planeNumber); 
            _biasMapV[key] = new TH1F(tmpName.c_str(),tmpName.c_str(),1,0,1);
            _biasMapV[key]->StatOverflows();
            _biasMapV[key]->SetDirectory(nullptr);
          }
          
          // Fill the bias map to detect residual drift over runs
          _biasMapU[key]->Fill(TE.GetHit().GetCoord()[0] - trk_u);   
          _biasMapV[key]->Fill(TE.GetHit().GetCoord()[1] - trk_v);   
        }
      }
    }  
    
    return;
  }
  
  
  //
  // Method called after each event to check the data processed
  //
  void GoeClusterCalibrator::check( LCEvent * ) {}
  
  //
  // Method called after all data processing
  //
  void GoeClusterCalibrator::end()
  {
    
    // CPU time end
    _timeCPU = clock()/1000 - _timeCPU;
    
    // Print message
    streamlog_out(MESSAGE3) << std::endl
                            << " "
                            << "Time per event: "
                            << std::setiosflags(std::ios::fixed | std::ios::internal )
                            << std::setprecision(3)
                            << _timeCPU/_nEvt
                            << " ms"
                            << std::endl
                            << std::setprecision(3)
                            << std::endl
                            << " "
                            << "Processor succesfully finished!"
                            << std::endl;
    
    // Make sure that all data was collected for a set of planes 
    // have the same protopixels. 
    bool consistencyTestPassed = true;
    std::vector<int> planeNumbersVec(_setOfPlaneNumbers.begin(), _setOfPlaneNumbers.end()); 
    for(auto i : planeNumbersVec ) {
      for(auto j : planeNumbersVec ) {    
        if ( TBDetector::Get(i).GetProtopixels() ==  TBDetector::Get(j).GetProtopixels() ) {
          streamlog_out(MESSAGE3) << "Plane " << i << " and plane " << j << " have identical map of protopixels." << std::endl;   
        } else {
          streamlog_out(MESSAGE3) << "Plane " << i << " and plane " << j << " have different map of protopixels!" << std::endl;   
          consistencyTestPassed = false;
        } 
      }  
    }  
    
    // Compute the average 2x2 covariance matrix for the 
    // estimated track states. We assume the track state
    // to be unbiased and will later subtract this number
    // from the residuals.  
    double trk_covU = _trackVarUHisto->GetMean();     
    double trk_covV = _trackVarVHisto->GetMean();  
    double trk_covUV = _trackCovUVHisto->GetMean();   
  
    // Compute the average incidence angles into the sensor
    // We assume a strongly collimated beam and the rms of the 
    // should be small (< few mrad) 
    double thetaU = _trackDuDwHisto->GetMean();
    double thetaV = _trackDvDwHisto->GetMean();

    // Check if there is a drift in the residuals in case 
    // data from more than one run is used. 
    for (const auto &entry: _biasMapU)
	{
      auto key_pair = entry.first;
      streamlog_out(MESSAGE3) << std::setprecision(5)
                              << "Run number " << key_pair.first << ", plane number " << key_pair.second << " has biasU="
                              << entry.second->GetMean() << "+/-" << entry.second->GetMeanError() << endl;
	}
    // Check if there is a drift in the residuals in case 
    // data from more than one run is used. 
    for (const auto &entry: _biasMapV)
	{
      auto key_pair = entry.first;
      streamlog_out(MESSAGE3) << std::setprecision(5)
                              << "Run number " << key_pair.first << ", plane number " << key_pair.second << " has biasV="
                              << entry.second->GetMean() << "+/-" << entry.second->GetMeanError() << endl;
	}
     
    streamlog_out(MESSAGE3) << std::setprecision(5)
                            << "Average telescope sigmaU=" << sqrt(trk_covU) << ", sigmaV=" << sqrt(trk_covV) << endl
                            << "Average track incidence angles on DUT are du/dw=" << thetaU  << ", dvdw=" << thetaV  << endl;
     
    // Number of tracks in the data sample
    const auto nEntries = m_rootTree->GetEntries();
    
    // Enumerate all types by unique name and count their
    // occurence in training data.
    vector< pair<string, float> > allTypes;
    
    for (int i = 0; i < nEntries; ++i) {
      m_rootTree->GetEntry(i);
      
      auto it = std::find_if(allTypes.begin(), allTypes.end(),
                             [&](const pair<string, float>& element) { return element.first == m_typeName;});
      
      //Shape name exists in vector
      if (it != allTypes.end()) {
        //increment key in map
        it->second++;
      }
      //Shape name does not exist
      else {
        //Not found, insert in vector
        allTypes.push_back(pair<string, int>(m_typeName, 1));
      }
    }
    
    // Make a list of good types having enough data for 
    // calibration. 
    
    // Vector with names of types having enough data
    vector< string > goodTypes;
    
    // Ratio of tracks in goodTypes relative to all tracks
    double coverage = 0.0;
    
    for (auto iter : allTypes) {
      auto name = iter.first;
      auto counter = iter.second;
      if (counter >=  _minClusters) {
        coverage += counter / nEntries;    
        goodTypes.push_back(name);
      } else {
        streamlog_out(MESSAGE3) << "  Unable to calibrate cluster type:  " << name << " because too few counts (" << counter  << ")" << endl;
      } 
    }     
    
    // Make 2D histogram for center of gravity residuals 
    // for all good types. 
    vector< pair< string, TH2D > > cogHistosVec;
    
    for (const string& name : goodTypes) {
      string cog_name = "cog_" + name;    
      TH2D cogHisto(cog_name.c_str(), cog_name.c_str(), 1, 0, 1, 1, 0, 1);
      cogHisto.SetDirectory(nullptr);
      cogHisto.StatOverflows();
      cogHistosVec.push_back(pair<string, TH2D>(name, cogHisto));  
    }     
    
    // Loop over the tree is to fill cog histograms
    for (int i = 0; i < nEntries; ++i) {
      m_rootTree->GetEntry(i); 
      
      auto it = std::find_if(cogHistosVec.begin(), cogHistosVec.end(),
                        [&](const pair<string, TH2D>& element) { return element.first == m_typeName;});
      
      //Item exists in maps
      if (it != cogHistosVec.end()  ) {
        it->second.Fill(m_cogResidualU, m_cogResidualV);   
      }
    }
    
    // Great, at this point we know the spatial resolution for the 
    // center of gravity method on all goodTypes. This is our baseline
    // to judge the resolution of eta shapes later. 
    
    // Map key is cluster type and value is center of gravity loss 
    map< string, double > cogLossMap;
    
    streamlog_out(MESSAGE3) << "Summary of center of gravity performance per type: "  << endl <<  endl;
    
    for (const auto &entry: cogHistosVec)
	{
      double covU = pow(entry.second.GetRMS(1), 2);
      double covV = pow(entry.second.GetRMS(2), 2);
      double covUV = entry.second.GetCovariance();
      
      covU =  std::max(covU - trk_covU, pow(_minSigmaU ,2));  
      covV =  std::max(covV - trk_covV, pow(_minSigmaV ,2));  
      covUV -= trk_covUV;    
      
      streamlog_out(MESSAGE3) << std::setprecision(5)
                              << "Name " << entry.first  << " entries=" << entry.second.GetEntries() 
                              << ", sigmaU=" << sqrt(covU) << ", sigmaV=" << sqrt(covV) << ", corrUV=" << covUV/sqrt(covU)/sqrt(covV) 
                              << ", loss=" << sqrt(covU)+sqrt(covV)
                              << std::endl;  
      
      TMatrixDSym HitCov(2);
      HitCov(0, 0) = covU;
      HitCov(1, 0) = covUV;
      HitCov(0, 1) = covUV;
      HitCov(1, 1) = covV;
          
      TMatrixDSymEigen HitCovE(HitCov);
      TVectorD eigenval = HitCovE.GetEigenValues();
      if (eigenval(0) <= 0 || eigenval(1) <= 0) {
        streamlog_out(ERROR3) << "Estimated covariance matrix not positive definite." << std::endl;
      }
      
      // We want to use this to judge the performance of eta calibrations. 
      cogLossMap[entry.first] = sqrt(covU)+sqrt(covV); 
	}
    
    // There are multiple ways to compute an eta variable for a given 
    // cluster type. These ways differ in the definition of the head 
    // and tail signal entering for computing eta=S_h/(S_t + S_h).   
    // The training data has computed a number of candidate eta 
    // definitions eta[i].
    //   
    // In the next step, we compute for each variable eta[i] the cluster 
    // shape corrections. We measure the performance of the corrections for 
    // on a cluster type T with a loss function L(T,eta[i]). The final 
    // set of cluster shape corrections is computed using the best eta 
    // variable eta[i_T] for each cluster type T.   
    
    // Map key is cluster type and value is vector of losses for different eta[i]
    map< string, vector<double> > lossMatrix;
    
    // Map key is cluster type and value is vector of eta binnings for different eta[i]
    map< string, vector< vector<double> > > etaBinningMap;
    
    // Map key is cluster type and value is vector of eta offset calibrations for different eta[i]
    map< string, vector< vector< vector<double> > > > etaOffsetCalibrationMap;
    
    // Loop over all collected eta variables and compute the losses
    for (int etaIndex = 0; etaIndex < 9; etaIndex++) {        
      
      // Vector with eta histograms 
      vector< pair<string, TH1D> > etaHistos;
      
      for (const string& name : goodTypes) {
        string eta_name =  "variant_" + to_string(etaIndex) + "_" + "eta_"  + name;     
        TH1D etaHisto(eta_name.c_str(), eta_name.c_str(), 301, 0, 1);
        etaHisto.SetDirectory(nullptr);
        etaHistos.push_back(pair<string, TH1D>(name, etaHisto));
      }     
      
      // Loop over the tree is to fill the eta histograms for
      // good types
      for (int i = 0; i < nEntries; ++i) {
        m_rootTree->GetEntry(i);
        auto it = std::find_if(etaHistos.begin(), etaHistos.end(),
                      [&](const pair<string, TH1D>& element) { return element.first == m_typeName;});
        //Item exists in map
        if (it != etaHistos.end()) {
          // increment key in map
          auto clusterEta = m_clusterEtaSS[etaIndex];
          if (thetaU * thetaV < 0) {clusterEta = m_clusterEtaOS[etaIndex];}
          it->second.Fill(clusterEta);
        }
      }
      
      // Vector for eta bin edges stored by type name 
      vector< pair< string, vector<double> > > etaBinEdgesVec;
      
      // Vector for offset histograms storing pairs of typename and a vector of 2d histos for different eta bins
      vector< pair< string, vector<TH2D> > > offsetHistosVec;
      
      for (auto iter : etaHistos) {
        auto name = iter.first;
        auto& histo = iter.second;
        int nClusters = histo.GetEntries();
        
        streamlog_out(MESSAGE3) << "Eta histogram " << histo.GetName() << " has " << nClusters << " entries" << std::endl; 
        
        // Try to split clusters into n bins with _minClusters clusters
        int nEtaBins  = std::max(int(nClusters / _minClusters), 1);
        nEtaBins =  std::min(nEtaBins, _maxEtaBins);  
        
        // We have to check for singular cases where eta distribution is concentrated is concentrated in 
        // less bins than required number of eta bins 
        int nFilledBins = 0;  
        for (auto ibin = 1; ibin <= histo.GetXaxis()->GetNbins(); ibin++) {  
          if (histo.GetBinContent(ibin) > 0) nFilledBins++;  
        }
        if (nFilledBins <= nEtaBins) {
          nEtaBins = 1; 
          streamlog_out(MESSAGE3) << "Eta histogram " << histo.GetName() << " is a delta spike" << std::endl;   
        }
        
        vector<double> etaBinEdges;
        vector< TH2D > offsetHistos;
        
        for (int i = 0; i < nEtaBins; i++) {
          // Position where to compute the quantiles in [0,1]
          double xq = double(i) / nEtaBins;
          // Double to contain the quantile
          double yq = 0;
          histo.GetQuantiles(1, &yq, &xq);
          streamlog_out(MESSAGE3) << " Quantile at xq =" << xq << " is yq=" << yq << std::endl;
          etaBinEdges.push_back(yq);
          
          string offsetname =  "variant_" + to_string(etaIndex) + "_" + "E" + std::to_string(i) + name;         
          TH2D offsetHisto(offsetname.c_str(), offsetname.c_str(), 1, 0, 1, 1, 0, 1);
          offsetHisto.StatOverflows();
          offsetHisto.SetDirectory(nullptr);
          offsetHistos.push_back(offsetHisto);
        }
        etaBinEdgesVec.push_back(pair< string, vector<double> >(name, etaBinEdges));
        offsetHistosVec.push_back(pair< string, vector<TH2D> >(name, offsetHistos));
        
        etaBinningMap[name].push_back(vector<double>(etaBinEdges));
      }
      
      // Loop over the tree is to fill offset histograms
      for (int i = 0; i < nEntries; ++i) {
        m_rootTree->GetEntry(i); 
        
        auto it = std::find_if(offsetHistosVec.begin(), offsetHistosVec.end(),
                          [&](const pair<string, vector<TH2D>>& element) { return element.first == m_typeName;});
        
        auto it2 = std::find_if(etaBinEdgesVec.begin(), etaBinEdgesVec.end(),
                          [&](const pair<string, vector<double>>& element) { return element.first == m_typeName;});
        
        //Item exists in maps
        if (it != offsetHistosVec.end()  && it2 != etaBinEdgesVec.end() ) {
          auto clusterEta = m_clusterEtaSS[etaIndex];
          if (thetaU * thetaV < 0) {clusterEta = m_clusterEtaOS[etaIndex];}
          auto etaBin = PolyClusterDescriptor::computeEtaBin(clusterEta, it2->second);
          it->second.at(etaBin).Fill(m_positionOffsetU, m_positionOffsetV);
        }
      }
      
      // Compute the losses for all good types 
      for (auto iter : offsetHistosVec) {
        auto name = iter.first;
        auto& histovec = iter.second;
        
        double lossU = 0;
        double lossV = 0;
        double norm = 0;
                
        vector< vector<double> > offsetCalibrations;

        // Loop over eta bins
        for (auto& histo : histovec) {
            
          double offsetU = histo.GetMean(1);
          double offsetV = histo.GetMean(2);
          double covU = pow(histo.GetRMS(1), 2);
          double covV = pow(histo.GetRMS(2), 2);
          double covUV = histo.GetCovariance();
          double error_offsetU = histo.GetMeanError(1);
          double error_offsetV = histo.GetMeanError(2);
          double error_covU = 2 * histo.GetRMS(1) * histo.GetRMSError(1); 
          double error_covV = 2 * histo.GetRMS(2) * histo.GetRMSError(2); 
          
          covU =  std::max(covU - trk_covU, pow(_minSigmaU ,2));  
          covV =  std::max(covV - trk_covV, pow(_minSigmaV ,2));  
          covUV -= trk_covUV;    
          
          streamlog_out(MESSAGE3) << std::setprecision(5)
                                  << "Name " << histo.GetName()  << " entries=" << histo.GetEntries() 
                                  << ", sigmaU=" << sqrt(covU) << ", sigmaV=" << sqrt(covV) << ", corrUV=" << covUV/sqrt(covU)/sqrt(covV) << std::endl;  
           
          TMatrixDSym HitCov(2);
          HitCov(0, 0) = covU;
          HitCov(1, 0) = covUV;
          HitCov(0, 1) = covUV;
          HitCov(1, 1) = covV;
          
          TMatrixDSymEigen HitCovE(HitCov);
          TVectorD eigenval = HitCovE.GetEigenValues();
          if (eigenval(0) <= 0 || eigenval(1) <= 0) {
            streamlog_out(ERROR3) << "Estimated covariance matrix not positive definite." << std::endl;
          }

          vector<double> etaBinOffsetCalibration;
          etaBinOffsetCalibration.push_back(histo.GetEntries());
          etaBinOffsetCalibration.push_back(offsetU);
          etaBinOffsetCalibration.push_back(offsetV);
          etaBinOffsetCalibration.push_back(covU);
          etaBinOffsetCalibration.push_back(covV);
          etaBinOffsetCalibration.push_back(covUV);
          etaBinOffsetCalibration.push_back(error_offsetU);
          etaBinOffsetCalibration.push_back(error_offsetV);
          etaBinOffsetCalibration.push_back(error_covU);
          etaBinOffsetCalibration.push_back(error_covV);
          
          // Append the calibration 
          offsetCalibrations.push_back(etaBinOffsetCalibration); 
          
          // Compute weighed mean of covariance matrix diagonals
          lossU += histo.GetEntries() * covU;
          lossV += histo.GetEntries() * covV;
          norm += histo.GetEntries();
        }
         
        // Append the normalized loss, i.e. sum of average sigmaU and sigmaV 
        lossMatrix[name].push_back( sqrt(lossU/norm) + sqrt(lossV/norm) );
        
        // Append the offset calibrations 
        etaOffsetCalibrationMap[name].push_back( offsetCalibrations );
      }
    }
    
    // The matrix of losses is now computed. Next we search for the 
    // best eta index i_T for each type T. 
    
    // Map key is cluster type and value is index i_T of best eta[i_T]
    map< string, int > bestEtaMap;
    
    for (auto const& x : lossMatrix) { 
      auto smallest = std::min_element( x.second.begin(), x.second.end() );
      if (smallest == x.second.end()) {
        streamlog_out(ERROR3) << "Cannot find minimum loss for " << x.first << ". Very odd." << endl;   
      }
      bestEtaMap[x.first] = std::distance(x.second.begin(), smallest); 
         
      streamlog_out(MESSAGE3) << std::setprecision(5)
                              << "Name " << x.first  << " has cog_loss=" << cogLossMap[x.first]  << ", eta_loss=" << *smallest << " at eta index=" << bestEtaMap[x.first] <<  endl;
      streamlog_out(MESSAGE3) << std::setprecision(5)
                              << "Name " << x.first  << " has reference loss=" << x.second[0]  <<  endl;
    }
    
    // Count total number of good cluster types
    int nTypes = (int) etaOffsetCalibrationMap.size();
      
    // Count total number of eta shapes 
    int nShapes = 0; 
    for (auto const& x : etaOffsetCalibrationMap) {
      nShapes += (int)x.second[bestEtaMap[x.first]].size();
    }
    
    // Close collector file
    _rootCollectorOutputFile->Write();
    _rootCollectorOutputFile->Close();
    delete _rootCollectorOutputFile;   
    
    if (nShapes > 0 and consistencyTestPassed) {
      
      streamlog_out(MESSAGE3) << "Create the clusterDB ... " << endl; 
      
      TFile * _rootClusterDBFile = new TFile( _clusterDBFileName.c_str(),"recreate");
      _rootClusterDBFile->cd("");
      
      // Book histograms for clusterDB
      
      string histoName;  
      
      histoName = "hDB_Coverage";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",1,0,1);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("coverage [%]");
      _histoMap[histoName]->SetBinContent( 1, 100 * coverage );
      _histoMap[histoName]->GetXaxis()->SetBinLabel( 1, "cluster found in clusterDB" );
      
      histoName = "hDB_Types";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nTypes,0,nTypes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("type");  
      
      histoName = "hDB_EtaIndex";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nTypes,0,nTypes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("eta index");  
      
      histoName = "hDB_Weight";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nShapes,0,nShapes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("label weight");  
      
      histoName = "hDB_U";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nShapes,0,nShapes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("offset u [mm]");  
      
      histoName = "hDB_V";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nShapes,0,nShapes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("offset v [mm]");  
      
      histoName = "hDB_Sigma2_U";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nShapes,0,nShapes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("sigma2 offset u [mm^2]");        
      
      histoName = "hDB_Sigma2_V";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nShapes,0,nShapes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("sigma2 offset v [mm^2]"); 
      
      histoName = "hDB_Cov_UV";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nShapes,0,nShapes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("covariance u-v [mm^2]"); 
      
       
      // Compute the moments of the offset histograms 
      int offsetBin = 0;
      int typeBin = 0; 
      
      for (auto const& x : etaOffsetCalibrationMap) {
        typeBin++;
        auto name = x.first;
        auto& calvec = x.second[bestEtaMap[x.first]];
        int typeCounter = 0; 
        
        // Loop over eta bins
        int etaBin = -1; 
        for (auto& cal : calvec) {
          // Compute offset moments
          offsetBin++;
          etaBin++; 
          auto shapeName = "E" + std::to_string(etaBin) + name;
          
          // Read back the calibration constants
          int counter = (int)cal[0];
          double offsetU = cal[1];
          double offsetV = cal[2];
          double covU = cal[3];
          double covV = cal[4];
          double covUV = cal[5];
          double error_offsetU = cal[6];
          double error_offsetV = cal[7];     
          double error_covU = cal[8];
          double error_covV = cal[9];          

          typeCounter += counter;
           
          streamlog_out(MESSAGE3) << std::setprecision(5)
                                  << "Name " << shapeName  << " entries=" << counter << ", posU=" << offsetU << ", posV=" << offsetV 
                                  << ", sigmaU=" << sqrt(covU) << ", sigmaV=" << sqrt(covV) << ", corrUV=" << covUV/sqrt(covU)/sqrt(covV) << std::endl;
          
          TMatrixDSym HitCov(2);
          HitCov(0, 0) = covU;
          HitCov(1, 0) = covUV;
          HitCov(0, 1) = covUV;
          HitCov(1, 1) = covV;
          
          TMatrixDSymEigen HitCovE(HitCov);
          TVectorD eigenval = HitCovE.GetEigenValues();
          if (eigenval(0) <= 0 || eigenval(1) <= 0) {
            streamlog_out(ERROR3) << "Estimated covariance matrix not positive definite." << std::endl;
          }
          
          // Store calibration result   
          histoName = "hDB_Weight";
          _histoMap[histoName]->SetBinContent( offsetBin, counter );
          _histoMap[histoName]->GetXaxis()->SetBinLabel( offsetBin, shapeName.c_str() );
          
          histoName = "hDB_U";
          _histoMap[histoName]->SetBinContent( offsetBin, offsetU );
          _histoMap[histoName]->SetBinError( offsetBin, error_offsetU );
          _histoMap[histoName]->GetXaxis()->SetBinLabel( offsetBin, shapeName.c_str() );
          
          histoName = "hDB_V"; 
          _histoMap[histoName]->SetBinContent( offsetBin, offsetV );
          _histoMap[histoName]->SetBinError( offsetBin, error_offsetV );
          _histoMap[histoName]->GetXaxis()->SetBinLabel( offsetBin, shapeName.c_str() );
               
          histoName = "hDB_Sigma2_U";
          _histoMap[histoName]->SetBinContent( offsetBin, covU );
          _histoMap[histoName]->SetBinError( offsetBin, error_covU );
          _histoMap[histoName]->GetXaxis()->SetBinLabel( offsetBin, shapeName.c_str() );
          
          histoName = "hDB_Sigma2_V";
          _histoMap[histoName]->SetBinContent( offsetBin, covV );
          _histoMap[histoName]->SetBinError( offsetBin, error_covV );
          _histoMap[histoName]->GetXaxis()->SetBinLabel( offsetBin, shapeName.c_str() );  
        
          histoName = "hDB_Cov_UV";
          _histoMap[histoName]->SetBinContent( offsetBin, covUV );
          _histoMap[histoName]->SetBinError( offsetBin, 0 );
          _histoMap[histoName]->GetXaxis()->SetBinLabel( offsetBin, shapeName.c_str() );  
        }

        // Add cluster type 
        histoName = "hDB_Types";
        _histoMap[histoName]->SetBinContent( typeBin, typeCounter );
        _histoMap[histoName]->GetXaxis()->SetBinLabel( typeBin, name.c_str() ); 

        // Add optimized eta index 
        histoName = "hDB_EtaIndex";
        _histoMap[histoName]->SetBinContent( typeBin, bestEtaMap[x.first] );
        _histoMap[histoName]->GetXaxis()->SetBinLabel( typeBin, name.c_str() ); 
      }
      
      for (auto const& x : etaBinningMap) {
        auto name = x.first;
        auto& etaBinEdges = x.second[bestEtaMap[x.first]];
        
        // Add eta bin edges for type
        TVectorD DB_etaBinEdges( etaBinEdges.size() );
        for ( size_t iBin=0; iBin < etaBinEdges.size(); iBin++) {
          DB_etaBinEdges[iBin] = etaBinEdges[iBin];
        }
        DB_etaBinEdges.Write(  string("DB_etaBinEdges_"+name).c_str() );
      }
         
      TVectorD DB_periods( 2 );
      DB_periods[0] = _vCellPeriod;
      DB_periods[1] = _uCellPeriod;
      DB_periods.Write("DB_periods");
      
      TVectorD DB_angles( 2 );
      DB_angles[0] = thetaU;
      DB_angles[1] = thetaV;
      DB_angles.Write("DB_angles");
      
      TVectorD DB_telcov( 3 );
      DB_telcov[0] = trk_covU;
      DB_telcov[1] = trk_covV;
      DB_telcov[2] = trk_covUV;
      DB_telcov.Write("DB_telcov");
      
      if (not planeNumbersVec.empty()) {
        // Loop over all protopixel for one sensor to which the cluseterDB is going to 
        // be used. 
        for (auto protopixel : TBDetector::Get(planeNumbersVec[0]).GetProtopixels())  {
          auto pixeltype = protopixel.first;
          auto points = protopixel.second;
          
          // Add eta bin edges for type
          TVectorD DB_protopixel( 2*points.size() );
          int iBin = 0;
          for (auto point: points){
            DB_protopixel[iBin] = std::get<0>(point);
            DB_protopixel[iBin+1] = std::get<1>(point);
            iBin+=2;
          }
          DB_protopixel.Write( string("DB_protopixel_"+std::to_string(pixeltype)).c_str() ); 
        }
      } else {
        streamlog_out(ERROR3) << "No information about protopixels is available. This is strange! "  << endl;   
      }
      
      streamlog_out(MESSAGE3) << "Created clusterDB with coverage " << 100 * coverage << " percent on training data sample." << std::endl; 
      
      streamlog_out(MESSAGE3) << "ClusterDB written to file " << _clusterDBFileName << std::endl; 
      
      // Close clusterDB  file
      _rootClusterDBFile->Write();
      _rootClusterDBFile->Close();
      delete _rootClusterDBFile;
      
    }
  }
  
  //
  // Method printing processor parameters
  //
  void GoeClusterCalibrator::printProcessorParams() const 
  {
    
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "GoeClusterCalibrator Development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
    
  }
  
} // Namespace



