// PixelClusterizer implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "PixelClusterizer.h"

// Include TBTools 
#include "DEPFET.h" 
#include "TBDetector.h"

// Include ROOT classes
#include <TFile.h>

#include <iomanip>
using namespace std::string_literals;
// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;

namespace depfet {

//
// Instantiate this object
//
PixelClusterizer aPixelClusterizer ;

//
// Constructor
//
PixelClusterizer::PixelClusterizer() : Processor("PixelClusterizer"),_inputDecodeHelper(""),
    _orginalOutputEncoderHelper(DEPFET::ZSCLUSTERDEFAULTENCODING),_clusterOutputEncoderHelper(DEPFET::ZSCLUSTERDEFAULTENCODING)
{

// Processor description
   _description = "PixelClusterizer: Looking for clusters in zs pixel data" ;

//   
// First of all, we need to register the input/output collections
   
    registerInputCollection (LCIO::TRACKERDATA, "SparseDataCollectionName",
                            "Name of input sparsified pixel data collection",
                            _sparseDataCollectionName, string("zsdata"));
    
    registerOutputCollection (LCIO::TRACKERPULSE, "ClusterCollectionName",
                            "Name of the output cluster collection",
                            _clusterCollectionName, string("zscluster"));
    
    registerProcessorParameter( "SparseZSCut","Threshold for zero suppression",
                               _sparseZSCut, static_cast<float > (0));

    registerProcessorParameter( "SparseSeedCut","Threshold for seed pixel signal",
                               _sparseSeedCut, static_cast<float > (0));
   
    registerProcessorParameter( "SparseClusterCut","Threshold for cluster signal",
                               _sparseClusterCut, static_cast<float > (0) ); 
   
    registerProcessorParameter( "AcceptDiagonalCluster","0: common side; 1: common corner; 3: max. one missing pixel",
                               m_acceptDiagonalClusters , static_cast<int > (1) ); 

    registerProcessorParameter( "DifferenceTimeCut","Max absolute time difference between neighbor digits",
                               _absoluteNeighborTimeCut, static_cast<float > (0) ); 
  
    registerProcessorParameter("NoiseDBFileName",
                               "This is the name of the ROOT file with the status mask (add .root)",
                               _noiseDBFileName, static_cast< string > ( "NoiseDB.root" ) ); 
    
}

//
// Method called at the beginning of data processing
//
void PixelClusterizer::init() {
   
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   _hitElements = 4;
   
   _dummyCollectionName = "original_data_"+_clusterCollectionName;
   
   // Open noiseDB file 
   TFile * noiseDBFile = new TFile(_noiseDBFileName.c_str(), "READ");
   if (!noiseDBFile->IsZombie()) { 
     for(int ipl=0;ipl<TBDetector::GetInstance().GetNSensors();ipl++)  { 
       int sensorID = TBDetector::Get(ipl).GetSensorID();  
       string histoName = "hDB_sensor"+to_string(sensorID) + "_mask";
       if ( (TH2F *) noiseDBFile->Get(histoName.c_str()) != nullptr) {
         _DB_Map_Mask[sensorID] = (TH2F *) noiseDBFile->Get(histoName.c_str());  
         _DB_Map_Mask[sensorID]->SetDirectory(0);
       }  
     }
     // Close root  file
     noiseDBFile->Close();
   }
   delete noiseDBFile;
           
   // Print set parameters
   printProcessorParams();
   
   // CPU time start
   _timeCPU = clock()/1000;
   
}

//
// Method called for each run
//
void PixelClusterizer::processRunHeader(LCRunHeader * run)
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
void PixelClusterizer::processEvent(LCEvent * evt)
{
   
   // Print event number
   if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                   << (evt->getEventNumber())
                                                                   << std::endl << std::endl;
   
   // More detailed event numbering for testing
   streamlog_out(MESSAGE2) << std::endl << "Starting with Event Number " << evt->getEventNumber()  << std::endl;  
    
   //
   // Open collections
   try {
     
     // Output collection containing cluster data 
     LCCollectionVec * clusterCollection = new LCCollectionVec(LCIO::TRACKERPULSE) ;
     
     // Build clusters, i.e. adjecant groups of sparsified pixels passing certain signal cuts  
     clusterize( evt , clusterCollection  ); 
     
     // Add clusterCollection to event
     evt->addCollection( clusterCollection, _clusterCollectionName );
             	    
   } catch(DataNotAvailableException &e){}  
   _nEvt ++ ;
}


//
// Method called after each event to check the data processed
//
void PixelClusterizer::check( LCEvent * )
{
}

//
// Method called after all data processing
//
void PixelClusterizer::end()
{
   
   // CPU time end
   _timeCPU = clock()/1000 - _timeCPU;
   
   // Print message
   streamlog_out(MESSAGE) << std::endl
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
   
   
   
}


//
// Method printing processor parameters
//
void PixelClusterizer::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "PixelClusterizer Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}

// Called by the processEvent() once for the  
// pixel detector.  

void PixelClusterizer::clusterize( LCEvent * evt , LCCollectionVec * clusterCollection  ) 
{  
    
  // Open zero suppressed pixel data  
  LCCollectionVec * Pix_collection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_sparseDataCollectionName)); 
  // Helper class for decoding pixel data 
  CellIDDecoder<TrackerDataImpl> PixelID( Pix_collection,&_inputDecodeHelper );
   
  // The original data collection contains the sparse pixel data for 
  // each accpeted cluster.
  LCCollectionVec * originalDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);
  // Helper class for encoding clsuters  
  CellIDEncoder<TrackerDataImpl> originalDataEncoder( DEPFET::ZSCLUSTERDEFAULTENCODING, originalDataCollection,&_orginalOutputEncoderHelper );
  CellIDEncoder<TrackerPulseImpl> clusterEncoder(DEPFET::ZSCLUSTERDEFAULTENCODING, clusterCollection,&_clusterOutputEncoderHelper );

  // Clustering loop over pixel detectors 
  for (unsigned int iDet = 0; iDet < Pix_collection->size(); iDet++) { 
     
    // Get zs pixels from next pixel detector   
    TrackerDataImpl * pixModule = dynamic_cast<TrackerDataImpl* > ( Pix_collection->getElementAt(iDet) );

    // Sensor ID for pixel detector
    int sensorID = PixelID( pixModule ) ["sensorID"s];
    
    // Read geometry info for sensor 
    int ipl = TBDetector::GetInstance().GetPlaneNumber(sensorID);      
    const Det& adet = TBDetector::Get(ipl);

    // Get min channel numbers
    int minUCell = adet.GetMinUCell();
    int minVCell = adet.GetMinVCell();    
    // Get max channel numbers 
    int maxUCell = adet.GetMaxUCell();   
    int maxVCell = adet.GetMaxVCell(); 
    
    // List of firing pixels. Each pixel has a iU, iV, charge and time 
    FloatVec pixVector = pixModule->getChargeValues();
    int npixels = pixVector.size()/_hitElements; 
        
    // All firing pixels will be accumulated in groups of geometrical adjecant 
    // pixels, which are maintained in a vector of groups. These groups are dynamic
    // and can be merged, splitted and created in clustering.
    
    // Vector of pixel groups     
    Pix_GroupVector pixGroups;
    
    // Max. number of groups
    pixGroups.reserve(npixels);
        
    // Loop over zspixels and build groups of 
    // adjecant pixels. 
    for (int iPix = 0; iPix < npixels; iPix++) 
    {   
      
      int iU = static_cast<int> (pixVector[iPix * _hitElements]);
      int iV = static_cast<int> (pixVector[iPix * _hitElements + 1]);
      float charge =  pixVector[iPix * _hitElements + 2];  
      float time =  pixVector[iPix * _hitElements + 3];     
       
      // Try to get status code for pixel 
      float status = 0; 
      if ( _DB_Map_Mask.find(sensorID) != _DB_Map_Mask.end() ) {
        // Here we use the same numbering convention as in the HotPixelKiller processor to map iu, iv to a bin in the mask. 
        status = _DB_Map_Mask[sensorID]->GetBinContent(iU-minUCell+1, iV-minVCell+1); 
      }
      
      // Print detailed pixel summary, for testing/debugging only !!! 
      streamlog_out(MESSAGE1) << "Pixel Nr. " << iPix << " on sensor " << sensorID  
                              << std::endl;  
      streamlog_out(MESSAGE1) << "   iU:" << iU << ", iV:" << iV
                              << ", charge:" << charge
                              << ", time:" << time
                              << std::endl;
      
      // If pixel signal below threshold, skip it in clusterization
      if ( charge < _sparseZSCut ) {
        streamlog_out(MESSAGE2) << "  Signal below ZS cut. Skipping it." << std::endl; 
        continue;
      }
       
      // If a pixel is out of range, skip it in clusterization
      if ( iU < minUCell || iU > maxUCell || iV < minVCell || iV > maxVCell ) {
        streamlog_out(MESSAGE2) << "  Invalid pixel address found. Skipping it." << std::endl; 
        continue;
      }
      
      // Skip masked pixels 
      if ( status  != 0 ) {
        streamlog_out(MESSAGE2) << "  Bad pixel found. Skipping it." << std::endl; 
        continue;
      }
      
      // Loop on all existing pixel groups, until you find that current pixel 
      // is a neighbour of a group.
      bool found= false;
         	
      Pix_GroupVector::iterator firstGroup = pixGroups.begin();
      Pix_GroupVector::iterator lastGroup  = pixGroups.end();    
              
      while( !found && firstGroup!= lastGroup)
      {
        
        if ( areNeighbours( *firstGroup, iU, iV, ipl, time) )
        {
           
          // If pixel is a duplicate of one in the cluster, do not add it.   
          if(!isDuplicated( *firstGroup, iU, iV )){
            
            // Add this pixel to this pixel group 
            (*firstGroup).push_back(iU);
            (*firstGroup).push_back(iV);
            (*firstGroup).push_back(charge);
            (*firstGroup).push_back(time);
            
            // See if iU/iV is a neighbour to any other groups, if yes perform merging 
            checkForMerge(iU, iV, ipl, time, firstGroup, lastGroup);
              
          } else {
            streamlog_out(MESSAGE2) << "  A pixel duplicate found. Skipping it." << std::endl; 
          }
          found = true; 
        }
        ++firstGroup;
      }
      
      // If pixel is isolated, seed a new candidate cluster. 
      if(!found)
      {
        FloatVec newGroup;
        newGroup.reserve(6);
        newGroup.push_back(iU);
        newGroup.push_back(iV);
        newGroup.push_back(charge);
        newGroup.push_back(time);
        
        pixGroups.push_back(std::move(newGroup));
        
      }
       
    }  
    
    // Now, we make a quality selection among pixel groups 
    streamlog_out(MESSAGE2) << std::endl << "Summary of cluster candidates: "  << std::endl << std::endl;
    
      
    
    // Count pixel groups 
    int groupNumber = 0;
    
    // Counts good clusters
    int clusterID = 0; 
    
    for( Pix_GroupVector::iterator group = pixGroups.begin() ;
     	group!= pixGroups.end() ; ++group) 
    {
      
      groupNumber++; 
      
      // If cluster is empty, i.e. it has been merged with another, 
      // do not attempt to make cluster.
      if ((*group).size()> 0)
      {
          
        // We have found a new cluster candidate. We loop over all pixels  
        // to calculate S/N ratios. 
        
        int cluQuality = 0;
        float clusterSignal = 0; 
        float seedSignal = 0; 
        
        
        // Prepare a TrackerData to store original data of cluster
        TrackerDataImpl* sparseCluster = new TrackerDataImpl ; 
        
        int npixels = (*group).size()/_hitElements; 
        auto& chargeVec=sparseCluster->chargeValues();
        chargeVec.reserve((*group).size());
        for ( int index=0; index < npixels; index++)
        {
           
          int iU = static_cast<int> ( (*group)[index * _hitElements]);
          int iV = static_cast<int> ( (*group)[index * _hitElements + 1]);
          float charge = ( (*group)[index * _hitElements + 2]);
          float time = ( (*group)[index * _hitElements + 3]);
                   
          clusterSignal += charge;
         
          if (charge > seedSignal ) {
            seedSignal = charge; 
          }
             
          // Store pixel data int EUTelescope format 
          chargeVec.push_back( iU );
          chargeVec.push_back( iV );
          chargeVec.push_back( charge );
          chargeVec.push_back( time );
          
        }
            
        // Use a global ADU threshold 
        float seedThreshold = _sparseSeedCut; 
        float clusterThreshold = _sparseClusterCut;
          
        // Verify if the cluster candidates is a good cluster
        if ( ( seedSignal >= seedThreshold ) && 
             ( clusterSignal >= clusterThreshold ) ) {
          
                
 	      // New accpeted cluster
          clusterID++; 
            
          streamlog_out(MESSAGE2) << " Stored cluster on sensor " << sensorID << " having total charge " << clusterSignal << std::endl;
          static auto idx_orig_sensorID=originalDataEncoder.index("sensorID"s); //find the address ONCE.
          static auto idx_orig_clusterID=originalDataEncoder.index("clusterID"s);
          static auto idx_orig_sparsePixelType=originalDataEncoder.index("sparsePixelType"s);
          static auto idx_orig_quality=originalDataEncoder.index("quality"s);
          // Ok good cluster ... save it   
          originalDataEncoder[idx_orig_sensorID] = sensorID;
          originalDataEncoder[idx_orig_clusterID] = 0;
          originalDataEncoder[idx_orig_sparsePixelType] = static_cast<int> (kSimpleSparsePixel);
          originalDataEncoder[idx_orig_quality] = cluQuality;
          originalDataEncoder.setCellID( sparseCluster );
          originalDataCollection->push_back( sparseCluster );
                             
          static auto idx_clust_sensorID=clusterEncoder.index("sensorID"s);
          static auto idx_clust_clusterID=clusterEncoder.index("clusterID"s);
          static auto idx_clust_sparsePixelType=clusterEncoder.index("sparsePixelType"s);
          static auto idx_clust_quality=clusterEncoder.index("quality"s);
          TrackerPulseImpl* zsPulse = new TrackerPulseImpl;
          clusterEncoder[idx_clust_sensorID]  = sensorID;
          clusterEncoder[idx_clust_clusterID] = 0;
          clusterEncoder[idx_clust_sparsePixelType] = static_cast<int> (kSimpleSparsePixel);
          clusterEncoder[idx_clust_quality] = cluQuality;
          clusterEncoder.setCellID( zsPulse );
          
          zsPulse->setCharge( clusterSignal );
          zsPulse->setQuality( cluQuality );
          zsPulse->setTrackerData( sparseCluster );
          clusterCollection->push_back( zsPulse );          
          
           
          
        } else {
          // Cluster candidate thrown away, so clean up memory    
          sparseCluster->chargeValues().clear(); 
          delete sparseCluster;
        }
        
      } // Non empty cluster candidate 
    } // Cluster candidate loop 
    
    // 
    // Free used memory 
     
    for( Pix_GroupVector::iterator group = pixGroups.begin() ;
      	group!= pixGroups.end() ; ++group) 
    {
      (*group).clear(); 
    }
    pixGroups.clear(); 
    
  } // End detector loop 
  
  // Add original data collection to event
  evt->addCollection( originalDataCollection, _dummyCollectionName );
  
}



// Checks if any other pixel group (apart from base group) neighbours iU/iV. 
// If so, merge with base group.  
 
void PixelClusterizer::checkForMerge( int iU, int iV,  int planeNumber, float time,
 Pix_GroupVector::iterator baseGroup,
 Pix_GroupVector::iterator lastGroup) 
{
  
  // First non-tested group
  Pix_GroupVector::iterator nextGroup(baseGroup+1);
   
  for (; nextGroup!= lastGroup; ++nextGroup)
  {              
    if (areNeighbours( *nextGroup, iU, iV, planeNumber, time ))
    {
      // Merge these pixel groups
      int npixels = (*nextGroup).size()/_hitElements; 
       
      for ( int index=0; index < npixels; index++)      
      {
        float iUtmp = (*nextGroup)[index * _hitElements];
        float iVtmp = (*nextGroup)[index * _hitElements + 1];
        float chargetmp =  (*nextGroup)[index * _hitElements + 2];    
        float timetmp =  (*nextGroup)[index * _hitElements + 3];     
        // Copy pixel to base group 
        (*baseGroup).push_back(iUtmp);
        (*baseGroup).push_back(iVtmp);
        (*baseGroup).push_back(chargetmp);
        (*baseGroup).push_back(timetmp);
      }
      (*nextGroup).clear();
       
    }
  }
}


// This method is called inside the clusterize() method in order to 
// determine if the pixel cell at address (iU,iV) should be added 
// to the candidate cluster passed as first argument.  

bool PixelClusterizer::areNeighbours( FloatVec &group, int iU, int iV, int planeNumber, float time ) 
{   
  int npixels = group.size()/_hitElements; 
  
  for ( int index=0; index < npixels; index++)
  {
           
    int iU1 = static_cast<int> (group[index * _hitElements]);
    int iV1 = static_cast<int> (group[index * _hitElements + 1]);
    int time1 = static_cast<int> (group[index * _hitElements + 3]);
 
    if(TBDetector::Get(planeNumber).areNeighbors(iV1, iU1, iV, iU) && std::abs(time-time1) <= _absoluteNeighborTimeCut ) return true;               
  }
    
  return false;
}


// This method is called inside the clusterize() method in order to 
// determine if the pixel cell with addess iU/iV is already part 
// of cluster candidate passed as first argument.  
                                        
bool PixelClusterizer::isDuplicated( FloatVec &group, int iU, int iV ) 
{ 
  bool duplicate = false;
  int npixels = group.size()/_hitElements; 
   
  for ( int index=0; index < npixels; index++)
  {
    int iU1 = static_cast<int> (group[index * _hitElements]);
    int iV1 = static_cast<int> (group[index * _hitElements + 1]);  
     
    // Duplicate?? 
    if(iV1 == iV && iU1 == iU){
      duplicate = true;
      break; 
    }
  }  
  return duplicate;
}

} // Namespace



