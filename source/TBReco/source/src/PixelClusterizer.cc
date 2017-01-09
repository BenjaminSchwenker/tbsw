// PixelClusterizer implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "PixelClusterizer.h"

// Include DEPFETTrackTools 
#include "DEPFET.h" 
#include "MatrixDecoder.h"

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
PixelClusterizer::PixelClusterizer() : Processor("PixelClusterizer")
{

// Processor description
   _description = "PixelClusterizer: Looking for clusters in zs pixel data" ;

//   
// First of all, we need to register the input/output collections
   
   registerInputCollection (LCIO::TRACKERDATA, "SparseDataCollectionName",
                            "Name of input sparsified pixel data collection",
                            _sparseDataCollectionName, string("zsdata"));
   
   registerInputCollection (LCIO::TRACKERRAWDATA, "StatusCollectionName",
                           "Input status collection name",
                           _statusCollectionName, string("status")); 
    
   registerOutputCollection (LCIO::TRACKERPULSE, "ClusterCollectionName",
                            "Name of the output cluster collection",
                            _clusterCollectionName, string("zscluster"));
    
   registerProcessorParameter( "SparseZSCut","Threshold for zero suppression",
                               _sparseZSCut, static_cast<float > (0));

   registerProcessorParameter( "SparseSeedCut","Threshold for seed pixel signal",
                               _sparseSeedCut, static_cast<float > (5));
   
   registerProcessorParameter( "SparseClusterCut","Threshold for cluster signal",
                               _sparseClusterCut, static_cast<float > (7) ); 
   
   registerProcessorParameter( "AcceptDiagonalCluster","0: common side; 1: common corner; 3: max. one missing pixel",
                               m_acceptDiagonalClusters , static_cast<int > (1) ); 
    
}

//
// Method called at the beginning of data processing
//
void PixelClusterizer::init() {
   
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   
   // Read detector constants from gear file
   _detector.ReadGearConfiguration();    
   
   _dummyCollectionName = "original_data_"+_clusterCollectionName;
   
   // The status data is not yet initialized
   _isStatusReady = false;
              
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
   
   // First of all we need to be sure that status data is initialized
   if ( !_isStatusReady ) {
     initializeStatus( evt ) ;
   }
   
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
void PixelClusterizer::check( LCEvent * evt )
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
  CellIDDecoder<TrackerDataImpl> PixelID( Pix_collection );  
  
  //Open status data
  LCCollectionVec * statusCollection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_statusCollectionName));
  // Helper class for decoding status data 
  CellIDDecoder<TrackerDataImpl> statusDecoder( statusCollection ); 
  
  // The original data collection contains the sparse pixel data for 
  // each accpeted cluster.
  LCCollectionVec * originalDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);
  
  // Clustering loop over pixel detectors 
  for (unsigned int iDet = 0; iDet < Pix_collection->size(); iDet++) { 
     
    // Get zs pixels from next pixel detector   
    TrackerDataImpl * pixModule = dynamic_cast<TrackerDataImpl* > ( Pix_collection->getElementAt(iDet) );
     
    // Get pixel status matrix 
    TrackerRawDataImpl * status = dynamic_cast<TrackerRawDataImpl*> (statusCollection->getElementAt( iDet ));  
      
    // DAQ ID for pixel detector
    int sensorID = PixelID( pixModule ) ["sensorID"];

    // Read geometry info for sensor 
    int ipl = _detector.GetPlaneNumber(sensorID);      
    Det& adet = _detector.GetDet(ipl);
  
    int noOfXPixels = adet.GetNColumns(); 
    int noOfYPixels = adet.GetNRows();

    MatrixDecoder matrixDecoder( noOfXPixels, noOfYPixels); 
    
    // Get max channel numbers 
    int maxCol = adet.GetNColumns() - 1; 
    int maxRow = adet.GetNRows() - 1;
    
    // List of firing pixels. Each pixel has a col, row and charge 
    FloatVec pixVector = pixModule->getChargeValues();
    int npixels = pixVector.size()/3; 
        
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
      
      int col = static_cast<int> (pixVector[iPix * 3]);
      int row = static_cast<int> (pixVector[iPix * 3 + 1]);
      float charge =  pixVector[iPix * 3 + 2];     
      
      // Print detailed pixel summary, for testing/debugging only !!! 
      streamlog_out(MESSAGE1) << "Pixel Nr. " << iPix << " on sensor " << sensorID  
                              << std::endl;  
      streamlog_out(MESSAGE1) << "   column:" << col << ", row:" << row
                              << ", charge:" << charge
                              << std::endl;
      
      // If pixel signal below threshold, skip it in clusterization
      if ( charge < _sparseZSCut ) {
        streamlog_out(MESSAGE2) << "  Signal below ZS cut. Skipping it." << std::endl; 
        continue;
      }
       
      // If a pixel is out of range, skip it in clusterization
      if ( col < 0 || col > maxCol || row < 0 || row > maxRow ) {
        streamlog_out(MESSAGE2) << "  Invalid pixel address found. Skipping it." << std::endl; 
        continue;
      }
      
      // If a pixel is bad, skip it in clusterization
      if ( status->getADCValues() [ matrixDecoder.getIndexFromXY(col,row) ]   != 0 ) {
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
        
        if ( areNeighbours( *firstGroup, col, row, m_acceptDiagonalClusters ) )
        {
           
          // If pixel is a duplicate of one in the cluster, do not add it.   
          if(!isDuplicated( *firstGroup, col, row )){
            
            // Add this pixel to this pixel group 
            (*firstGroup).push_back(col);
            (*firstGroup).push_back(row);
            (*firstGroup).push_back(charge);
            
            // See if col/row is a neighbour to any other groups, if yes perform merging 
            checkForMerge(col, row, firstGroup, lastGroup);
              
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
        newGroup.push_back(col);
        newGroup.push_back(row);
        newGroup.push_back(charge);
        
        pixGroups.push_back(newGroup);
        
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
        
        int npixels = (*group).size()/3; 
        
        for ( int index=0; index < npixels; index++)
        {
           
          int col = static_cast<int> ( (*group)[index * 3]);
          int row = static_cast<int> ( (*group)[index * 3 + 1]);
          float charge = ( (*group)[index * 3 + 2]);
                   
          clusterSignal += charge;
         
          if (charge > seedSignal ) {
            seedSignal = charge; 
          }
             
          // Store pixel data int EUTelescope format 
          sparseCluster->chargeValues().push_back( col );
          sparseCluster->chargeValues().push_back( row );
          sparseCluster->chargeValues().push_back( charge );   
          
        }
            
        // Use a global ADU threshold 
        float seedThreshold = _sparseSeedCut; 
        float clusterThreshold = _sparseClusterCut;
          
        // Verify if the cluster candidates is a good cluster
        if ( ( seedSignal >= seedThreshold ) && 
             ( clusterSignal >= clusterThreshold ) ) {
        
                
 	  // New accpeted cluster
          clusterID++; 
          
          // Ok good cluster ... save it 
          CellIDEncoder<TrackerDataImpl> originalDataEncoder( DEPFET::ZSCLUSTERDEFAULTENCODING, originalDataCollection );
          CellIDEncoder<TrackerPulseImpl> clusterEncoder(DEPFET::ZSCLUSTERDEFAULTENCODING, clusterCollection ); 
           
          // Set the ID for this zsCluster
          originalDataEncoder["sensorID"] = sensorID;
          originalDataEncoder["clusterID"] = 0;
          originalDataEncoder["sparsePixelType"] = static_cast<int> (kSimpleSparsePixel);
          originalDataEncoder["quality"] = cluQuality;
          originalDataEncoder.setCellID( sparseCluster );
          originalDataCollection->push_back( sparseCluster );
                             
          TrackerPulseImpl* zsPulse = new TrackerPulseImpl;
          clusterEncoder["sensorID"]  = sensorID;
          clusterEncoder["clusterID"] = 0;
          clusterEncoder["sparsePixelType"] = static_cast<int> (kSimpleSparsePixel);
          clusterEncoder["quality"] = cluQuality;
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



// Checks if any other pixel group (apart from base group) neighbours col/row. 
// If so, merge with base group.  
 
void PixelClusterizer::checkForMerge( int col, int row,
 Pix_GroupVector::iterator baseGroup,
 Pix_GroupVector::iterator lastGroup) 
{
  
  // First non-tested group
  Pix_GroupVector::iterator nextGroup(baseGroup+1);
   
  for (; nextGroup!= lastGroup; ++nextGroup)
  {              
    if (areNeighbours( *nextGroup, col, row, m_acceptDiagonalClusters ))
    {
      // Merge these pixel groups
      int npixels = (*nextGroup).size()/3; 
       
      for ( int index=0; index < npixels; index++)      
      {
        float coltmp = (*nextGroup)[index * 3];
        float rowtmp = (*nextGroup)[index * 3 + 1];
        float chargetmp =  (*nextGroup)[index * 3 + 2];     
        // Copy pixel to base group 
        (*baseGroup).push_back(coltmp);
        (*baseGroup).push_back(rowtmp);
        (*baseGroup).push_back(chargetmp);
      }
      (*nextGroup).clear();
       
    }
  }
}


// This method is called inside the clusterize() method in order to 
// determine if the pixel cell at address (col,row) should be added 
// to the candidate cluster passed as first argument.  
//
// Different clustering strategies are foreseen: 
// 
// m_accept Clusters
//   = 0: Accept pixels which have a side in common with a pixel cell 
//        in the list
//   = 1: A common corner suffices (default setting)
//   = 2: Max distance is a missing pixel in a row or a column
//   = 3: Max distance is a missing diagonal pixel 


bool PixelClusterizer::areNeighbours( FloatVec &group, int col, int row, int m_accept ) 
{   
    
  bool match=false;
  int npixels = group.size()/3; 
  
  for ( int index=0; index < npixels; index++)
  {
           
    int col1 = static_cast<int> (group[index * 3]);
    int row1 = static_cast<int> (group[index * 3 + 1]);
    
    int deltarow = abs(row-row1);
    int deltacol = abs(col-col1);
          
    // A side in common
    if(deltacol+deltarow < 2) match = true;
     
    // A corner in common 
    if(m_accept == 1 && deltacol == 1 
                                       && deltarow == 1) match = true;
     
    // max distance is 2 pixels (includes cases with missing pixels)
    if(m_accept == 2 && deltacol+deltarow < 3 ) match = true;
     
    // max distance is 2 pixels along a diagonal (includes cases with missing pixels)
    if(m_accept == 3 && deltacol < 3 && deltarow < 3) match = true;
                      
  }
    
  return match;
}


// This method is called inside the clusterize() method in order to 
// determine if the pixel cell with addess col/row is already part 
// of cluster candidate passed as first argument.  
                                        
bool PixelClusterizer::isDuplicated( FloatVec &group, int col, int row ) 
{ 
  bool duplicate = false;
  int npixels = group.size()/3; 
   
  for ( int index=0; index < npixels; index++)
  {
    int col1 = static_cast<int> (group[index * 3]);
    int row1 = static_cast<int> (group[index * 3 + 1]);  
     
    // Duplicate?? 
    if(row1 == row && col1 == col){
      duplicate = true;
      break; 
    }
  }  
  return duplicate;
}


// This method is called to initialize the status information,
// namely linking status data  to pixel data.
 
void  PixelClusterizer::initializeStatus( LCEvent * event )  {
  
  streamlog_out( MESSAGE2 ) << "Initializing status" << endl;
  
  try {
         
    // Open pixel data  
    LCCollectionVec * Pix_collection = dynamic_cast < LCCollectionVec * > (event->getCollection(_sparseDataCollectionName)); 
    CellIDDecoder<TrackerDataImpl> PixelID( Pix_collection ); 
    
    //Open status data
    LCCollectionVec * statusCollection = dynamic_cast < LCCollectionVec * > (event->getCollection(_statusCollectionName)); 
    CellIDDecoder<TrackerRawDataImpl> statusDecoder( statusCollection ); 
    
    // We are assuming that the pixel and status collections 
    // are aligned according to the sensorID. In other words, we are 
    // assuming the element i-th in the all the collections corresponds 
    // to the same detector. Let's test this ...

    streamlog_out( MESSAGE3 ) << "Found status collection " << endl;
    
    for ( size_t iDet = 0 ; iDet < statusCollection->size(); ++iDet ) {

      streamlog_out( MESSAGE3 ) << "Found iDet " << iDet << endl;
       
      TrackerDataImpl * pixel = dynamic_cast<TrackerDataImpl* > ( Pix_collection->getElementAt(iDet) )  ;
      TrackerRawDataImpl * status = dynamic_cast < TrackerRawDataImpl * >(statusCollection->getElementAt(iDet));
      
      int dataID = static_cast<int> ( PixelID(pixel)["sensorID"] );
      int statusID = static_cast<int> ( statusDecoder(status)["sensorID"] );
            
      if (dataID != statusID) {
        streamlog_out(ERROR3) << "Status collection not aligned to detector data!" << std::endl << std::endl;      
        exit(-1); 
      }
      
    }
    
    _isStatusReady = true;
    
  } catch (  lcio::DataNotAvailableException ) {
    streamlog_out( MESSAGE2 ) << "Unable to to find status collection" << endl;
    _isStatusReady = false;
  }
   
}

} // Namespace



