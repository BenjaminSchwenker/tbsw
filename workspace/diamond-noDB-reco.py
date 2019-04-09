# -*- coding: utf-8 -*-

# script to analyse diamond testbeam data (planar). Not using DB cluster

import tbsw 
import multiprocessing
import os

# Path to steering files
steerfiles = 'steering-files/diamond-tb/'

# parameters to set for an analysis
biasVoltage = 300 # V
beamEnergy = 4.6 # GeV
gearBaseName = 'gear_dia_planar_batch'
gearBatch = 2 # number of the batch

usePseudoPixel = False
useClusterDB = False


#parameters
# group of: (rawfile, bias voltage, beam energy, gearfile batch, maxrecord number)
runConfiguration = [
                    ('/work2/H.beck/testbeamNov2017/data/100V/run000998.raw', 100, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/100V/run001001.raw', 100, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/100V/run001002.raw', 100, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/100V/run001007.raw', 100, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/100V/run001011.raw', 100, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/100V/run001012.raw', 100, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/100V/run001028.raw', 100, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/100V/run001034.raw', 100, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/150V/run001084.raw', 150, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/150V/run001086.raw', 150, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/150V/run001087.raw', 150, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/150V/run001088.raw', 150, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/150V/run001090.raw', 150, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/150V/run001091.raw', 150, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/150V/run001093.raw', 150, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/150V/run001094.raw', 150, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/150V/run001095.raw', 150, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/150V/run001097.raw', 150, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/150V/run001098.raw', 150, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/150V/run001099.raw', 150, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/200V/run0010.raw', 200, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/200V/run0010.raw', 200, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/200V/run0010.raw', 200, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/200V/run0010.raw', 200, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/200V/run0010.raw', 200, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/200V/run0010.raw', 200, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/200V/run0010.raw', 200, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/200V/run0010.raw', 200, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000732.raw', 300, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000733.raw', 300, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000735.raw', 300, 4.6, 2, 100000),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000736.raw', 300, 4.6, 2, 100000),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000738.raw', 300, 4.6, 2, 200000),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000739.raw', 300, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000740.raw', 300, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000742.raw', 300, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000744.raw', 300, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000747.raw', 300, 4.6, 2, 100000),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000748.raw', 300, 4.6, 2, 180000),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000750.raw', 300, 4.6, 2, 210000),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000752.raw', 300, 4.6, 2, 140000),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000753.raw', 300, 4.6, 2, 200000),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000754.raw', 300, 4.6, 2, 150000),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000755.raw', 300, 4.6, 2, 220000),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000756.raw', 300, 4.6, 2, 190000),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000758.raw', 300, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000759.raw', 300, 4.6, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000763.raw', 300, 4.6, 2, 230000),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000764.raw', 300, 4.6, 2, 210000),
                    ('/work2/H.beck/testbeamNov2017/data/300V/run000766.raw', 300, 4.6, 2, 130000),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000961.raw', 500, 4.0, 2, 90000),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000963.raw', 500, 4.0, 2, 100000),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000964.raw', 500, 4.0, 2, 80000),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000965.raw', 500, 4.0, 2, 60000),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000966.raw', 500, 4.0, 2, 80000),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000968.raw', 500, 4.0, 2, 60000),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000969.raw', 500, 4.0, 2, 60000),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000970.raw', 500, 4.0, 2, 120000),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000976.raw', 500, 4.0, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000978.raw', 500, 4.0, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000980.raw', 500, 4.0, 2, 130000),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000981.raw', 500, 4.0, 2, 110000),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000982.raw', 500, 4.0, 2, 120000),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000985.raw', 500, 4.0, 2, 130000),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000987.raw', 500, 4.0, 2, 110000),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000988.raw', 500, 4.0, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000989.raw', 500, 4.0, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000991.raw', 500, 4.0, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/500V/run000992.raw', 500, 4.0, 2, -1),
                    ('/work2/H.beck/testbeamNov2017/data/600V/run001101.raw', 600, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/600V/run001103.raw', 600, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/600V/run001112.raw', 600, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/600V/run001115.raw', 600, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/600V/run001116.raw', 600, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/600V/run001117.raw', 600, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/600V/run001118.raw', 600, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/600V/run001119.raw', 600, 4.0, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/600V/run001121.raw', 600, 4.6, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/600V/run001122.raw', 600, 4.6, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/600V/run001132.raw', 600, 4.6, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/600V/run001133.raw', 600, 4.6, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/600V/run001134.raw', 600, 4.6, 4, -1),
                    ('/work2/H.beck/testbeamNov2017/data/600V/run001135.raw', 600, 4.6, 4, -1),
                   ]

# List of runs to be processed
runlist = [runConfig[0] for runConfig in runConfiguration if runConfig[1] == biasVoltage and runConfig[2]== beamEnergy and runConfig[3]==gearBatch][0]

nEventslist = [runConfig[4] for runConfig in runConfiguration if runConfig[1] == biasVoltage and runConfig[2]== beamEnergy and runConfig[3]==gearBatch][0]

# Run which is used to calibrate all other runs
calibrun = [runConfig[0] for runConfig in runConfiguration if runConfig[1] == biasVoltage and runConfig[2]== beamEnergy and runConfig[3]==gearBatch and runConfig[4]==-1][0]

# Min number of clusters or cut away
minCluster = 200

# Filter Channel: 20 + running number of UBSPix readout dut assignment
filterChannelRef = 1
filterChannelDia = 2

singleHitSeeding = "0 1" # needed because variable later used

if gearBatch == 2:
  ignoreIDsRef = "0 1 2 3 4 5 21" # FEI4ClusterCalibrator
  ignoreIDsDia = "0 1 2 3 4 5 22" # DiamondClusterCalibrator
  filterChannelRef = 2 # FEI4Unpacker
  filterChannelDia = 1 # DiamondUnpacker

elif gearBatch == 4:
  ignoreIDsRef = "0 1 2 3 4 5 22" # FEI4ClusterCalibrator
  ignoreIDsDia = "0 1 2 3 4 5 21" # DiamondClusterCalibrator

# end of global variables

def create_calibration_path(Env, rawfile, gearfile):
  """
  Returns a list of tbsw path objects to calibrate the tracking telescope
  """
  
  # Calibrations are organized in a sequence of calibration paths.
  # The calibration paths are collected in a list for later execution
  calpaths = []

  # Create path for detector level masking of hot channels
  mask_path = Env.create_path('mask_path')
  mask_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1})

  rawinput = tbsw.Processor(name="RawInputProcessor", proctype="EudaqInputProcessor")
  rawinput.param("FileNames", rawfile)
  mask_path.add_processor(rawinput)

  geo = tbsw.Processor(name="Geo",proctype="Geometry")
  geo.param("AlignmentDBFilePath", "localDB/alignmentDB.root")
  geo.param("ApplyAlignment", "true")
  geo.param("OverrideAlignment", "true")
  mask_path.add_processor(geo)
  
  if usePseudoPixel: # not set up
    mask_path = add_unpackersPseudo(mask_path)
    diasorter = tbsw.Processor(name="DiamondRawHitSorter", proctype="DiamondRawHitSorter")
    diasorter.param("PixelType", diamond_pixeltype)
    mask_path.add_processor(diasorter)
  else:
    mask_path = add_unpackers(mask_path)

  # Noisy pixel 
  m26hotpixelkiller = tbsw.Processor(name="M26HotPixelKiller",proctype="HotPixelKiller",)
  m26hotpixelkiller.param("InputCollectionName", "zsdata_m26")
  m26hotpixelkiller.param("MaxOccupancy", 0.001)
  m26hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-M26.root")
  m26hotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(m26hotpixelkiller)

  fei4hotpixelkiller = tbsw.Processor(name="FEI4HotPixelKiller", proctype="HotPixelKiller")
  fei4hotpixelkiller.param("InputCollectionName", "zsdata_fei4")
  fei4hotpixelkiller.param("MaxOccupancy", 0.001)
  fei4hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-FEI4.root")
  fei4hotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(fei4hotpixelkiller)

  diahotpixelkiller = tbsw.Processor(name="DiamondHotPixelKiller", proctype="HotPixelKiller")
  diahotpixelkiller.param("InputCollectionName", "zsdata_dia")
  diahotpixelkiller.param("MaxOccupancy", 0.001)
  diahotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-DIA.root")
  diahotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(diahotpixelkiller)

  # Add path for masking
  calpaths.append(mask_path)

  # Create path for detector level creation of clusters
  clusterizer_path = Env.create_path('clusterizer_path')
  clusterizer_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  -1})
  clusterizer_path.add_processor(rawinput)
  clusterizer_path.add_processor(geo)
  
  if usePseudoPixel: #not set up
    clusterizer_path = add_unpackersPseudo(clusterizer_path)
    diasorter = tbsw.Processor(name="DiamondRawHitSorter", proctype="DiamondRawHitSorter")
    diasorter.param("PixelType", diamond_pixeltype) # in TH2Poly layout diamond_pixeltype no longer needed, need a more general sorter, which sorts all pixel in one pseudo declaration (no start again a t 0,0 for next type)
    clusterizer_path.add_processor(diasorter)
  else:
    clusterizer_path = add_unpackers(clusterizer_path)

  clusterizer_path = add_clusterizers(clusterizer_path)

  lciooutput = tbsw.Processor(name="LCIOOutput",proctype="LCIOOutputProcessor")
  lciooutput.param("LCIOOutputFile","tmp.slcio")
  lciooutput.param("LCIOWriteMode","WRITE_NEW")
  clusterizer_path.add_processor(lciooutput)

  # Finished with path for clusterizers
  calpaths.append(clusterizer_path)

  # Create path for pre alignmnet and dqm based on hits
  correlator_path = Env.create_path("correlator_path")
  correlator_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': "tmp.slcio" , 'Verbosity': "MESSAGE2"})
  correlator_path.add_processor(geo)
  
  correlator_path = add_hitmakers(correlator_path) # only cogHit in the start, cluster calibration later (if wanted)


  hitdqm = tbsw.Processor(name="RawDQM",proctype="RawHitDQM")
  hitdqm.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_dia")
  hitdqm.param("RootFileName","RawDQM.root")
  correlator_path.add_processor(hitdqm)

  correlator = tbsw.Processor(name="TelCorrelator", proctype="Correlator")
  correlator.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_dia")
  correlator.param("OutputRootFileName","XCorrelator.root")
  correlator.param("ReferencePlane","4")
  correlator.param("ParticleCharge","-1")
  correlator.param("ParticleMass","0.000511")
  correlator.param("ParticleMomentum", beamEnergy)
  correlator_path.add_processor(correlator)

  # Finished with path for hit based pre alignment
  calpaths.append(correlator_path)

  # Create path for pre alignment with loose cut track sample
  prealigner_path = Env.create_path('prealigner_path')
  prealigner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': "tmp.slcio" })
  prealigner_path.add_processor(geo)

  prealigner_path = add_hitmakers(prealigner_path)

  trackfinder_loosecut = tbsw.Processor(name="AlignTF_LC",proctype="FastTracker")
  trackfinder_loosecut.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_dia")
  trackfinder_loosecut.param("ExcludeDetector", "")
  trackfinder_loosecut.param("MaxTrackChi2", 10000000)
  trackfinder_loosecut.param("MaximumGap", 1)
  trackfinder_loosecut.param("MinimumHits",7)
  trackfinder_loosecut.param("OutlierChi2Cut", 100000000)
  trackfinder_loosecut.param("ParticleCharge","-1")
  trackfinder_loosecut.param("ParticleMass","0.000511")
  trackfinder_loosecut.param("ParticleMomentum", beamEnergy)
  trackfinder_loosecut.param("SingleHitSeeding", singleHitSeeding)
  trackfinder_loosecut.param("MaxResidualU","0.5")
  trackfinder_loosecut.param("MaxResidualV","0.5")
  prealigner_path.add_processor(trackfinder_loosecut)

  prealigner = tbsw.Processor(name="PreAligner",proctype="KalmanAligner")
  prealigner.param('ErrorsShiftX' , '0 10 10 10 10 10 10 0')
  prealigner.param('ErrorsShiftY' , '0 10 10 10 10 10 10 0')
  prealigner.param('ErrorsShiftZ' , '0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsAlpha'  , '0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsBeta'   , '0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsGamma'  , '0 0.01 0.01 0.01 0.01 0.01 0.01 0')
  prealigner_path.add_processor(prealigner)

  # Finished with path for prealigner
  calpaths.append(prealigner_path)


  # Create path for alignment with tight cut track sample
  aligner_path = Env.create_path('aligner_path')
  aligner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': "tmp.slcio" })
  aligner_path.add_processor(geo)

  aligner_path = add_hitmakers(aligner_path)

  # no dia in this, later, as in old script
  trackfinder_tightcut = tbsw.Processor(name="AlignTF_TC",proctype="FastTracker")
  trackfinder_tightcut.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_dia")
  trackfinder_tightcut.param("ExcludeDetector", "")
  trackfinder_tightcut.param("MaxTrackChi2", 100)
  trackfinder_tightcut.param("MaximumGap", 1)
  trackfinder_tightcut.param("MinimumHits",7)
  trackfinder_tightcut.param("OutlierChi2Cut", 10)
  trackfinder_tightcut.param("ParticleCharge","-1")
  trackfinder_tightcut.param("ParticleMass","0.000511")
  trackfinder_tightcut.param("ParticleMomentum", beamEnergy)
  trackfinder_tightcut.param("SingleHitSeeding", singleHitSeeding)
  trackfinder_tightcut.param("MaxResidualU","0.5")
  trackfinder_tightcut.param("MaxResidualV","0.5")
  aligner_path.add_processor(trackfinder_tightcut)

  aligner = tbsw.Processor(name="Aligner",proctype="KalmanAligner")
  aligner.param('ErrorsShiftX' , '0 10 10 10 10 10 10 0' )
  aligner.param('ErrorsShiftY' , '0 10 10 10 10 10 10 0')
  aligner.param('ErrorsShiftZ' , '0 10 10 10 10 10 10 0')
  aligner.param('ErrorsAlpha'  , '0 0 0 0 0 0 0 0') # errors for dut, ref? TODO
  aligner.param('ErrorsBeta'   , '0 0 0 0 0 0 0 0')
  aligner.param('ErrorsGamma'  , '0 0.01 0.01 0.01 0.01 0.01 0.01 0')
  aligner_path.add_processor(aligner)

  # Finished with path for aligner
  # Repeat this 3x
  calpaths.append(aligner_path)
  calpaths.append(aligner_path)
  calpaths.append(aligner_path)

  # Create path for some track based dqm using current calibrations
  dqm_path = Env.create_path('dqm_path')
  dqm_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': "tmp.slcio" })
  dqm_path.add_processor(geo)
  dqm_path = add_hitmakers(dqm_path)
  dqm_path.add_processor(trackfinder_tightcut)

  teldqm = tbsw.Processor(name="TelescopeDQM", proctype="TrackFitDQM")
  teldqm.param("RootFileName","TelescopeDQM_basic.root")
  dqm_path.add_processor(teldqm)

  # Finished with path for teldqm
  calpaths.append(dqm_path)
           
  return calpaths

def create_reco_path(Env, rawfile, gearfile, maxrecoNo):
  """
  Returns a list of tbsw path objects for reconstruction of a test beam run
  """

  reco = Env.create_path('reco')
  
  reco.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxrecoNo})

  rawinput = tbsw.Processor(name="RawInputProcessor",proctype="EudaqInputProcessor")
  rawinput.param('FileNames', rawfile)
  reco.add_processor(rawinput)

  geo = tbsw.Processor(name="Geo",proctype="Geometry")
  geo.param("AlignmentDBFilePath", "localDB/alignmentDB.root")
  geo.param("ApplyAlignment", "true")
  geo.param("OverrideAlignment", "true")
  reco.add_processor(geo)

  if usePseudoPixel: # not set up to use pseudo pixel
    reco = add_unpackersPseudo(reco)
    diasorter = tbsw.Processor(name="DiamondRawHitSorter", proctype="DiamondRawHitSorter")
    diasorter.param("PixelType", diamond_pixeltype)
    reco.add_processor(diasorter)
  else: 
    reco = add_unpackers(reco)
   
  reco = add_clusterizers(reco)

  #trackfinder RecoTF
  trackfinder = tbsw.Processor(name="RecoTF", proctype="FastTracker")
  trackfinder.param("InputHitCollectionNameVec","hit_m26 hit_fei4")
  trackfinder.param("ExcludeDetector", "3")
  trackfinder.param("MaxTrackChi2", "100")
  trackfinder.param("MaximumGap", "1")
  trackfinder.param("MinimumHits","6")
  trackfinder.param("OutlierChi2Cut", "20")
  trackfinder.param("ParticleCharge","-1")
  trackfinder.param("ParticleMass","0.000511")
  trackfinder.param("ParticleMomentum", beamEnergy)
  trackfinder.param("SingleHitSeeding", singleHitSeeding)
  trackfinder.param("MaxResidualU","0.4")
  trackfinder.param("MaxResidualV","0.4")
  reco.add_processor(trackfinder)

  diamond_analyser = tbsw.Processor(name="DiamondAnalyzer", proctype="PixelDUTAnalyzer")
  diamond_analyser.param("HitCollection", "hit_dia")
  diamond_analyser.param("DigitCollection","zsdata_dia")
  diamond_analyser.param("TrackCollection", "tracks")
  diamond_analyser.param("DUTPlane", "3")
  diamond_analyser.param("MaxResidualU", "0.4")
  diamond_analyser.param("MaxResidualV", "0.4")
  diamond_analyser.param("ReferencePlane", "4")
  diamond_analyser.param("RootFileName", "Histos-DIA.root")
  reco.add_processor(diamond_analyser)

  fei4_analyser = tbsw.Processor(name="FEI4Analyzer", proctype="PixelDUTAnalyzer")
  fei4_analyser.param("HitCollection", "hit_fei4")
  fei4_analyser.param("DigitCollection","zsdata_fei4")
  fei4_analyser.param("TrackCollection", "tracks") 
  fei4_analyser.param("DUTPlane", "4")
  fei4_analyser.param("MaxResidualU", "0.2")
  fei4_analyser.param("MaxResidualV", "0.2")
  fei4_analyser.param("ReferencePlane", "-1")
  fei4_analyser.param("RootFileName", "Histos-FEI.root")
  reco.add_processor(fei4_analyser)

  return [reco]


def add_unpackers(path):
  """
  Adds unpackers to the path
  """

  m26unpacker = tbsw.Processor(name="M26Unpacker",proctype="NIUnpacker")
  m26unpacker.param("InputCollectionName", "NI")
  m26unpacker.param("OutputCollectionName","zsdata_m26")
  path.add_processor(m26unpacker)

  fei4unpacker = tbsw.Processor(name="FEI4Unpacker", proctype="USBPixUnpacker")
  fei4unpacker.param("OutputCollectionName", "zsdata_fei4")
  fei4unpacker.param("FilterChannel", filterChannelRef)
  path.add_processor(fei4unpacker)

  diaunpacker = tbsw.Processor(name="DiamondUnpacker", proctype="USBPixUnpacker")
  diaunpacker.param('OutputCollectionName','zsdata_dia')
  diaunpacker.param("FilterChannel", filterChannelDia)
  path.add_processor(diaunpacker)

  return path

def add_unpackersPseudo(path):
  """
  Adds unpackers to the path using DB pixel alignment
  """

  m26unpacker = tbsw.Processor(name="M26Unpacker",proctype="NIUnpacker")
  m26unpacker.param("InputCollectionName", "NI")
  m26unpacker.param("OutputCollectionName","zsdata_m26")
  path.add_processor(m26unpacker)

  fei4unpacker = tbsw.Processor(name="FEI4Unpacker", proctype="USBPixUnpacker")
  fei4unpacker.param("OutputCollectionName", "zsdata_fei4")
  fei4unpacker.param("FilterChannel", filterChannelRef)
  path.add_processor(fei4unpacker)

  diaunpacker = tbsw.Processor(name="DiamondUnpacker", proctype="USBPixUnpacker")
  diaunpacker.param('OutputCollectionName','zsdata_raw_dia')
  diaunpacker.param("FilterChannel", filterChannelDia)
  path.add_processor(diaunpacker)

  return path

def add_clusterizers(path):
  """
  Adds clusterizers to the path
  """

  m26clust = tbsw.Processor(name="M26Clusterizer",proctype="PixelClusterizer")
  m26clust.param("NoiseDBFileName","localDB/NoiseDB-M26.root")
  m26clust.param("SparseDataCollectionName","zsdata_m26")
  m26clust.param("ClusterCollectionName","zscluster_m26")
  m26clust.param("SparseClusterCut",0)
  m26clust.param("SparseSeedCut", 0)
  m26clust.param("SparseZSCut", 0)
  path.add_processor(m26clust)

  fei4clust = tbsw.Processor(name="FEI4Clusterizer",proctype="PixelClusterizer")
  fei4clust.param("NoiseDBFileName","localDB/NoiseDB-FEI4.root")
  fei4clust.param("SparseDataCollectionName","zsdata_fei4")
  fei4clust.param("ClusterCollectionName","zscluster_fei4")
  fei4clust.param("SparseClusterCut",0)
  fei4clust.param("SparseSeedCut", 0)
  fei4clust.param("SparseZSCut", 0)
  path.add_processor(fei4clust)

  diaclust = tbsw.Processor(name="DiamondClusterizer",proctype="PixelClusterizer")
  diaclust.param("NoiseDBFileName","localDB/NoiseDB-DIA.root")
  diaclust.param("SparseDataCollectionName","zsdata_dia")
  diaclust.param("ClusterCollectionName","zscluster_dia")
  diaclust.param("SparseClusterCut",0)
  diaclust.param("SparseSeedCut", 0)
  diaclust.param("SparseZSCut", 0)
  path.add_processor(diaclust)

  return path


def add_hitmakers(path):
  """
  Add CoG hitmakers to the path
  """

  m26hitmaker = tbsw.Processor(name="M26CogHitMaker",proctype="CogHitMaker")
  m26hitmaker.param("ClusterCollection","zscluster_m26")
  m26hitmaker.param("HitCollectionName","hit_m26")
  #m26hitmaker.param("SigmaUCorrections", "0.698 0.31 0.315")
  #m26hitmaker.param("SigmaVCorrections", "0.698 0.31 0.315")
  path.add_processor(m26hitmaker)

  fei4hitmaker = tbsw.Processor(name="FEI4CogHitMaker",proctype="CogHitMaker")
  fei4hitmaker.param("ClusterCollection","zscluster_fei4")
  fei4hitmaker.param("HitCollectionName","hit_fei4")
  #fei4hitmaker.param("SigmaUCorrections", "1.0 0.5 0.3")
  #fei4hitmaker.param("SigmaVCorrections", "1.0 0.5 0.3")
  path.add_processor(fei4hitmaker)

  diahitmaker = tbsw.Processor(name="DiamondCogHitMaker",proctype="CogHitMaker")
  diahitmaker.param("ClusterCollection","zscluster_dia")
  diahitmaker.param("HitCollectionName","hit_dia")
  #diahitmaker.param("SigmaUCorrections", "0.8 0.3 0.3")
  #diahitmaker.param("SigmaVCorrections", "0.8 0.3 0.3")
  path.add_processor(diahitmaker)

  return path

def add_hitmakersDB(path):
  """
  Add cluster shape hitmakers to the path (requiring clusterDBs)
  """

  m26goehitmaker = tbsw.Processor(name="M26GoeHitMaker",proctype="GoeHitMaker")
  m26goehitmaker.param("ClusterCollection","zscluster_m26")
  m26goehitmaker.param("HitCollectionName","goehit_m26")
  m26goehitmaker.param("ClusterDBFileName","localDB/clusterDB-M26.root")
  path.add_processor(m26goehitmaker)

  fei4goehitmaker = tbsw.Processor(name="FEI4GoeHitMaker",proctype="GoeHitMaker")
  fei4goehitmaker.param("ClusterCollection","zscluster_fei4")
  fei4goehitmaker.param("HitCollectionName","goehit_fei4")
  fei4goehitmaker.param("ClusterDBFileName","localDB/clusterDB-FEI4.root")
  path.add_processor(fei4goehitmaker)

  diagoehitmaker = tbsw.Processor(name="DiamondGoeHitMaker",proctype="GoeHitMaker")
  diagoehitmaker.param("ClusterCollection","zscluster_dia")
  diagoehitmaker.param("HitCollectionName","goehit_dia")
  diagoehitmaker.param("ClusterDBFileName","localDB/clusterDB-DIA.root")
  diagoehitmaker.param("UseCenterOfGravityFallback","true")
  #diagoehitmaker.param("SigmaUCorrections", "0.8 0.3 0.3")
  #diagoehitmaker.param("SigmaVCorrections", "0.8 0.3 0.3")
  path.add_processor(diagoehitmaker)

  return path

# hitmakers only for mimosa and ref
def add_hitmakers_nodia(path):
  """
  Add CoG hitmakers to the path
  """

  m26hitmaker = tbsw.Processor(name="M26CogHitMaker",proctype="CogHitMaker")
  m26hitmaker.param("ClusterCollection","zscluster_m26")
  m26hitmaker.param("HitCollectionName","hit_m26")
  #m26hitmaker.param("SigmaUCorrections", "0.698 0.31 0.315")
  #m26hitmaker.param("SigmaVCorrections", "0.698 0.31 0.315")
  path.add_processor(m26hitmaker)

  fei4hitmaker = tbsw.Processor(name="FEI4CogHitMaker",proctype="CogHitMaker")
  fei4hitmaker.param("ClusterCollection","zscluster_fei4")
  fei4hitmaker.param("HitCollectionName","hit_fei4")
  #fei4hitmaker.param("SigmaUCorrections", "1.0 0.5 0.3")
  #fei4hitmaker.param("SigmaVCorrections", "1.0 0.5 0.3")
  path.add_processor(fei4hitmaker)

  return path

def add_hitmakersDB_nodia(path):
  """
  Add cluster shape hitmakers to the path (requiring clusterDBs)
  """

  m26goehitmaker = tbsw.Processor(name="M26GoeHitMaker",proctype="GoeHitMaker")
  m26goehitmaker.param("ClusterCollection","zscluster_m26")
  m26goehitmaker.param("HitCollectionName","goehit_m26")
  m26goehitmaker.param("ClusterDBFileName","localDB/clusterDB-M26.root")
  path.add_processor(m26goehitmaker)

  fei4goehitmaker = tbsw.Processor(name="FEI4GoeHitMaker",proctype="GoeHitMaker")
  fei4goehitmaker.param("ClusterCollection","zscluster_fei4")
  fei4goehitmaker.param("HitCollectionName","geohit_fei4")
  fei4goehitmaker.param("ClusterDBFileName","localDB/clusterDB-FEI4.root")
  path.add_processor(fei4goehitmaker)

  return path


def calibrate(params):

  rawfile, gearfile = params

  # Tag for calibration data
  caltag = os.path.splitext(os.path.basename(rawfile))[0] 
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag
  # containing all calibration data.
  CalObj = tbsw.Calibration(steerfiles=steerfiles, name=caltag + '-cal')

  # Create list of calibration paths
  calpath = create_calibration_path(CalObj, rawfile, gearfile)

  # Run the calibration steps
  CalObj.calibrate(paths=calpath,ifile=rawfile,caltag=caltag)


def reconstruct(params):

  calibfile, rawfile, gearfile, maxrecoNo = params

  # Tag for calibration data
  caltag = os.path.splitext(os.path.basename(calibfile))[0] 

  # Reconsruct the rawfile using caltag. Resulting root files are
  # written to folder root-files/
  recotag = os.path.splitext(os.path.basename(rawfile))[0] 
  
  RecObj = tbsw.Reconstruction(steerfiles=steerfiles, name=recotag + '-reco' )

  # Create reconstuction path
  recopath = create_reco_path(RecObj, rawfile, gearfile, maxrecoNo)

  # Run the reconstuction
  RecObj.reconstruct(paths=recopath,ifile=rawfile,caltag=caltag)


if __name__ == '__main__':
  count = multiprocessing.cpu_count()
  pool = multiprocessing.Pool(processes=count)

  try:
    import psutil
    parent = psutil.Process()
    parent.nice(18)
    for child in parent.children():
      child.nice(18)
    print("All %d processes set to minimum priority"%(count))
  except:
    print("Could not change process priority. Please install psutil")
    pass
 
  # calibrate with the first rawfile of the runlist and use this calibration for the other runs to save time, because the setup did not change between runs
  params = [ (calibrun, '{}{}.xml'.format(gearBaseName, gearBatch))]
  pool.map(calibrate, params)

  # reconstruct runs based on the one calibration per pixel type
  #params = [(calibrun, rawfile, '{}{}.xml'.format(gearBaseName, gearBatch), maxrecoNo) for rawfile, maxrecoNo in zip(runlist, nEventslist)]
  params = [(calibrun, runlist, '{}{}.xml'.format(gearBaseName, gearBatch), nEventslist)]
  pool.map(reconstruct, params)
