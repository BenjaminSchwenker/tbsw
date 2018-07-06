"""
This is an example script to demonstrate how TBSW can be used to create an X0 image from a simulated
test beam experiment.

The script below simulates a test beam experiment where charged tracks from a monoenergetic beam 
cross a a misaligned pixel telescope containing six Mimosa 26 detector planes. Two data sets are 
simulated. A first 'air' run is done with no additional scatterer between the telescope arms. The 
'air' run is used to calibrate the telescope. In a second 'aluminium' run, a aluminium plate with 
a well known thickness profile is inserted in between the telescope arms. This second run is used 
to compute a X0 image from the reconstructed scattering angles. The known comparison between the 
reconstructed X0 image and the a priori known image is used to calibrate the beam energy and the 
angular resolution of the telescope. This second step completes the calibration of the telescope
for X0 imaging. 

Author: Ulf Stolzenberg <ulf.stolzenberg@phys.uni-goettingen.de>  
"""

from tbsw import *
import multiprocessing

import tbsw.x0imaging.X0Calibration

# Path to steering files 
# Steeringfiles are xml files and define details of the simulation like how many events are produced
# or how M26 sensors are digitized. XML parameters can be adjusted using any test editor

# Steerfiles for the telescope calibration
steerfiles_cali = 'steering-files/x0-tboct16-2GeV-air/'

# Steerfiles for the angle reconstruction (can be the same directory as telescope calibration steerfiles)
steerfiles_reco = 'steering-files/x0-tboct16-2GeV-air/'

# Steerfiles for the x0calibration/x0imaging (can be the same directory as telescope calibration steerfiles)
steerfiles_x0 = 'steering-files/x0-tboct16-2GeV-air/'

# Nominal Beam energy
beamenergy=2.0

# cal tags
# telescope calibration cal tag (typically named after telescope setup, beam energy etc.)
caltag='40mm-spacing-2GeV'

# x0 calibration cal tag
x0tag='05mmalu'

# Name of the gearfile, which describes the telescope setup 
gearfile = 'gear.xml'

# Alignment DB file name
alignmentdb_filename='alignmentDB.root'

# Number of iterations during target alignment
# Set to 0 or negative integer to disable target alignment
targetalignment_iterations=0

# File names and lists of filenames for the different steps 

# global path to raw files
rawfile_path='/work1/rawdata/tboct16/'

# raw file used during telescope calibration (best use data with scattering target)
# The calibration has to be done for every telescope setup, beam energy and m26 threshold settings
cali_run='run006973.raw'
rawfile_cali = rawfile_path + cali_run

# raw file used for target alignment (only useful with a thick (X/X0 > 5 %) scattering target)
TA_run='run006958.raw'
rawfile_TA = rawfile_path + TA_run

# List of runs, which are used as input for the scattering angle reconstruction
# The angle reconstruction step is essential and every run, that will be used later during the x0 calibration or x0 imaging steps, must be listed
RunList_reco = [
		    'run006958.raw',
		    'run006959.raw',
		    'run006960.raw',
		    'run006973.raw',
          ]

RawfileList_reco = [rawfile_path+x for x in RunList_reco]

# List of runs, which are input for the x0 calibration
# Typically runs with various different materials and thicknesses have to be used to achieve a sensible calibration
# The different measurement regions and other options have to be set in the x0.cfg file in the steer files directory
RunList_x0cali = [
		    'run006958.raw',
		    'run006973.raw',
          ]

RawfileList_x0cali = [rawfile_path+x for x in RunList_x0cali]

# List of runs, which are input for the first x0 image
# Use only runs, with exactly the same target material and positioning
RunList_x0image = [
		    'run006958.raw',
		    'run006959.raw',
		    'run006960.raw',
          ]

RawfileList_x0image = [rawfile_path+x for x in RunList_x0image]

# Number of events ...
# for telescope calibration
nevents_cali = 700000

# for target alignment
nevents_TA = 1000000

# for angle reconstruction (-1 use all available events)
nevents_reco = -1

# Processor settings and sequence during telescope calibration
def create_calibration_path(Env, rawfile, gearfile):
  """
  Returns a list of tbsw path objects to calibrate the tracking telescope
  """
  
  hotpixelkiller = Env.create_path('hotpixelkiller')
  hotpixelkiller.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000})  
  hotpixelkiller.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  hotpixelkiller.add_processor(name="M26Unpacker")
  hotpixelkiller.add_processor(name="M26HotPixelKiller")
  
  clusterizer = Env.create_path('clusterizer')
  clusterizer.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali}) 
  clusterizer.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  clusterizer.add_processor(name="M26Unpacker")
  clusterizer.add_processor(name="M26Clusterizer")
  clusterizer.add_processor(name="LCIOOutput")

  correlator = Env.create_path('correlator')
  correlator.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" }) 
  correlator.add_processor(name="M26CogHitMaker")
  correlator.add_processor(name="RawDQM") 
  correlator.add_processor(name="TelCorrelator")
  
  kalman_aligner_1 = Env.create_path('kalman_aligner_1')
  kalman_aligner_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_1.add_processor(name="M26CogHitMaker")
  kalman_aligner_1.add_processor(name="AlignTF_LC")
  kalman_aligner_1.add_processor(name="PreAligner")
  
  kalman_aligner_2 = Env.create_path('kalman_aligner_2')
  kalman_aligner_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_2.add_processor(name="M26CogHitMaker")
  kalman_aligner_2.add_processor(name="AlignTF_TC")
  kalman_aligner_2.add_processor(name="TelAligner")
  
  telescope_dqm = Env.create_path('telescope_dqm')
  telescope_dqm.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm.add_processor(name="M26CogHitMaker")
  telescope_dqm.add_processor(name="AlignTF_TC")
  telescope_dqm.add_processor(name="TelescopeDQM")
  
  cluster_calibration_1 = Env.create_path('cluster_calibration_1')
  cluster_calibration_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" }) 
  cluster_calibration_1.add_processor(name="M26CogHitMaker") 
  cluster_calibration_1.add_processor(name="AlignTF_TC")
  cluster_calibration_1.add_processor(name="M26ClusterCalibrator")
  
  kalman_aligner_3 = Env.create_path('kalman_aligner_3')
  kalman_aligner_3.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_3.add_processor(name="M26GoeHitMaker")
  kalman_aligner_3.add_processor(name="AlignTF_TC")
  kalman_aligner_3.add_processor(name="TelAligner")
  
  cluster_calibration_2 = Env.create_path('cluster_calibration_2')
  cluster_calibration_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" })  
  cluster_calibration_2.add_processor(name="M26GoeHitMaker") 
  cluster_calibration_2.add_processor(name="AlignTF_TC")
  cluster_calibration_2.add_processor(name="M26ClusterCalibrator")

  correlator2 = Env.create_path('correlator2')
  correlator2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })
  correlator2.add_processor(name="M26GoeHitMaker")  
  correlator2.add_processor(name="RawDQM", params={'RootFileName': 'RawDQM2.root'})
  correlator2.add_processor(name="TelCorrelator", params={'OutputRootFileName': 'XCorrelator2.root'})
  
  kalman_aligner_4 = Env.create_path('kalman_aligner_4')
  kalman_aligner_4.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" }) 
  kalman_aligner_4.add_processor(name="M26GoeHitMaker")   
  kalman_aligner_4.add_processor(name="AlignTF_LC")
  kalman_aligner_4.add_processor(name="PreAligner")
  
  kalman_aligner_5 = Env.create_path('kalman_aligner_5')
  kalman_aligner_5.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" }) 
  kalman_aligner_5.add_processor(name="M26GoeHitMaker")  
  kalman_aligner_5.add_processor(name="AlignTF_TC")
  kalman_aligner_5.add_processor(name="TelAligner") 
  
  telescope_dqm2 = Env.create_path('telescope_dqm2')
  telescope_dqm2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" }) 
  telescope_dqm2.add_processor(name="M26GoeHitMaker")  
  telescope_dqm2.add_processor(name="AlignTF_TC")
  telescope_dqm2.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM2.root'})
  
  
  # create sequence of calibration paths 
  calpath= [ hotpixelkiller ,
			 clusterizer, 
             correlator, 
             kalman_aligner_1, 
             kalman_aligner_2, 
             kalman_aligner_2, 
             kalman_aligner_2, 
             telescope_dqm, 
             cluster_calibration_1,
             kalman_aligner_3, 
             kalman_aligner_3, 
             kalman_aligner_3,  
             cluster_calibration_2, 
             cluster_calibration_2, 
             cluster_calibration_2, 
             cluster_calibration_2, 
             cluster_calibration_2, 
             cluster_calibration_2, 
             cluster_calibration_2, 
             cluster_calibration_2, 
             correlator2, 
             kalman_aligner_4, 
             kalman_aligner_5, 
             kalman_aligner_5, 
             kalman_aligner_5, 
             telescope_dqm2, 
           ]
  
  return calpath

# Processor settings and sequence during angle reconstruction
def create_reco_path(Env, rawfile, gearfile, numberofevents):
  """
  Returns a list of tbsw path objects to reconstruct a test beam run 
  """
  
  reco = Env.create_path('reco')
  reco.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : numberofevents}) 
  reco.add_processor(name="RawInputProcessor", params={'FileName': rawfile})
  reco.add_processor(name="M26Unpacker")
  reco.add_processor(name="M26Clusterizer")
  reco.add_processor(name="M26GoeHitMaker")
  reco.add_processor(name="DownstreamFinder")
  reco.add_processor(name="UpstreamFinder")
  reco.add_processor(name="X0Imager")
    
  return [ reco ]


# Perform the telescope calibration
def calibrate(params):
  """
  Calibrates an misaligned tracking telescope from run data. 
  Creates a folder localDB/caltag in workspace containing 
  calibration results. 
  Creates a folder tmp-runs/name-sim/ and populates it with 
  Marlin steering and logfiles.  
  """ 
  
  rawfile, steerfiles, gearfile, caltag = params
  print(gearfile)
  
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = Calibration(steerfiles=steerfiles, name=caltag + '-cal') 
  CalObj.set_beam_momentum(beamenergy)
  
  # Create list of calibration steps 
  calpath = create_calibration_path(CalObj, rawfile, gearfile)
  
  # Run the calibration steps 
  CalObj.calibrate(path=calpath,ifile=rawfile,caltag=caltag)



# Perform the angle reconstruction of a single run
def reconstruct(params):

  rawfile, steerfiles, gearfile, caltag = params

  # Set cal tag that includes run name
  name = os.path.splitext(os.path.basename(rawfile))[0] + '-' + caltag
  
  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=name + '-reco' )
  RecObj.set_beam_momentum(beamenergy)

  # Create reconstuction path
  recopath = create_reco_path(RecObj, rawfile, gearfile, nevents_reco)  

  # Use caltag of the last target alignment iteration
  iteration_string='-target-alignment-it'+str(targetalignment_iterations-1)
  localcaltag=caltag+iteration_string

  if targetalignment_iterations < 1:
    localcaltag=caltag

  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile,caltag=localcaltag) 

# Perform target alignment
def targetalignment(params):
  """
  Starts the scattering angle reconstruction and vertex fit on the central target
  plane. Afterwards the mean vertex z position is set as the new target z position in
  the aligment DB file and the calibration tag is exported
    :@params:       consists of rawfile, steerfiles, gearfile, caltag, iteration
    :@rawfile:      Input file for the reconstruction 
    :@BE            Nominal beam energy of the run
    :@nevents       Number of events
    :@steerfiles:   Directory with the steering files for the reconstruction
    :@gearfile:     Name of the gear file
    :@caltag:       calibration tag for the reconstruction
    :@iteration:    Target alignment iteration counter  
    :author: benjamin.schwenker@phys.uni-goettinge.de  
  """ 

  rawfile, steerfiles, calibrationtag, iteration = params

  if rawfile == None:
    return None

  if iteration == None:
    return None

  prev_iteration_string='-target-alignment-it'+str(iteration-1)
  curr_iteration_string='-target-alignment-it'+str(iteration)

  localcaltag=caltag+prev_iteration_string

  if iteration == 0:
    localcaltag=calibrationtag

  newcaltag=calibrationtag+curr_iteration_string

  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=newcaltag )

  # Set Beam energy
  RecObj.set_beam_momentum(beamenergy)
  
  recopath = create_reco_path(RecObj, rawfile, gearfile, nevents_TA)

  # Run the reconstuction  
  RecObj.reconstruct(path=reco,ifile=rawfile,caltag=localcaltag) 

  # Read the vertex position and save it in the alignmentDB
  dbname=RecObj.create_dbfilename("alignmentDB.root")
  treename=RecObj.get_rootfilename('X0')
  save_targetpos(treename,dbname)
  RecObj.export_caltag(newcaltag)

# Perform x0 calibration
def xx0calibration(params):

  x0caltag, RunList, steerfiles, calibrationtag, delete = params

  # Base filename of the X0 root file
  basefilename='X0cal-merge-'+x0caltag

  # Total path of X0 root file
  filename='root-files/'+basefilename+'.root'

  # Merge the root trees in the root files directory
  tbsw.x0imaging.X0Calibration.merge_rootfile(filename=filename,RunList=RunList,caltag=caltag)

  # Generate a uncalibrated X/X0 image
  tbsw.x0imaging.X0Calibration.x0imaging(filename=filename,caltag='',deletetag=delete,steerfiles=steerfiles,nametag='Uncalibrated')

  # Path to uncalibrated X0 image file
  imagefilename='/root-files/'+basefilename+'-UncalibratedX0image.root'

  # Do a calibration of the angle resolution
  tbsw.x0imaging.X0Calibration.x0calibration(filename=filename,imagefilename=imagefilename,caltag=x0caltag,steerfiles=steerfiles)

# Generate x0 image
def xx0image(params):

  x0caltag, RunList, steerfiles, calibrationtag, delete, listnametag = params

  # Base filename of the X0 root file
  basefilename='X0-'+listnametag+'-merge-'+x0caltag

  # Total path of X0 root file
  filename='root-files/'+basefilename+'.root'

  # Merge the root trees in the root files directory
  tbsw.x0imaging.X0Calibration.merge_rootfile(filename=filename,RunList=RunList,caltag=caltag)

  # Do a calibration of the angle resolution
  tbsw.x0imaging.X0Calibration.x0imaging(filename=filename,caltag=x0caltag,deletetag=deletetag,steerfiles=steerfiles_reco,nametag='Calibrated')

  
if __name__ == '__main__':


  # Calibrate the telescope 
  params_cali = ( rawfile_cali, steerfiles_cali, gearfile, caltag)
  print(params_cali)
  #calibrate( params_cali )


  # Target alignment
  for it in range(0,targetalignment_iterations):
    params_TA = (rawfile_TA, steerfiles_reco, caltag, it)
    print "The parameters for the target alignment are: " 
    print params_TA

    #targetalignment(params_TA)


  # Angle reconstruction
  params_reco=[(x, steerfiles_reco, gearfile, caltag) for x in RawfileList_reco]
  print "The parameters for the reconstruction are: " 
  print params_reco

  count = multiprocessing.cpu_count()
  pool = multiprocessing.Pool(processes=count)
  #pool.map(reconstruct, params_reco)


  # start x0 calibration
  deletetag='1'
  params_x0cali = ( x0tag, RawfileList_x0cali, steerfiles_x0, caltag, deletetag)
  xx0calibration(params_x0cali)


  # Generate a calibrated X/X0 image
  nametag='image1'
  params_x0image = ( x0tag, RawfileList_x0image, steerfiles_x0, caltag, deletetag, nametag)
  xx0image(params_x0image)

