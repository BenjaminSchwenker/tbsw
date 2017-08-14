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

# Path to steering files 
# Steeringfiles are xml files and define details of the simulation like how many events are produced
# or how M26 sensors are digitized. XML parameters can be adjusted using any test editor
steerfiles_cali = 'steering-files/x0-tboct16-2GeV-air/'
steerfiles_reco = 'steering-files/x0-tboct16-2GeV-air/'
# Gearfile for runs 
gearfile = 'gear.xml'
#cal tag
caltag='air'
# File name for raw data 
rawfile_cali = '/work1/rawdata/DESY_Oktober16/2GeV_air/run006973.raw'
# Nominal Beam energy
beamenergy=2.0

RunList_reco = [
		    '/work1/rawdata/DESY_Oktober16/2GeV_05mmalu/run006958.raw',
		    '/work1/rawdata/DESY_Oktober16/2GeV_05mmalu/run006959.raw',
		    '/work1/rawdata/DESY_Oktober16/2GeV_05mmalu/run006960.raw',
		    '/work1/rawdata/DESY_Oktober16/2GeV_05mmalu/run006961.raw',
          ]
# Number of events to simulate 
nevents_cali = 700000
nevents_reco = -1


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
  kalman_aligner_1.add_processor(name="PreAligner", params={'RootFileName': 'KalmanAlign-iteration-1.root'})
  
  kalman_aligner_2 = Env.create_path('kalman_aligner_2')
  kalman_aligner_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_2.add_processor(name="M26CogHitMaker")
  kalman_aligner_2.add_processor(name="AlignTF_TC")
  kalman_aligner_2.add_processor(name="TelAligner", params={'RootFileName': 'KalmanAlign-iteration-2.root'})
  
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
  kalman_aligner_3.add_processor(name="TelAligner", params={'RootFileName': 'KalmanAlign-iteration-3.root'})
  
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
  kalman_aligner_4.add_processor(name="PreAligner", params={'RootFileName': 'KalmanAlign-iteration-4.root'})
  
  kalman_aligner_5 = Env.create_path('kalman_aligner_5')
  kalman_aligner_5.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" }) 
  kalman_aligner_5.add_processor(name="M26GoeHitMaker")  
  kalman_aligner_5.add_processor(name="AlignTF_TC")
  kalman_aligner_5.add_processor(name="TelAligner", params={'RootFileName': 'KalmanAlign-final.root'}) 
  
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


def create_reco_path(Env, rawfile, gearfile):
  """
  Returns a list of tbsw path objects to reconstruct a test beam run 
  """
  
  reco = Env.create_path('reco')
  reco.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_reco}) 
  reco.add_processor(name="RawInputProcessor", params={'FileName': rawfile})
  reco.add_processor(name="M26Unpacker")
  reco.add_processor(name="M26Clusterizer")
  reco.add_processor(name="M26GoeHitMaker")
  reco.add_processor(name="DownstreamFinder")
  reco.add_processor(name="UpstreamFinder")
  reco.add_processor(name="X0Imager")
    
  return [ reco ]


def reconstruct(params):

  
  rawfile, steerfiles, gearfile, caltag = params

  # Set cal tag that includes run name
  name = os.path.splitext(os.path.basename(rawfile))[0] + '-' + caltag
  
  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=name + '-reco' )
  RecObj.set_beam_momentum(beamenergy)

  # Create reconstuction path
  recopath = create_reco_path(RecObj, rawfile, gearfile)  

  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile,caltag=caltag) 

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
  


# Function which starts the imaging script
def x0imaging(filename,caltag,deletetag):

  flags='./tbsw/x0imaging/GenerateImage.py -i '+filename+' -f '+imagecfgfilename+' -c '+caltag+' -d '+`deletetag`
  print('Starting X0 imaging')
  print(flags)
  subprocess.call(flags, shell=True)

  return None

# Function which starts the x0 calibration script
def x0calibration(filename,imagefilename,caltag):

  flags='./tbsw/x0imaging/X0Calibration.py -i '+filename+' -f '+calibrationcfgfilename+' -m '+imagefilename+' -c '+caltag
  print('Starting X0 calibration')
  print(flags)
  subprocess.call(flags, shell=True)

  return None

# Function which merges the result root files
def merge_rootfile(filename,RunList_reco):

  flags='hadd '+filename+' '
  for run in RunList_reco:
    name=os.path.splitext(os.path.basename(run))[0]
    flags=flags+'root-files/X0-'+name+'-air-reco.root '

  if os.path.isfile(filename):
    os.remove(filename)
  subprocess.call(flags, shell=True)

  return None

  
if __name__ == '__main__':

  # Get current directory
  cwdir = os.getcwd()

  params_cali = ( rawfile_cali, steerfiles_cali, gearfile, caltag)

  print(params_cali)

  # Calibrate the telescope 
  calibrate( params_cali )

  params_reco=[(x, steerfiles_reco, gearfile, caltag) for x in RunList_reco]
  print "The parameters for the reconstruction are: " 
  print params_reco

  count = multiprocessing.cpu_count()
  pool = multiprocessing.Pool(processes=count)
  pool.map(reconstruct, params_reco)

  imagecfgfilename='steering-files/x0-tb/image.cfg'
  calibrationcfgfilename='steering-files/x0-tb/x0calibration.cfg'
  deletetag=1

  # Base filename of the X0 root file
  basefilename='X0-merge-05mmalu'

  # Total path of X0 root file
  filename='root-files/'+basefilename+'.root'

  # Total path if the different kinds of X0 image files
  imagefile=cwdir+'/root-files/'+basefilename+'-X0image.root'
  uncalibratedimagefile=cwdir+'/root-files/'+basefilename+'-Uncalibrated-X0image.root'
  calibratedimagefile=cwdir+'/root-files/'+basefilename+'-Calibrated-X0image.root'

  # Merge the root trees in the root files directory
  merge_rootfile(filename,RunList_reco)

  # Total path if the different kinds of X0 image files
  imagefile=cwdir+'/root-files/'+basefilename+'-X0image.root'

  # Generate a uncalibrated X/X0 image
  x0imaging(filename,'',deletetag)

  # Rename the image file
  if os.path.isfile(uncalibratedimagefile):
     os.remove(uncalibratedimagefile) 
  if os.path.isfile(imagefile):
     shutil.move(imagefile,uncalibratedimagefile)

  # Rename the directory in tmp-runs, in which the image was generated
  if os.path.isdir(uncaltmpdir):
     shutil.rmtree(uncaltmpdir) 
  if os.path.isdir(tmpdir):
     shutil.copytree(tmpdir, uncaltmpdir)

  # Do a calibration of the angle resolution
  x0calibration(filename,imagefile,air_caltag)

  # Generate a calibrated X/X0 image
  x0imaging(filename,air_caltag,deletetag)

  # Rename the image file
  if os.path.isfile(calibratedimagefile):
     os.remove(calibratedimagefile) 
  if os.path.isfile(imagefile):
     shutil.move(imagefile,calibratedimagefile)

  # Rename the directory in tmp-runs, in which the image was generated
  if os.path.isdir(caltmpdir):
     shutil.rmtree(caltmpdir) 
  if os.path.isdir(tmpdir):
     shutil.copytree(tmpdir, caltmpdir)

