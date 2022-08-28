"""
Some helper code to define processor paths for X/X0 studies

:author: ulf.stolzenberg@phys.uni-goettinge.de  
"""

 

from tbsw.tbsw import Simulation, Calibration, Reconstruction
import tbsw.X0Calibration as X0Calibration
import tbsw.GenerateImage as GenerateImage 
import tbsw.path_utils as path_utils 
import tbsw.DQMplots as DQMplots 
import tbsw.gear as gear

import os


def simulate(rawfile_air=None, rawfile_alu_list=None, steerfiles=None, gearfile="gear.xml", nevents_air=1000000, nevents_alu=6000000, nprocesses=2, mean_list=[0.0,0.0,0.0,0.0,0.0,0.0], sigma_list=[0.5,0.5,1.0,0.3,0.3,1.5], sensorexception_list=[5,0,11],  modeexception_list=[''], beamenergy=2.0): 
  """
  Generates air and aluminium  slcio files from a simulated test beam experiment
  Creates a folder in tmp-runs/ and populates it with 
  Marlin steering and logfiles. The slcio files are stored in the base directory of 
  the workspace 
  :@rawfile_air:       			Name of the air slcio file to be used   
  :@rawfile_alu_list:       	Names of the aluminium slcio file to be used   
  :@steerfiles     				Name of the directory with the geometry information
  :@gearfile       				Name of the gear file in the steering files directory
  :@nevents_air      			Number of air events
  :@nevents_alu      			Number of alu events
  :@nprocesses      			Number individual aluminium runs
  :@mean_list		  			List with gaussian mean values of the 6 alignment parameters x,y,z,alpha,beta,gamma
  :@mean_list		  			List with gaussian sigma values of the 6 alignment parameters x,y,z,alpha,beta,gamma,
 								the alignment parameter of individual sensor planes is determined from the corresponding gaussian
  :@sensorexception_list		List with sensors that are not artificially misaligned
  :@sensorexception_list		List with alignment parameters that are not artificially misaligned
  :@beamenergy     			Beam energy of the particle beam in GeV
  """ 

  if rawfile_air == None:
    print("Rawfile air name missing! Skip simulation.")
    return None

  if rawfile_alu_list == None:
    print("Rawfile alu list name missing! Skip simulation.")
    return None

  if steerfiles == None:
    print("Steerfiles name missing! Skip simulation.")
    return None
  
  # Create tmpdir to hold all steerfiles and log files 
  SimObj = Simulation(steerfiles=steerfiles, name='x0-sim' )
  
  # Get path to gearfile in simulation environment
  localgearfile = SimObj.get_filename(gearfile)
  
  # Misalign the telescope by applying random shifts to constants in local gearfile
  gear.randomize_telescope(gearfile=localgearfile, mean_list=mean_list, sigma_list=sigma_list, sensorexception_list=sensorexception_list, modeexception_list=modeexception_list)
  
  # Create a gearfile for the run with an Al plate 
  # Only the Al plate is added in between the telescope arms
  # The misaligned positions of the telescope planes remain
  # the same. 
  gearfile_alu = "gear_alu.xml"
  localgearfile_alu = SimObj.copy_file(gearfile, gearfile_alu)
  
  # Replace the air layer  by a Al plate with 0.5mm thickness 
  gear.set_parameter(gearfile=localgearfile_alu, sensorID=11, parametername='thickness', value=0.50)
  gear.set_parameter(gearfile=localgearfile_alu, sensorID=11, parametername='radLength', value=89.0)
  gear.set_parameter(gearfile=localgearfile_alu, sensorID=11, parametername='atomicNumber', value=13)
  gear.set_parameter(gearfile=localgearfile_alu, sensorID=11, parametername='atomicMass', value=27)  
  
  # List to populate with simulation pathes (=jobs)
  simpaths = []
   
  # Create path to simulate air run (just air between telescope arms) 
  simpath_air = path_utils.create_x0sim_path(SimObj, 'sim_air', rawfile_air, gearfile, nevents_air,  beamenergy)
  simpaths.append(simpath_air)
  
  # Create path to simulate alu run (Al plate between telescope arms)
  for nruns in range(0,nprocesses):
    simname='sim_alu_run{:d}'.format(nruns+1)
    simpath_alu = path_utils.create_x0sim_path(SimObj, simname, rawfile_alu_list[nruns], gearfile_alu, nevents_alu/nprocesses, beamenergy)
    simpaths.append(simpath_alu)   
  
  # Run simulation to create rawfile with simulated digits 
  SimObj.simulate(paths=simpaths,caltag='x0-sim-truthdb')  


# Perform the telescope calibration
def calibrate(rawfile=None, steerfiles=None, caltag="default", gearfile="gear.xml", nevents=50000, Use_clusterDB=True, beamenergy=2.0, mcdata=False, Use_LongTelescopeCali=True, UseOuterPlanesForClusterDB=False, csvdata=False):
  """
  Calibrates an misaligned tracking telescope from run data. 
  Creates a folder localDB/caltag in workspace containing 
  calibration results. 
  Creates a folder tmp-runs/name-sim/ and populates it with 
  Marlin steering and logfiles. 
  :@rawfile:       					Name of the raw file to be used for the telescope calibration  
  :@steerfiles     					Name of the directory with the geometry information
  :@caltag         					Name of the calibration tag, used to write/read to calbration information
  :@gearfile       					Name of the gear file in the steering files directory
  :@nevents       					Number of events used during the telescope calibration
  :@Use_clusterDB  					Switch to enable/disable generation of a clusterDB
  :@beamenergy     					Beam energy of the particle beam in GeV
  :@mcdata         					Switch to indicate Monte Carlo or real test beam data (True: MC data, False: Real TB data)
  :@Use_LongTelescopeCali  			Switch to enable/disable step by step correlation between upstream telescope tracks and downstream hits on individual sensors
  :@UseOuterPlanesForClusterDB		Switch to enable/disable usage of the outer telescope planes during the cluster calibration
  :author: ulf.stolzenberg@phys.uni-goettingen.de   
  """ 

  if rawfile == None:
    print("Rawfile name missing! Skip telescope calibration.")
    return None

  if steerfiles == None:
    print("Steerfiles name missing! Skip telescope calibration.")
    return None

  
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = Calibration(steerfiles=steerfiles, name=caltag + '-cal') 
  
  # Create list of calibration steps 
  calpaths = path_utils.create_x0analysis_calibration_paths(CalObj, rawfile, gearfile, nevents, Use_clusterDB, beamenergy, mcdata, Use_LongTelescopeCali, UseOuterPlanesForClusterDB, csvdata=csvdata)
  
  # Run the calibration steps 
  CalObj.calibrate(paths=calpaths,ifile=rawfile,caltag=caltag)

  # Create DQM plots 
  DQMplots.calibration_DQMPlots(name=caltag + '-cal')


# Perform the angle reconstruction of a single run
def reconstruct(rawfile=None, steerfiles=None, caltag='default', gearfile='gear.xml', nevents=-1, Use_SingleHitSeeding=True, Use_clusterDB=True, beamenergy=2.0, mcdata=False, csvdata=False):

  """
  Calibrates an misaligned tracking telescope from run data. 
  Creates a folder localDB/caltag in workspace containing 
  calibration results. 
  Creates a folder tmp-runs/name-sim/ and populates it with 
  Marlin steering and logfiles. 
  :@params					Consists of:
  :@rawfile:       			Name of the raw file to be used for angle reconstruction  
  :@steerfiles     			Name of the directory with the geometry information
  :@caltag         			Name of the calibration tag, used to read to calbration information
  :@gearfile       			Name of the gear file in the steering files directory
  :@nevents   				Number of events used
  :@Use_SingleHitSeeding	Enable single hit seeding along beam axis
  :@Use_clusterDB  			Switch to enable/disable generation of a clusterDB
  :@beamenergy     			Beam energy of the particle beam in GeV
  :@mcdata         			Switch to indicate Monte Carlo or real test beam data (True: MC data, False: Real TB data)
  :author: ulf.stolzenberg@phys.uni-goettingen.de   
  """ 
  
  if rawfile == None:
    print("Rawfile name missing! Skip angle reconstruction.")
    return None
  
  if steerfiles == None:
    print("Steerfiles name missing! Skip angle reconstruction.")
    return None
  
  # Set cal tag that includes run name
  name = os.path.splitext(os.path.basename(rawfile))[0] + '-' + caltag
  
  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=name )
  
  # Create reconstuction path
  recopath = path_utils.create_anglereco_path(RecObj, rawfile, gearfile, nevents, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata, csvdata=csvdata)  
  
  # Run the reconstuction  
  RecObj.reconstruct(paths=recopath,ifile=rawfile,caltag=caltag) 


# Create angle reconstruction DQM plots
def reconstruction_DQM(rawfile=None,caltag="default"):
  """
  Creates DQM plots for the scattering angle reconstruction of a raw file. 
  :@rawfile:       			Name of the raw file to be used for the telescope calibration  
  :@caltag         			Name of the calibration tag, used to write/read to calbration information
  :author: ulf.stolzenberg@phys.uni-goettingen.de   
  """ 

  # Set cal tag that includes run name
  name = os.path.splitext(os.path.basename(rawfile))[0] + '-' + caltag
  
  # Create DQM plots
  DQMplots.anglereco_DQMPlots(filepath='root-files/X0-{}.root'.format(name)) 


def targetalignment(rawfile=None, steerfiles=None, iteration=None, caltag="default", gearfile="gear.xml", nevents=1000000, Use_SingleHitSeeding=False, Use_clusterDB=True, beamenergy=2.0, mcdata=False, csvdata=False):
  """ 
  Starts the scattering angle reconstruction and vertex fit on the central target
  plane. Afterwards the mean vertex z position is set as the new target z position in
  the aligment DB file and the calibration tag is exported
    :@rawfile:      		Input file for the reconstruction 
    :@steerfiles    		Name of the directory with the geometry information
    :@iteration:    		Target alignment iteration counter 
    :@caltag:       		Name of the calibration tag, used to write/read to calbration information
    :@gearfile:     		Name of the gear file
    :@nevents       		Number of events
    :@Use_SingleHitSeeding 	Enable single hit seeding along beam axis
    :@Use_clusterDB 		Switch to enable/disable generation of a clusterDB
    :@beamenergy    		Beam energy of the particle beam in GeV
    :@mcdata        		Switch to indicate Monte Carlo or real test beam data (True: MC data, False: Real TB data)
    :author: ulf.stolzenberg@phys.uni-goettingen.de 
  """ 
   
  if rawfile == None:
    print("Rawfile name missing! Skip target alignment.")
    return None

  if steerfiles == None:
    print("Steerfiles name missing! Skip target alignment.")
    return None

  if iteration == None:
    print("Iteration counter missing! Skip target alignment.")
    return None
  
  if iteration == 0:
    oldcaltag=caltag
  else: 
    oldcaltag=caltag+'-target-align-it-{:d}'.format(iteration-1)
     
  newcaltag=caltag+'-target-align-it-{:d}'.format(iteration)
  
  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name='x0-reco-targetalign-it-{:d}'.format(iteration) )
  
  # Create reconstuction path
  recopath = path_utils.create_anglereco_path(RecObj, rawfile, gearfile, nevents, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata, csvdata=csvdata)  
  
  # Run the reconstuction  
  RecObj.reconstruct(paths=recopath,ifile=rawfile,caltag=oldcaltag)  
  
  # Read the vertex position and save it in the alignmentDB
  dbname=RecObj.create_dbfilename("alignmentDB.root")
  treename=RecObj.get_rootfilename('X0')
  gear.save_targetpos(treename,dbname)
  RecObj.export_caltag(newcaltag)


# Perform x0 calibration
def xx0calibration(RunList=None, steerfiles=None, caltag="default", x0caltag=None):
  """ 
  Starts the radiation length calibration procedure, which is used to determine
  calibration parameters such as the global offset of the telescope angle resolution, the global offset of the multiple
  scattering width and beam energy gradients. Creates a folder localDB/x0caltag in workspace containing 
  calibration results. Creates a folder tmp-runs/X0Calibration-x0caltag/ and populates it with 
  DQM Files and logfiles. 
    :@RunList:      List of runs used in the radiation length calibration 
    :@steerfiles    Name of the directory with the geometry information
    :@x0caltag:     Name of the x0 calibration tag, used to write/read to x0 calbration information
    :@caltag:       Name of the calibration tag, used to write/read to calbration information
    :author: ulf.stolzenberg@phys.uni-goettingen.de 
  """ 

  if RunList == None:
    print("RunList name missing! Skip x0 calibration.")
    return None

  if steerfiles == None:
    print("Steerfiles name missing! Skip x0 calibration.")
    return None

  if x0caltag == None:
    x0caltag=caltag

  # Create list with input root files from list of input raw files
  RootFileList_x0cali=[]
  X0Calibration.CreateRootFileList(rawlist=RunList,rootlist=RootFileList_x0cali, caltag=caltag)

  # Generate a uncalibrated X/X0 image
  imagenametag='X0image-calitarget-Uncalibrated'
  GenerateImage.x0imaging(rootfilelist=[RootFileList_x0cali[0]],caltag='',steerfiles=steerfiles,nametag=imagenametag)

  # Path to uncalibrated X0 image file
  imagefilename='/root-files/'+imagenametag+'.root'

  # Do a calibration of the angle resolution
  X0Calibration.x0calibration(rootfilelist=RootFileList_x0cali,imagefile=imagefilename,caltag=x0caltag,steerfiles=steerfiles)

  nametag='X0Calibration-'+x0caltag
  DQMplots.x0calibration_DQMPlots(nametag=nametag)


# Generate x0 image
def xx0image(RunList=None, steerfiles=None, caltag="default", listname="X0image", x0caltag=None):
  """ 
  Starts the radiation length imaging procedure. Produces radiation length images and populates root/ results/
  and tmp with folders according to the selected listname.
    :@RunList:      List of runs used in the radiation length calibration 
    :@steerfiles    Name of the directory with the geometry information
    :@x0caltag:     Name of the x0 calibration tag, used to write/read to x0 calbration information
    :@caltag:       Name of the calibration tag, used to write/read to calbration information
    :@listname:     Name of the radiation length image
    :author: ulf.stolzenberg@phys.uni-goettingen.de 
  """ 

  if RunList == None:
    print("RunList name missing! Skip x0 imaging.")
    return None

  if steerfiles == None:
    print("Steerfiles name missing! Skip x0 imaging.")
    return None

  if x0caltag == None:
    x0caltag=caltag

  RootFileList_x0image=[]
  X0Calibration.CreateRootFileList(rawlist=RunList,rootlist=RootFileList_x0image, caltag=caltag)

  if listname=='':
    print("No image name found. Using default naming scheme!")
    listname='X0image'

  # Determine name of image
  if caltag=='':
    nametag=listname+'-Uncalibrated'
  else:
    nametag=listname+'-Calibrated-'+x0caltag

  # Do a calibration of the angle resolution
  GenerateImage.x0imaging(rootfilelist=RootFileList_x0image,caltag=x0caltag,steerfiles=steerfiles,nametag=nametag)

  DQMplots.x0image_Plots(nametag=nametag)


