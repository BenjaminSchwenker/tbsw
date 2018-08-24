"""
Some helper code to manipulate gear files

:author: ulf.stolzenberg@phys.uni-goettinge.de  
"""

import os
import shutil
import subprocess
import glob
import xml.etree.ElementTree
import random
import math

from ROOT import TFile, TH1F, gROOT

try:
    from itertools import zip_longest as zip_longest
except:
    from itertools import izip_longest as zip_longest

def getCurrentValue(event,branchname):
  return getattr(event,branchname)

# Function, which saves the vertex z mean position from a given file as the new nominal position of the target in the alignment DB file
def save_targetpos(treefilename=None,dbfilename=None): 
  """
  Updates the alignment DB with the expected target position from the vertex fit
    :@filename:    Name of file with MSCTree which includes the vertex parameters 
    :@dbfilename   Name of the db root file 
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """ 

  if treefilename == None:
    return None

  if dbfilename == None:
    return None

  gROOT.SetBatch(1)
  rawfile = gROOT.FindObject( treefilename )

  # if tree file is open, close it
  if rawfile:
    rawfile.Close()

  # Open tree file in read mode
  rawfile = TFile( treefilename, 'READ' )
  tree=rawfile.Get('MSCTree')
  
  # Calculate mean value of the vertex z position
  count=0
  mean=0
 
  for event in tree:
    mean += getCurrentValue(event,"vertex_z")
    count +=1

  if count == 0:
    raise ValueError('Tree is empty!')

  # Round mean value to first decimal (100 microns)
  mean=round(mean/count,1)

  print "New target position is"
  print mean

  Modify_AlignmentDBFile(dbfilename=dbfilename, planenumber=3, mode='z', value=mean)


def Create_AlignmentDBFile_From_Gear(gearfile=None, truthdbfilename=None):
  """
  Writes sensor position and angles from gear file into truthdb root file 
    :@gearfile:       gear file to be copied into alignment root file  
    :@truthdbfilename      Name of the output root file 
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """   

  if gearfile == None:
    return None

  if truthdbfilename == None:
    return None

  # Open db file
  dbfile = TFile( truthdbfilename, 'RECREATE', 'alignment parameters from ' + gearfile )

  # Define lists of alignment parameters
  id_list = []
  xpos_list = []
  ypos_list = []
  zpos_list = []
  xrot_list = []
  yrot_list = []
  zrot_list = []

  tree = xml.etree.ElementTree.parse(gearfile)
  root = tree.getroot()

  # Read out the alignment parameters
  for detectors in root.findall('detectors'): 
    for detector in detectors.findall('detector'):
      for layers in detector.findall('layers'):
        for layer in layers.findall('layer'):

          for sensitive in layer.findall('sensitive'):
            xpos_list.append(float(sensitive.get('positionX')))
            ypos_list.append(float(sensitive.get('positionY')))
            zpos_list.append(float(sensitive.get('positionZ')))
            xrot_list.append(float(sensitive.get('alpha')))
            yrot_list.append(float(sensitive.get('beta')))
            zrot_list.append(float(sensitive.get('gamma')))
            id_list.append(int(sensitive.get('ID')))

  # Sort z position list and the corresponding sensor id list
  zpos_list2, id_list2 = (list(t) for t in zip(*sorted(zip(zpos_list, id_list))))

  # get number of planes
  nentries=len(id_list2)

  # ID histogram
  hSensorID = TH1F("hSensorID","",nentries,0,nentries)
  hSensorID.SetTitle("")
  hSensorID.GetXaxis().SetTitle("plane")
  hSensorID.GetYaxis().SetTitle("Sebsor ID") 

  # X position histogram
  hPositionX = TH1F("hPositionX","",nentries,0,nentries)
  hPositionX.SetTitle("")
  hPositionX.GetXaxis().SetTitle("plane")
  hPositionX.GetYaxis().SetTitle("position x [mm]") 

  # X position histogram
  hPositionY = TH1F("hPositionY","",nentries,0,nentries)
  hPositionY.SetTitle("")
  hPositionY.GetXaxis().SetTitle("plane")
  hPositionY.GetYaxis().SetTitle("position y [mm]")

  # Z position histogram
  hPositionZ = TH1F("hPositionZ","",nentries,0,nentries)
  hPositionZ.SetTitle("")
  hPositionZ.GetXaxis().SetTitle("plane")
  hPositionZ.GetYaxis().SetTitle("position z [mm]")

  # alpha rotation histogram
  hRotationAlpha = TH1F("hRotationAlpha","",nentries,0,nentries)
  hRotationAlpha.SetTitle("")
  hRotationAlpha.GetXaxis().SetTitle("plane")
  hRotationAlpha.GetYaxis().SetTitle("rotation alpha [rad]") 

  # beta rotation histogram
  hRotationBeta = TH1F("hRotationBeta","",nentries,0,nentries)
  hRotationBeta.SetTitle("")
  hRotationBeta.GetXaxis().SetTitle("plane")
  hRotationBeta.GetYaxis().SetTitle("rotation beta [rad]")

  # gamma rotation histogram
  hRotationGamma = TH1F("hRotationGamma","",nentries,0,nentries)
  hRotationGamma.SetTitle("")
  hRotationGamma.GetXaxis().SetTitle("plane")
  hRotationGamma.GetYaxis().SetTitle("rotation gamma [rad]")

  # Loop over sensor ids
  for bin,sensid in enumerate(id_list2):

    # Find list index for this sensor id
    index = id_list.index(sensid)
  
    # Fill histograms
    hSensorID.SetBinContent(bin+1,id_list[index])
    hPositionX.SetBinContent(bin+1,xpos_list[index])
    hPositionY.SetBinContent(bin+1,ypos_list[index])
    hPositionZ.SetBinContent(bin+1,zpos_list[index])
    hRotationAlpha.SetBinContent(bin+1,xrot_list[index]/180*3.1415)   # angles in gear file are given in degree -> change to rad
    hRotationBeta.SetBinContent(bin+1,yrot_list[index]/180*3.1415)    # angles in gear file are given in degree -> change to rad
    hRotationGamma.SetBinContent(bin+1,zrot_list[index]/180*3.1415)   # angles in gear file are given in degree -> change to rad
  
  dbfile.Write()
  dbfile.Close()


def set_parameter(gearfile=None, sensorID=None, parametername=None, value=None):
  """
  Overrides value field in all sensors with a specific sensor ID in gearfile
    :@gearfile:       gear file to be overwritten  
    :@sensorID        Sensor ID with the parameter to be modified
    :@parametername   Parameter to be modified     
    :@value:          insert string value in field value 
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """   
  tree = xml.etree.ElementTree.parse(gearfile)
  root = tree.getroot() 
  
  for detectors in root.findall('detectors'): 
    for detector in detectors.findall('detector'):
      for layers in detector.findall('layers'):
        for layer in layers.findall('layer'):

          for ladder in layer.findall('ladder'):
            ID=ladder.get('ID')
            if ID==str(sensorID):
              print('[INFO] Changing plane with ID '+ID)
              ladder.set(parametername, str(value))

          for sensitive in layer.findall('sensitive'):
            ID=sensitive.get('ID')
            if ID==str(sensorID):
              sensitive.set(parametername, str(value))

  tree.write(gearfile) 


def add_offset(gearfile=None, sensorID=None, parametername=None, value=None):
  """
  Adds value to parameter in all sensors with a specific sensor ID in gearfile
    :@gearfile:       gear file to be overwritten  
    :@sensorID        Sensor ID with the parameter to be modified
    :@parametername   Parameter to be modified     
    :@value:          insert string value in field value 
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """   
  tree = xml.etree.ElementTree.parse(gearfile)
  root = tree.getroot() 
  
  for detectors in root.findall('detectors'): 
    for detector in detectors.findall('detector'):
      for layers in detector.findall('layers'):
        for layer in layers.findall('layer'):
          for ladder in layer.findall('ladder'):
            ID=ladder.get('ID')
            if ID==str(sensorID):
              #print('Changing ladder with ID '+ID)
              if ladder.get(parametername) is None:
                #print('Parameter with name '+parametername+' doesnt exist in ladder!')
                pass
              else:
                laddervalue=float(value)+float(ladder.get(parametername))
                ladder.set(parametername, str(laddervalue))


          for sensitive in layer.findall('sensitive'):
            ID=sensitive.get('ID')
            if ID==str(sensorID):
              #print('Changing sensitive volume with ID '+ID)
              if sensitive.get(parametername) is None:
                #print('Parameter with name '+parametername+' doesnt exist in sensitive!')
                pass
              else:
                sensvalue=float(value)+float(sensitive.get(parametername))
                sensitive.set(parametername, str(sensvalue))


  tree.write(gearfile) 


def randomize_gearparameter(gearfile=None, sensorID=None, parametername=None, mean=None, sigma=None):
  """
  Overrides position and orientation of a single sensor in gearfile with gaussian random variable
    :@gearfile:          gear file to be overwritten  
    :@sensorID           Sensor ID with the parameter to be modified
    :@parametername      Parameter to be modified     
	:@mean				 Mean of the gaussian random distribution
    :@sigma              Sigma of the gaussian random distribution
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """   

  # Calculate random shift for this parameter
  value=random.gauss(mean,sigma)
  add_offset(gearfile=gearfile, sensorID=sensorID, parametername=parametername, value=value)


def randomize_telescope(gearfile=None, mean_list=[0.0,0.0,0.0,0.0,0.0,0.0], sigma_list=[0.0,0.0,0.0,0.0,0.0,0.0], sensorexception_list=[], modeexception_list=[]):
  """
  Overrides position and orientation of all sensors in gearfile with gaussian random variables
    :@gearfile:               gear file to be overwritten  
	:@mean_list               List of mean values, should contain 6 entries: 1 for each alignment paramter [mean_x,mean_y,mean_z,mean_alpha,mean_beta,mean_gamma]
    :@sigma_list              List of sigma values, should contain 6 entries: 1 for each alignment parameter [sigma_x,sigma_y,sigma_z,sigma_alpha,sigma_beta,sigma_gamma]
    :@sensorexception_list    List of sensor ids, which will not be misaligned
    :@modeexception_list      List of misalignment modes, which are not used
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """  

  tree = xml.etree.ElementTree.parse(gearfile)
  root = tree.getroot()  

  # This is the default mode list: Containing all the positions and angles
  default_mode_list=['positionX','positionY','positionZ','alpha','beta','gamma']

  # Combine the default_mode_list and the mean and sigma lists to tuples
  parameter_tuple_list=zip_longest(default_mode_list, mean_list, sigma_list, fillvalue=0.0)

  # Run over the tuple list and remove tuples, that have a first element that is also present in the modeexception_list
  for modeexception in modeexception_list:
      parameter_tuple_list = filter(lambda x: x[0] != modeexception,parameter_tuple_list)

  for detectors in root.findall('detectors'): 
    for detector in detectors.findall('detector'):
      for layers in detector.findall('layers'):
        for layer in layers.findall('layer'):

          for ladder in layer.findall('sensitive'):
            ID=ladder.get('ID')

            # Check whether the sensor exception list contains this id
            if sensorexception_list.count(int(ID)) < 1:
              print('[INFO] Misalign position of plane with ID '+ID)
              # Loop over the alignment parameters in the default mode list
              for parameter_tuple in parameter_tuple_list:
                #Randomize the mode for this specific mode and sensor ID
                randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername=parameter_tuple[0], mean=parameter_tuple[1], sigma=parameter_tuple[2])


def set_globalparameter(gearfile=None, parametername=None, value=None):
  """
  Overrides all parameters with a certain name on all planes
    :@gearfile:          gear file to be overwritten  
    :@parametername   Parameter to be modified     
    :@value:          insert string value in field value 
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """ 

  tree = xml.etree.ElementTree.parse(gearfile)
  root = tree.getroot() 
  
  for detectors in root.findall('detectors'): 
    for detector in detectors.findall('detector'):
      for layers in detector.findall('layers'):
        for layer in layers.findall('layer'):

          for ladder in layer.findall('ladder'):
            ladder.set(parametername, str(value))

          for sensitive in layer.findall('sensitive'):
            sensitive.set(parametername, str(value))

  tree.write(gearfile) 


def add_weakmode(gearfile=None, parametername=None, value=None):
  """
  Adds a weak mode to a certain parameter on all planes
  The offset on each plane depends on the z position
  
    :@gearfile:          gear file to be overwritten  
    :@parametername      Parameter to be modified     
    :@value:             max offset on last plane
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """ 

  tree = xml.etree.ElementTree.parse(gearfile)
  root = tree.getroot() 

  maxvalue=value

  # List of z positions
  zPositions=[]
  IDs=[]
  # loop over ladders to find the z positions
  for detectors in root.findall('detectors'): 
    for detector in detectors.findall('detector'):
      for layers in detector.findall('layers'):
        for layer in layers.findall('layer'):

          for ladder in layer.findall('ladder'):
            laddervalue=float(value)+float(ladder.get(parametername))
            zPositions.append(float(ladder.get('positionZ')))
            IDs.append(ladder.get('ID'))

  # Get the max z value in the z position list
  maxzpos=max(zPositions)

  # Calculate the parameter change per mm in z direction
  slope=float(value)/maxzpos

  values=[]

  # Calculate the values for the other planes
  for zPosition in zPositions:
    values.append(float(zPosition*slope))

  planenumber=0
  for sensID in IDs:
    set_parameter(gearfile=gearfile, sensorID=sensID, parametername=parametername, value=values[planenumber])
    planenumber=planenumber+1
    

def rotate_telescope(gearfile=None, parametername=None, valueZ=None, valueAngle=None):
  """
  Rotates the center of the telescope planes in the x-z or y-z plane by an angle valueAngle
  The plane is determined by the parametername, in case the parameter name isn't positionX or positionY
  nothing happens. The axis is located on the z axis at valueZ. 
  
    :@gearfile:          gear file to be overwritten  
    :@parametername      Parameter, which defines the rotation plane    
    :@valueZ:            z position of rotation axis
    :valueAngle          Rotation angle in degrees
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """ 

  tree = xml.etree.ElementTree.parse(gearfile)
  root = tree.getroot() 

  # Read out positions of telescope before rotation
  xDistances=[]
  yDistances=[]
  zPositions=[]
  IDs=[]

  # angle in radians
  valueAngle = math.radians(valueAngle)

  # loop over ladders to find the z positions, IDs
  for detectors in root.findall('detectors'): 
    for detector in detectors.findall('detector'):
      for layers in detector.findall('layers'):
        for layer in layers.findall('layer'):
          for ladder in layer.findall('ladder'):
            xDistances.append(float(ladder.get('positionX')))
            yDistances.append(float(ladder.get('positionY')))
            zPositions.append(float(ladder.get('positionZ')))
            IDs.append(ladder.get('ID'))

  # Calculate the distance between the rotaion axis position valueZ and the z positions of the sensors
  zDistances=[]
  for zPosition in zPositions:
     zDistances.append(zPosition-valueZ)

  planenumber=0

  # Calculate the positions of the  sensor center
  if(parametername == 'positionX'):
    set_globalparameter(gearfile=gearfile, parametername='alpha', value=valueAngle)
    for zDistance in zDistances:

      # Calculate position after rotation in XZ plane
      z=math.cos(valueAngle)*zDistance-math.sin(valueAngle)*xDistances[planenumber]
      x=math.sin(valueAngle)*zDistance+math.cos(valueAngle)*xDistances[planenumber]

      # Set the calculated parameters in the gear file
      set_parameter(gearfile=gearfile, sensorID=IDs[planenumber], parametername='positionX', value=x)
      set_parameter(gearfile=gearfile, sensorID=IDs[planenumber], parametername='positionZ', value=z+valueZ)

      planenumber=planenumber+1
    
    
  # Calculate the positions of the  sensor center    
  elif(parametername == 'positionY'):
    set_globalparameter(gearfile=gearfile, parametername='beta', value=valueAngle)
    for zDistance in zDistances:

      # Calculate position after rotation in YZ plane
      z=math.cos(valueAngle)*zDistance-math.sin(valueAngle)*xDistances[planenumber]
      y=math.sin(valueAngle)*zDistance+math.cos(valueAngle)*xDistances[planenumber]

      # Set the calculated parameters in the gear file
      set_parameter(gearfile=gearfile, sensorID=IDs[planenumber], parametername='positionY', value=y)
      set_parameter(gearfile=gearfile, sensorID=IDs[planenumber], parametername='positionZ', value=z+valueZ)

      planenumber=planenumber+1
  else:
    print('Rotation plane not well defined!')


def Modify_AlignmentDBFile(dbfilename=None, planenumber=None, mode=None, value=None):
  """
  Change single position or angle in alignment DB file
    :@dbfilename      Name of the output root file 
    :@planenumber     Planenumber to be modified (starting with 0)
    :@mode            Mode to be modified ['x','y','z','alpha','beta','gamma']
    :@value           New value in mm or rad
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """

  if dbfilename == None:
    return None

  if planenumber == None:
    return None

  if mode == None:
    return None

  if value == None:
    return None
    
  if os.path.isfile(dbfilename):
    dbfile = TFile( dbfilename, 'UPDATE' )
  else: 
    raise ValueError('alignment DB ('+dbfilename+') file not found') 

  # Get access to histogram  

  if mode=='x':
    histo = dbfile.Get("hPositionX")
    dbfile.Delete("hPositionX;1")

  elif mode=='y':
    histo = dbfile.Get("hPositionY")
    dbfile.Delete("hPositionY;1")

  elif mode=='z':
    histo = dbfile.Get("hPositionZ")
    dbfile.Delete("hPositionZ;1")

  elif mode=='alpha':
    histo = dbfile.Get("hRotationAlpha")
    dbfile.Delete("hRotationAlpha:1")

  elif mode=='beta':
    histo = dbfile.Get("hRotationBeta")
    dbfile.Delete("hRotationBeta;1")

  elif mode=='gamma':
    histo = dbfile.Get("hRotationGamma")
    dbfile.Delete("hRotationGamma;1")

  else:
    return None

  histo.SetBinContent(planenumber+1,value)
  
  dbfile.Write()
  dbfile.Close()


     
