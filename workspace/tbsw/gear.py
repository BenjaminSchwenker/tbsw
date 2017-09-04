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

from ROOT import TFile, TH1F


def Create_AlignmentDBFile_From_Gear(gearfile=None, dbfilename="alignmentDB_simulation.root"):
  """
  Overrides value field in all sensors with a specific sensor ID in gearfile
    :@gearfile:       gear file to be copied into alignment root file  
    :@dbfilename      Name of the output root file 
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """   

  if gearfile == None:
    return None

  # Open db file
  dbfile = TFile( dbfilename, 'RECREATE', 'alignment parameters from ' + gearfile )

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


def randomize_telescope(gearfile=None, mean_pos=None, sigma_pos=None, mean_rot=None, sigma_rot=None):
  """
  Overrides position and orientation of all sensors in gearfile with gaussian random variables
    :@gearfile:          gear file to be overwritten  
	:@mean_pos		     Mean of the gaussian random distribution, which determine the plane positions
    :@sigma_pos          Sigma of the gaussian random distribution, which determine the plane positions
	:@mean_rot		     Mean of the gaussian random distribution, which determine the plane rotations
    :@sigma_rot          Sigma of the gaussian random distribution, which determine the plane rotations
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """  

  tree = xml.etree.ElementTree.parse(gearfile)
  root = tree.getroot()  

  for detectors in root.findall('detectors'): 
    for detector in detectors.findall('detector'):
      for layers in detector.findall('layers'):
        for layer in layers.findall('layer'):

          for ladder in layer.findall('sensitive'):
            ID=ladder.get('ID')
            print('[INFO] Misalign position of plane with ID '+ID)
            randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='positionX', mean=mean_pos, sigma=sigma_pos)
            randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='positionY', mean=mean_pos, sigma=sigma_pos)
            randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='positionZ', mean=mean_pos, sigma=sigma_pos)

            randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='alpha', mean=mean_rot, sigma=sigma_rot)
            randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='beta', mean=mean_rot, sigma=sigma_rot)
            randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='gamma', mean=mean_rot, sigma=sigma_rot)


def randomize_singleparameter_telescope(gearfile=None, mean=None, sigma=None, parametername=None):
  """
  Overrides parameter of all sensors in gearfile with gaussian random variables
    :@gearfile:          gear file to be overwritten  
	:@mean		         Mean of the gaussian random distribution
    :@sigma              Sigma of the gaussian random distribution
    :@parametername      Name of the parameter, which will be changes
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """  

  tree = xml.etree.ElementTree.parse(gearfile)
  root = tree.getroot()  

  for detectors in root.findall('detectors'): 
    for detector in detectors.findall('detector'):
      for layers in detector.findall('layers'):
        for layer in layers.findall('layer'):

          for ladder in layer.findall('sensitive'):
            ID=ladder.get('ID')
            print('[INFO] Misalign parameter: '+parametername+' of plane with ID '+ID)
            randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername=parametername, mean=mean, sigma=sigma)


def randomize_innertelescope(gearfile=None, mean_pos=None, sigma_pos=None, mean_rot=None, sigma_rot=None):
  """
  Overrides position and orientation of all sensors, except the sensors with ID 0 and 5 in gearfile with gaussian random variables
    :@gearfile:          gear file to be overwritten  
	:@mean_pos		     Mean of the gaussian random distribution, which determine the plane positions
    :@sigma_pos          Sigma of the gaussian random distribution, which determine the plane positions
	:@mean_rot		     Mean of the gaussian random distribution, which determine the plane rotations
    :@sigma_rot          Sigma of the gaussian random distribution, which determine the plane rotations
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """  

  tree = xml.etree.ElementTree.parse(gearfile)
  root = tree.getroot()  

  for detectors in root.findall('detectors'): 
    for detector in detectors.findall('detector'):
      for layers in detector.findall('layers'):
        for layer in layers.findall('layer'):

          for ladder in layer.findall('sensitive'):
            ID=ladder.get('ID')
            if ID != '0' and ID != '5':
            	print('[INFO] Misalign position of plane with ID '+ID)
            	randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='positionX', mean=mean_pos, sigma=sigma_pos)
            	randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='positionY', mean=mean_pos, sigma=sigma_pos)
            	randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='positionZ', mean=mean_pos, sigma=sigma_pos)

            	randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='alpha', mean=mean_rot, sigma=sigma_rot)
            	randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='beta', mean=mean_rot, sigma=sigma_rot)
            	randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='gamma', mean=mean_rot, sigma=sigma_rot)


def randomize_telescope_withsensorexception(gearfile=None, mean_pos=None, sigma_pos=None, mean_rot=None, sigma_rot=None, exception_ID=None):
  """
  Overrides position and orientation of all sensors, except the sensor with ID 'exception_ID' in gearfile with gaussian random variables
    :@gearfile:          gear file to be overwritten  
	:@mean_pos		     Mean of the gaussian random distribution, which determine the plane positions
    :@sigma_pos          Sigma of the gaussian random distribution, which determine the plane positions
	:@mean_rot		     Mean of the gaussian random distribution, which determine the plane rotations
    :@sigma_rot          Sigma of the gaussian random distribution, which determine the plane rotations
    :@exception_ID       ID of sensor that will not be misaligned
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """  

  tree = xml.etree.ElementTree.parse(gearfile)
  root = tree.getroot()  

  for detectors in root.findall('detectors'): 
    for detector in detectors.findall('detector'):
      for layers in detector.findall('layers'):
        for layer in layers.findall('layer'):

          for ladder in layer.findall('sensitive'):
            ID=ladder.get('ID')
            if ID != str(exception_ID):
            	print('[INFO] Misalign position of plane with ID '+ID)
            	randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='positionX', mean=mean_pos, sigma=sigma_pos)
            	randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='positionY', mean=mean_pos, sigma=sigma_pos)
            	randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='positionZ', mean=mean_pos, sigma=sigma_pos)

            	randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='alpha', mean=mean_rot, sigma=sigma_rot)
            	randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='beta', mean=mean_rot, sigma=sigma_rot)
            	randomize_gearparameter(gearfile=gearfile, sensorID=ID, parametername='gamma', mean=mean_rot, sigma=sigma_rot)


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


def shift_telescope(gearfile=None, valueX=None, valueY=None, valueZ=None):
  """
  Adds an offset to a certain parameter on all planes
    :@gearfile:          gear file to be overwritten     
    :@valueX:            X shift
    :@valueY:            Y shift
    :@valueZ:            Z shift	
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """ 

  tree = xml.etree.ElementTree.parse(gearfile)
  root = tree.getroot() 
  
  for detectors in root.findall('detectors'): 
    for detector in detectors.findall('detector'):
      for layers in detector.findall('layers'):
        for layer in layers.findall('layer'):

          for ladder in layer.findall('ladder'):
            laddervalueX=float(valueX)+float(ladder.get('positionX'))
            ladder.set('positionX', str(laddervalueX))

            laddervalueY=float(valueY)+float(ladder.get('positionY'))
            ladder.set('positionY', str(laddervalueY))

            laddervalueZ=float(valueZ)+float(ladder.get('positionZ'))
            ladder.set('positionZ', str(laddervalueZ))

          for sensitive in layer.findall('sensitive'):
            sensvalueX=float(valueX)+float(sensitive.get('positionX'))
            sensitive.set('positionX', str(sensvalueX))

            sensvalueY=float(valueY)+float(sensitive.get('positionY'))
            sensitive.set('positionY', str(sensvalueY))

            sensvalueZ=float(valueZ)+float(sensitive.get('positionZ'))
            sensitive.set('positionZ', str(sensvalueZ))

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


    


  

  


     
