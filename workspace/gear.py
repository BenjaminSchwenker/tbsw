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
              print('Changing ladder with ID '+ID)
              ladder.set(parametername, str(value))

          for sensitive in layer.findall('sensitive'):
            ID=sensitive.get('ID')
            if ID==str(sensorID):
              print('Changing sensitive volume with ID '+ID)
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
              print('Changing ladder with ID '+ID)
              laddervalue=float(value)+float(ladder.get(parametername))
              ladder.set(parametername, str(laddervalue))

          for sensitive in layer.findall('sensitive'):
            ID=sensitive.get('ID')
            if ID==str(sensorID):
              print('Changing sensitive volume with ID '+ID)
              sensvalue=float(value)+float(sensitive.get(parametername))
              sensitive.set(parametername, str(sensvalue))

  tree.write(gearfile) 


def randomize_gearparameter(gearfile=None, sensorID=None, parametername=None, mean=None, sigma=None):
  """
  Overrides position and orientation of multiple sensors in gearfile with gaussian random variables
    :@gearfile:          gear file to be overwritten  
    :@sensorID           Sensor ID with the parameter to be modified
    :@parametername      Parameter to be modified     
	:@mean				 Mean of the gaussian random distribution
    :@sigma              Sigma of the gaussian random distribution
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """   

  # Calculate random shift for this parameter
  value=random.gauss(mean,sigma)
  print('Random value '+str(value))
  add_offset(gearfile=gearfile, sensorID=sensorID, parametername=parametername, value=value)


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
    

  

  


     
