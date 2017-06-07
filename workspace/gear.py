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

def override_gear(gearfile=None, sensorID=None, parametername=None, value=None):
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
            if ID==sensorID:
              print('Changing ladder with ID '+ID)
              ladder.set(parametername, str(value))

          for sensitive in layer.findall('sensitive'):
            ID=sensitive.get('ID')
            if ID==sensorID:
              print('Changing sensitive volume with ID '+ID)
              sensitive.set(parametername, str(value))

  tree.write(gearfile) 

def randomize_gearpars(gearfile=None, parameter_tuple=None):
  """
  Overrides position and orientation of multiple sensors in gearfile with gaussian random variables
    :@gearfile:          gear file to be overwritten  
    :@parameter_tuple    tuple of parameter names and its mean and sigma values, tuple should be of the form ([sensorID1(int),name1(string),mean1(float),sigma1(float)],[sensorID2,name2,mean2,sigma2]....)
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """   

  for par in parameter_tuple[:]: 

    # Read out parameter name and mean and sigma
    sensorID=par[0]
    parametername=par[1]
    mean=par[2]
    sigma=par[3]

    # Calculate random shift for this parameter
    value=random.gauss(mean,sigma)
    print('Random value '+str(value))
    override_gear(gearfile=gearfile, sensorID=sensorID, parametername=parametername, value=value)

def set_globaloffset(gearfile=None, parametername=None, value):
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
            ID=ladder.get('ID')
            print('Changing ladder with ID '+ID)
            ladder.set(parametername, str(value))

          for sensitive in layer.findall('sensitive'):
            ID=sensitive.get('ID')
            print('Changing sensitive volume with ID '+ID)
            sensitive.set(parametername, str(value))

  tree.write(gearfile) 


     
