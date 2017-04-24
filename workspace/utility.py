#!/usr/bin/python3

"""
Utility code to use TBSW in pyhton scripts 

:author: benjamin.schwenker@phys.uni-goettinge.de  
"""

import os
import shutil
import subprocess
import sys, getopt 
import glob

class Environment(object):
  """
  Class which implements an environment for executing Marlin with all needed 
  steering and config files. 
  
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self, tag, steerfiles):
    """
    Parameters
    ----------
    tag : string 
        name of the temporary folder holding the environment
    steerfiles : string 
        name of the folder holding steer and config files  
    """
    self.tag = tag
    self.fullpath = os.getcwd() 
    self.tmpdir = os.path.join(fullpath+'/tmp-runs',tag)  
    
    # create an empty tmpdir
    if not os.path.isdir(fullpath+'/tmp-runs'):
      os.mkdir(fullpath+'/tmp-runs')
    if os.path.isdir(self.tmpdir):
      shutil.rmtree(self.tmpdir)
    os.mkdir(self.tmpdir)   
    
    # copy the steer files 
    if not os.path.isdir(steerfiles):
      print ('[PRINT] error: no steerfiles found')        
    else: 
      shutil.copytree(steerfiles,self.tmpdir)
      
  def add_caltag(self, caltag):  
    # check that calibration files exist
    caldir = self.fullpath+'/cal-files/'+caltag   
    if not os.path.isdir(caldir):
      print ('[Print] warning: no calibration data found')  
    else: 
      # copy calibration files 
      shutil.copytree(caldir,self.tmpdir+'/cal-files')
  
  def create_caltag(self, caltag):
    # create new calibration tag    
    if not os.path.isdir(self.fullpath+'/cal-files'):
      os.mkdir(self.fullpath+'/cal-files')
    caldir = self.fullpath+'/cal-files/'+caltag  
    if os.path.isdir(caldir):
      shutil.rmtree(caldir)   
    os.mkdir(caldir)		
    
    # popule caltag with DB files    
    for dbfile in glob.glob('*DB*'): 
      shutil.copy(dbfile, os.path.join(caldir,dbfile))  
    # treat the gear file as a calibration file as well
    shutil.copy('gear.xml', os.path.join(caldir,'gear.xml'))  
  	                           
  def add_inputfile(self, inputfile):
    os.symlink( os.path.join(self.fullpath, inputfile), self.tmpdir+'/inputfilename') 

  def execute(self,path):
    # execute in tmp dir  
    os.chdir(tmpdir)
    
    for xmlfile in path:
      action = '/$MARLIN/bin/Marlin ' + xmlfile + ' > log-' + xmlfile + ' 2>&1'
      subprocess.call(action, shell=True)
    
    os.chdir(fullpath) 
  
  def cleanup(self):
    for tmpfile in glob.glob('tmp*'):
      os.remove(tmpfile)

  def copy_rootfiles(self):
    # store dqm files 
    if not os.path.isdir(fullpath+'/root-files'):
      os.mkdir(fullpath+'/root-files')
    
    for dqmfile in glob.glob('*.root'): 
      name = os.path.splitext(os.path.basename(dqmfile))[0]
      shutil.move(dqmfile, fullpath+'/root-files/'+name+'-'+runtag+'-'+caltag+'.root')  
    


def simulate(steerfiles, marlinpath, caltag, ofile):
  """
  Creates a lcio file containing N events. The events contain collections for 
  simulated tracks measured signals. 
  
  :@steerfile:  name of folder containing all steering files  
  :@marlinpath: list containing Marlin xml files that will be executed 
  :@caltag:     name of calibration tag (optional)
  :@ofile:      name of output lcio file
  
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """

  return None


def reconstruct(steerfiles, marlinpath, caltag, ifile):
  """
  Reconstructs an input file with raw data using a caltag for calibration. 
  
  :@steerfile:  name of folder containing all steering files 
  :@marlinpath: list containing Marlin xml files that will be executed 
  :@caltag:     name of calibration tag (optional)
  :@ifile:      name of input file with raw data
  
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """   

  return None


def calibrate(steerfiles, marlinpath, caltag, ifile):
  """
  Calibrate raw data from a test beam experiment. 
  
  :@steerfile:  name of folder containing all steering files 
  :@marlinpath: list containing Marlin xml files that will be executed 
  :@caltag:     name of calibration tag (optional)
  :@ifile:      name of input file with raw data
  
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """     

  return None

