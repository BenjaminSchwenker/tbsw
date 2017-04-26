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
  def __init__(self, name, steerfiles):
    """
    Parameters
    ----------
    name : string 
        name of the temporary folder holding the environment
    steerfiles : string 
        name of the folder holding steer and config files  
    """
    self.name = name
    self.runtag = ''
    self.caltag = ''
    self.cwdir = os.getcwd() 
    self.tmpdir = os.path.join(self.cwdir+'/tmp-runs',self.name)  
    
    # create an empty tmpdir
    if not os.path.isdir(self.cwdir+'/tmp-runs'):
      os.mkdir(self.cwdir+'/tmp-runs')
    if os.path.isdir(self.tmpdir):
      shutil.rmtree(self.tmpdir)
    os.mkdir(self.tmpdir)   
    
    # copy the steer files 
    if not os.path.isdir(steerfiles):
      print ('[PRINT] error: no steerfiles found')        
    else: 
      shutil.copytree(steerfiles,self.tmpdir)
      
  def add_caltag(self, caltag):  
    self.caltag = caltag
    # check that calibration files exist
    caldir = self.cwdir+'/cal-files/'+caltag   
    if not os.path.isdir(caldir):
      print ('[Print] warning: no calibration data found')  
    else: 
      # copy calibration files 
      shutil.copytree(caldir,self.tmpdir+'/cal-files')
    
  
  def create_caltag(self, caltag):
    self.caltag = caltag
    # create new calibration tag    
    if not os.path.isdir(self.cwdir+'/cal-files'):
      os.mkdir(self.cwdir+'/cal-files')
    caldir = self.cwdir+'/cal-files/'+caltag  
    if os.path.isdir(caldir):
      shutil.rmtree(caldir)   
    os.mkdir(caldir)		
    
    # popule caltag with DB files    
    for dbfile in glob.glob('*DB*'): 
      shutil.copy(dbfile, os.path.join(caldir,dbfile))  
    # treat the gear file as a calibration file as well
    shutil.copy('gear.xml', os.path.join(caldir,'gear.xml'))  
  	                           
  def link_input(self, inputfile):
    self.runtag = os.path.splitext(os.path.basename(inputfile))[0]
    os.symlink( os.path.join(self.cwdir, inputfile), self.tmpdir+'/inputfilename') 

  def run(self,path):
    # run Marlin in tmpdir  
    os.chdir(self.tmpdir)
    
    for xmlfile in path:
      logfile = os.path.splitext( os.path.basename(xmlfile))[0] + '.log'
      action = '/$MARLIN/bin/Marlin ' + xmlfile + ' > ' + logfile + ' 2>&1'
      subprocess.call(action, shell=True)
      print ('[Print] Executing Marlin ' + xmlfile + ' is done')    
    
    # remove tmp* files 
    for tmpfile in glob.glob('tmp*'):
      os.remove(tmpfile)
    
    # go back to workspace
    os.chdir(self.cwdir) 
  
  def copy_rootfiles(self):
    # cd into tmpdir  
    os.chdir(self.tmpdir)
     
    # copy root files 
    if not os.path.isdir(self.cwdir+'/root-files'):
      os.mkdir(self.cwdir+'/root-files')
    
    for rootfile  in glob.glob('*.root'): 
      name = os.path.splitext(os.path.basename(rootfile))[0]
      shutil.move(rootfile, self.cwdir+'/root-files/'+name+'-'+self.runtag+'-'+self.caltag+'.root')  
    


def simulate(steerfiles, path, caltag, ofile):
  """
  Creates a lcio file containing N events. The events contain collections for 
  simulated tracks measured signals. 
  
  :@steerfile:  name of folder containing all steering files  
  :@path:       list containing Marlin xml files that will be executed 
  :@caltag:     name of calibration tag (optional)
  :@ofile:      name of output lcio file
  
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  
  name = os.path.splitext(os.path.basename(ofile))[0] + '-' + caltag + '-sim'
  env = Environment(name=name, steerfiles=steerfiles)  
  env.add_caltag(caltag)
  env.run(path)
  
  src =  os.path.join(env.tmpdir, 'outputfile.slcio')	
  dest = os.path.join(env.cwdir, ofile)	
  shutil.move(src, dest)
  
  return None
  
   
def reconstruct(steerfiles, path, caltag, ifile):
  """
  Reconstructs an input file with raw data using a caltag for calibration. 
  
  :@steerfile:  name of folder containing all steering files 
  :@path:       list containing Marlin xml files that will be executed 
  :@caltag:     name of calibration tag (optional)
  :@ifile:      name of input file with raw data
  
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """   
  
  print ('[Print] Starting to process file ' + ifile + ' ...')  
  
  name = os.path.splitext(os.path.basename(ifile))[0] + '-' + caltag + '-reco'
  env = Environment(name=name, steerfiles=steerfiles)  
  env.add_caltag(caltag)
  env.link_input(ifile)
  env.run(path)
  env.copy_rootfiles()

  print ('[Print] Done processing file ' + ifile)  

  return None


def calibrate(steerfiles, path, caltag, ifile):
  """
  Calibrate raw data from a test beam experiment. 
  
  :@steerfile:  name of folder containing all steering files 
  :@marlinpath: list containing Marlin xml files that will be executed 
  :@caltag:     name of calibration tag (optional)
  :@ifile:      name of input file with raw data
  
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """   
  
  print ('[Print] Starting to process file ' + ifile + ' ...')    
  
  name = os.path.splitext(os.path.basename(ifile))[0] + '-' + caltag + '-cal'
  env = Environment(name=name, steerfiles=steerfiles)  
  env.link_input(ifile)
  env.run(path)
  env.create_caltag(caltag) 
  
  print ('[Print] Done processing file ' + ifile)  
     
  return None

