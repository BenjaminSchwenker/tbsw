#!/usr/bin/python3

"""
Utility code to use TBSW in pyhton scripts 

:author: benjamin.schwenker@phys.uni-goettinge.de  
"""

import os
import shutil
import subprocess
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
    
    # create tmp-runs if not exist
    if not os.path.isdir(self.cwdir+'/tmp-runs'):
      os.mkdir(self.cwdir+'/tmp-runs')

    # remove old tmpdir if exists 
    if os.path.isdir(self.tmpdir):
      shutil.rmtree(self.tmpdir)
    
    # create tmpdir (containing steer files)
    if os.path.isdir(steerfiles):
      shutil.copytree(steerfiles,self.tmpdir)   
    else: 
      raise ValueError('No steerfiles can be found')
     
  def add_caltag(self, caltag):  
    self.caltag = caltag
    # check that calibration files exist
    caldir = self.cwdir+'/cal-files/'+caltag   
    if os.path.isdir(caldir):
      # copy calibration files 
      shutil.copytree(caldir,self.tmpdir+'/cal-files')  
    else: 
      print ('[INFO] Caltag not found') 
    
  
  def create_caltag(self, caltag):
    self.caltag = caltag
    caldir = self.cwdir+'/cal-files/'+caltag
    
    # create folder cal-files if not exist    
    if not os.path.isdir(self.cwdir+'/cal-files'):
      os.mkdir(self.cwdir+'/cal-files')
    # overwrite caltag if exists     
    if os.path.isdir(caldir):
      shutil.rmtree(caldir)
    
    #create caltag and populate with DB files    
    os.mkdir(caldir)		 
    for dbfile in glob.glob(self.tmpdir + '/*DB*'): 
      shutil.copy(dbfile, os.path.join(caldir,os.path.basename(dbfile)))  
    
    print ('[INFO] Created new caltag ', caltag) 
     
  	                           
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
      print ('[INFO] Marlin ' + xmlfile + ' is done')    
    
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
  
  print ('[INFO] Starting to process file ' + ifile + ' ...')  
  
  name = os.path.splitext(os.path.basename(ifile))[0] + '-' + caltag + '-reco'
  env = Environment(name=name, steerfiles=steerfiles)  
  env.add_caltag(caltag)
  env.link_input(ifile)
  env.run(path)
  env.copy_rootfiles()

  print ('[INFO] Done processing file ' + ifile)  

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
  
  print ('[INFO] Starting to process file ' + ifile + ' ...')    
  
  name = os.path.splitext(os.path.basename(ifile))[0] + '-' + caltag + '-cal'
  env = Environment(name=name, steerfiles=steerfiles)  
  env.link_input(ifile)
  env.run(path)
  env.create_caltag(caltag) 
  
  print ('[INFO] Done processing file ' + ifile)  
     
  return None

