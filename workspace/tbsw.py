"""
Some helper code to use TBSW with pyhton scripts 

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
  def __init__(self, name='default', steerfiles=None):
    """
    Parameters
    ----------
    name : string 
        name of the temporary folder holding copies of config files  
    steerfiles : string 
        name of the folder holding config files  
    """
    if steerfiles==None:
      raise ValueError('Parameter steerfiles is missing')
    
    self.name = name
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
    
  def add_caltag(self, caltag):  
    # check that calibration files exist
    caldir = self.cwdir+'/cal-files/'+caltag   
    if os.path.isdir(caldir):
      # copy calibration files 
      shutil.copytree(caldir,self.tmpdir+'/cal-files')  
    else: 
      print ('[INFO] Caltag not found') 
    
  def create_caltag(self, caltag):
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
    if os.path.isfile(inputfile):
      os.symlink( os.path.join(self.cwdir, inputfile), self.tmpdir+'/inputfilename') 
    else: 
      raise ValueError('No input file found')
  
  def get_name(self, filename):
    localname = os.path.join(self.tmpdir,filename)
    if os.path.isfile(localname):
      return localname  
    else: 
      raise ValueError('No file found')   
    
  def copy_rootfiles(self):
    # cd into tmpdir  
    os.chdir(self.tmpdir)
     
    # copy root files 
    if not os.path.isdir(self.cwdir+'/root-files'):
      os.mkdir(self.cwdir+'/root-files')
    
    for rootfile  in glob.glob('*.root'): 
      basename = os.path.splitext(os.path.basename(rootfile))[0]
      shutil.move(rootfile, self.cwdir+'/root-files/'+basename+'-'+self.name+'.root')  
    

class Simulation(Environment):
  """
  Class to run test beam simulations
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self, steerfiles, name='sim'): 
    Environment.__init__(self, name=name, steerfiles=steerfiles)
    
  def simulate(self, path=['simulation.xml'], caltag=None, ofile='mc.slcio'):
    """
    Creates a lcio file called ofile containing simulated events. 
    :@path:       list containing Marlin xml files that will be executed 
    :@caltag:     name of calibration tag (optional)
    :@ofile:      name of output lcio file
    
    :author: benjamin.schwenker@phys.uni-goettinge.de  
    """  
    
    print ('[INFO] Starting to simulate ' + ofile + ' ...') 

    if not caltag==None:
      self.add_caltag(caltag)
    self.run(path)
    
    src =  os.path.join(self.tmpdir, 'outputfile.slcio')	
    dest = os.path.join(self.cwdir, ofile)	
    shutil.move(src, dest)
   
class Reconstruction(Environment):
  """
  Class to run test beam reconstruction using calibration data
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self, steerfiles, name='reco'): 
    Environment.__init__(self, name=name, steerfiles=steerfiles)
  
  def reconstruct(self, path=['reco.xml'], caltag=None, ifile=None):
    """
    Reconstructs an input file with raw data using a caltag for calibration. 
    :@path:       list containing Marlin xml files that will be executed 
    :@caltag:     name of calibration tag (optional)
    :@ifile:      name of input file with raw data
    
    :author: benjamin.schwenker@phys.uni-goettinge.de  
    """  
    
    if ifile==None:
      raise ValueError('Parameter ifile is missing') 

    if caltag==None:
      raise ValueError('Parameter caltag is missing') 
    
    print ('[INFO] Starting to reconstruct file ' + ifile + ' ...')  
     
    self.add_caltag(caltag)
    self.link_input(ifile)
    self.run(path)
    self.copy_rootfiles()
    
    print ('[INFO] Done processing file ' + ifile)  


class Calibration(Environment):
  """
  Class to run test beam calibration using calibration run
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self, steerfiles, name='cal'): 
    Environment.__init__(self, name=name, steerfiles=steerfiles)
  
  def calibrate(self, path=[], caltag='default', ifile=None):
    """
    Calibrate beam telescope using a calibration run  
    :@path: list containing Marlin xml files that will be executed 
    :@caltag:     name of calibration tag (optional)
    :@ifile:      name of input file with raw data
    
    :author: benjamin.schwenker@phys.uni-goettinge.de  
    """   

    if ifile==None:
      raise ValueError('Parameter ifile is missing') 
    
    print ('[INFO] Starting to calibrate file ' + ifile + ' ...')    
    
    self.link_input(ifile)
    self.run(path)
    self.create_caltag(caltag) 
    
    print ('[INFO] Done processing file ' + ifile)  
     
