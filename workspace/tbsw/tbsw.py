"""
Some helper code to use TBSW with pyhton scripts 

:author: benjamin.schwenker@phys.uni-goettinge.de  
"""

import os
import shutil
import subprocess
import glob
import xml.etree.ElementTree



class Path(object):
  """
  Class which implements an interface to create Marlin steer files. Templates for creating 
  Marling processors taken from steerfiles/processors.xml. Path objects create ready to use
  steerfiles. Template parameters can be locally modified on the fly.
  
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self, name='main', tmpdir=''):
    self.name = name
    self.tmpdir = tmpdir
    
    self.tree = xml.etree.ElementTree.parse(os.path.join(self.tmpdir,'processors.xml'))
    root = self.tree.getroot() 
    
    # remove all processors from tree 
    for processor in root.findall( 'processor' ):
      root.remove(processor)
    
    # add an empty execute tag
    root.append(xml.etree.ElementTree.Element("execute")) 
    
    # write new steer file 
    self.tree.write(os.path.join(self.tmpdir,self.name+'.xml'))
  
  def set_globals(self, params={}):
    global_tag = self.tree.getroot().find("./global") 
    
    # override attributes 
    for key, value in params.items():   
      parameter = global_tag.find("./parameter[@name='%s']" % str(key))
      parameter.set('value', str(value))      
      
    # write new steer file 
    self.tree.write(os.path.join(self.tmpdir,self.name+'.xml'))
    
  def add_processor(self,name='',params={}):
    # try to find processor with name in processor.xml
    basetree = xml.etree.ElementTree.parse(os.path.join(self.tmpdir,'processors.xml'))
    xpath="./processor[@name='%s']" % str(name)
    new_processor = basetree.getroot().find(xpath)
    
    # override attributes on the fly  
    for key, value in params.items():   
      parameter = new_processor.find("./parameter[@name='%s']" % str(key))
      parameter.set('value', str(value)) 
       
    # add processor to xml tree
    root = self.tree.getroot()
    root.append(new_processor)
    
    # add processor to execute tag   
    xml.etree.ElementTree.SubElement(root.find('execute'), tag='processor', attrib={'name' : str(name) })
    
    # write new steer file 
    self.tree.write(os.path.join(self.tmpdir,self.name+'.xml'))
    
  def get_steerfile(self): 
    return self.name+'.xml'
  
  def get_name(self): 
    return self.name
  

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
      raise ValueError('Steerfiles ', steerfiles, ' cannot be found. ', os.getcwd() )
    
    # create localDB folder 
    os.mkdir(self.tmpdir+'/localDB')
     
  def create_path(self, name='main'):
    return Path( name=name, tmpdir=self.tmpdir)
  
  def override(self, xpath='', value=''): 
    xmlfile = self.get_filename('processors.xml') 
    override_xml(xmlfile=xmlfile, xpath=xpath, value=value)
  
  def set_gearfile(self, name=''): 
    self.override(xpath="./global/parameter[@name='GearXMLFile']", value=name)  

  def set_beam_momentum(self, momentum):
    """
    Set momentum in processors.xml file  
    :@momentum:      momentum to be set
    
    :author: benjamin.schwenker@phys.uni-goettinge.de  
    """    
    xmlfile = self.get_filename('processors.xml')
    tree = xml.etree.ElementTree.parse(xmlfile)
    root = tree.getroot()

    print(' Changed beam energy to: ', momentum)
    

    for processor in root.findall('processor'): 
      for parameter in processor.findall('parameter'):
        name = parameter.get('name')
        value = parameter.get('value')
        if name=='ParticleMomentum':
		  parameter.set('value', str(momentum))
        if name=='BeamMomentum':
		  parameter.set('value', str(momentum))

    tree.write(xmlfile) 
    
  def run(self,pathlist):
    # run Marlin in tmpdir  
    os.chdir(self.tmpdir)
    
    for path in pathlist:
      xmlfile = path.get_steerfile()
      logfile = path.get_name() + '.log'
      action = '/$MARLIN/bin/Marlin ' + xmlfile + ' > ' + logfile + ' 2>&1'
      subprocess.call(action, shell=True)
      print ('[INFO] Marlin ' + xmlfile + ' is done')    
    
    # remove tmp* files 
    for tmpfile in glob.glob('tmp*'):
      os.remove(tmpfile)
    
    # go back to workspace
    os.chdir(self.cwdir) 
    
  def import_caltag(self, caltag):  
    # check that caltag exists
    caldir = self.cwdir+'/localDB/'+caltag   
    if os.path.isdir(caldir):  
      for dbfile in glob.glob(caldir+'/*'): 
        shutil.copy(dbfile, os.path.join(self.tmpdir+'/localDB',os.path.basename(dbfile)))  
    else: 
      print ('[INFO] Caltag not found') 
     
  def export_caltag(self, caltag):
    caldir = self.cwdir+'/localDB/'+caltag
    
    # create folder workspace/localDB if not exist    
    if not os.path.isdir(self.cwdir+'/localDB'):
      os.mkdir(self.cwdir+'/localDB')
    # overwrite caltag if exists     
    if os.path.isdir(caldir):
      shutil.rmtree(caldir)
    
    #create caltag and populate with DB files    
    os.mkdir(caldir)		 
    for dbfile in glob.glob(self.tmpdir + '/localDB/*'): 
      shutil.copy(dbfile, os.path.join(caldir,os.path.basename(dbfile)))  
    
    print ('[INFO] Created new caltag ', caltag) 
     	                           

  def copy_file(self, filename, copy_filename):
    localname = os.path.join(self.tmpdir,filename)
    if os.path.isfile(localname):
      shutil.copy(localname, self.tmpdir + '/' + copy_filename ) 
    else: 
      raise ValueError('No file found')  

  def create_dbfilename(self, filename):
    localdbfilename = os.path.join(self.tmpdir + '/localDB/',filename)
    return localdbfilename

  def get_filename(self, filename):
    localname = os.path.join(self.tmpdir,filename)
    if os.path.isfile(localname):
      return localname  
    else: 
      raise ValueError('No file found') 

  def copy_rootfiles(self):
    # create root-files if not exist 
    if not os.path.isdir(self.cwdir+'/root-files'):
      os.mkdir(self.cwdir+'/root-files')
    
    for rootfile  in glob.glob(self.tmpdir + '/*.root'): 
      basename = os.path.splitext(os.path.basename(rootfile))[0]
      shutil.move(rootfile, self.cwdir+'/root-files/'+basename+'-'+self.name+'.root')  

  def get_rootfilename(self,name):

    # if root files director doesn't exist there is nothing to do 
    if not os.path.isdir(self.cwdir+'/root-files'):
      raise ValueError('root-files dir does not exist!') 

    localname=self.cwdir+'/root-files/'+name+'-'+self.name+'.root'

    if os.path.isfile(localname):
      return localname  
    else: 
      raise ValueError('No file found')    

class Simulation(Environment):
  """
  Class to run test beam simulations
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self, steerfiles, name='sim'): 
    Environment.__init__(self, name=name, steerfiles=steerfiles)
    
  def simulate(self, path=[], caltag='default'):
    """
    Creates a lcio file(s) with simulated events. 
    :@path:       list of path objects 
    :@caltag:     name of calibration tag (optional)
    
    :author: benjamin.schwenker@phys.uni-goettinge.de  
    """  
    
    print ('[INFO] Starting to simulate a test beam run ...') 
     
    self.import_caltag(caltag)
    self.run(path)
    self.export_caltag(caltag)
   
class Reconstruction(Environment):
  """
  Class to run test beam reconstruction using calibration data
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self, steerfiles, name='reco'): 
    Environment.__init__(self, name=name, steerfiles=steerfiles)
  
  def reconstruct(self, path=[], caltag='default', ifile=''):
    """
    Reconstructs an input file with raw data using a caltag for calibration. 
    :@path:       list of path objects
    :@caltag:     name of calibration tag (optional)
    :@ifile:      name of input file with raw data
    
    :author: benjamin.schwenker@phys.uni-goettinge.de  
    """  
    
    print ('[INFO] Starting to reconstruct file ' + ifile + ' ...')  
     
    self.import_caltag(caltag)
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
  
  def calibrate(self, path=[], caltag='default', ifile=''):
    """
    Calibrate beam telescope using a calibration run  
    :@path: list containing Marlin xml files that will be executed 
    :@caltag:     name of calibration tag (optional)
    :@ifile:      name of input file with raw data
    
    :author: benjamin.schwenker@phys.uni-goettinge.de  
    """   
      
    print ('[INFO] Starting to calibrate file ' + ifile + ' ...')    
    
    self.run(path)
    self.export_caltag(caltag) 
    
    print ('[INFO] Done processing file ' + ifile)  


def override_xml(xmlfile='', xpath='', value=None):
  """
  Overrides value field in all nodes specified by xpath found in xmlfile
    :@xmlfile:    XML steering file to be overwritten  
    :@xpath:      xpath to list of XML nodes to be modified    
    :@value:      insert string value in field value 
    :author: benjamin.schwenker@phys.uni-goettinge.de  
  """   
  tree = xml.etree.ElementTree.parse(xmlfile)
  root = tree.getroot() 
  
  for node in root.findall(xpath): 
    node.set('value', str(value))   
  tree.write(xmlfile)




     
