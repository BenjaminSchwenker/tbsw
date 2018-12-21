"""
Some helper code to use Marlin with pyhton scripts 

:author: benjamin.schwenker@phys.uni-goettinge.de  
"""

import copy
import os
import shutil
import subprocess
import glob
import xml.etree.ElementTree


class Processor(object):
  """
  Processor class implements an interface to create a <processor> node in a Marlin steer file. 
  It is basically a simple container for the processor's name and type as well as all non default
  parameter values.  
   
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self, name, proctype, params={}):
    """  
    :@name: processor name string
    :@proctype: processor type string
    :@params: parameter dictionary
    """    
    self.params = {}
    self.params.update(params)
    
    self.name = name
    self.proctype = proctype

  def param(self, name, value):
    """  
    :@name: parameter name string
    :@value: parameter value (must be convertible to string) 
    """    
    self.params[name] = value
  
  def params(self, params):
    """  
    :@params: parameter dictionary
    """   
    self.params.update(params)
  
  def __str__(self):
    return "Processor name={} type={} override={}".format(self.name, self.proctype,str(self.params))
  
class Path(object):
  """
  Path class implements an interface to create Marlin steer files. The command Marlin -x creates  
  a xml file containing all available processor types and their steering parameters with default values.
  Default values for processor parameters can be overriden.  
   
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self, name='main', tmpdir=''):
    """  
    :@name: name of path 
    :@tmpdir: name of folder to put the steerfile  
    """   
    self.name = name
    self.tmpdir = tmpdir
    self.globals = {}
    self.processors = []      
  
  def set_globals(self, params):
    """  
    :@params: dictioanary for globals parameters of steerfile 
    """   
    self.globals.update(params)
      
  def add_processor(self, processor):
    """  
    :@processor: add processor to path
    """    
    self.processors.append(processor)    

  def get_steerfile(self): 
    """  
    Create xml steerfile for Marlin 
    """      
     
    # Create a xml tree from the template file processors.xml 
    # Never change the variable basetree 
    basetree = xml.etree.ElementTree.parse(os.path.join(self.tmpdir,'processors.xml'))
    
    # Create a xml tree from the template file processors.xml 
    # We will create the new steerfile from the variable tree
    tree = xml.etree.ElementTree.parse(os.path.join(self.tmpdir,'processors.xml'))
    root = tree.getroot() 
    
    # remove all processors from tree 
    for processor in root.findall( 'processor' ):
      root.remove(processor)
    
    # add an empty execute tag
    root.append(xml.etree.ElementTree.Element("execute")) 
    
    # override global parameters as needed
    global_node = root.find("./global")  
    for key, value in self.globals.items():   
      parameter = global_node.find("./parameter[@name='%s']" % str(key))
      parameter.set('value', str(value))      
        
    # create processor nodes   
    for processor in self.processors:
      # find processor node in basetree 
      xpath="./processor[@type='%s']" % str(processor.proctype)
      processor_node_base = basetree.getroot().find(xpath)
      
      # make a copy to insert into your new tree
      processor_node = copy.deepcopy(processor_node_base)      
      
      # set the processors name
      processor_node.set('name',str(processor.name))
      
      # find and remove all optional parameters 
      for parameter in processor_node.findall( "./parameter[@isOptional]" ):
        processor_node.remove(parameter)      
      
      # override processor parameter as needed 
      for key, value in processor.params.items():   
        parameter = processor_node.find("./parameter[@name='%s']" % str(key))
        parameter.set('value', str(value)) 
       
      # add processor to xml tree
      root.append(processor_node)
      
      # add processor to execute tag   
      xml.etree.ElementTree.SubElement(root.find('execute'), tag='processor', attrib={'name' : str(processor.name) })
      
    # write new steer file 
    tree.write(os.path.join(self.tmpdir,self.name+'.xml'))
    
    # return name of steering file 
    return self.name+'.xml'
    
  def get_name(self): 
    return self.name
   
class Environment(object):
  """
  Class which implements an environment for executing Marlin with all needed 
  steering files. 
  
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
      raise ValueError('Steerfiles {} cannot be found.'.format(steerfiles) )
    
    # create file processors.xml 
    subprocess.call('/$MARLIN/bin/Marlin -x > {}'.format(self.tmpdir+'/processors.xml'), shell=True)
    
    # create localDB folder 
    os.mkdir(self.tmpdir+'/localDB')
    
  def create_path(self, name='main'):
    return Path( name=name, tmpdir=self.tmpdir)
  
  def override(self, xpath='', value=''): 
    xmlfile = self.get_filename('processors.xml') 
    override_xml(xmlfile=xmlfile, xpath=xpath, value=value)
  
  def set_gearfile(self, name=''): 
    self.override(xpath="./global/parameter[@name='GearXMLFile']", value=name)  
  
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

    return self.tmpdir + '/' + copy_filename

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
    
  def simulate(self, paths=[], caltag='default'):
    """
    Creates a lcio file(s) with simulated events. 
    :@paths:      list of path objects 
    :@caltag:     name of calibration tag (optional)
    
    :author: benjamin.schwenker@phys.uni-goettinge.de  
    """  
    
    print ('[INFO] Starting to simulate a test beam run ...') 
     
    self.import_caltag(caltag)
    self.run(paths)
    self.export_caltag(caltag)
   
class Reconstruction(Environment):
  """
  Class to run test beam reconstruction using calibration data
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self, steerfiles, name='reco'): 
    Environment.__init__(self, name=name, steerfiles=steerfiles)
  
  def reconstruct(self, paths=[], caltag='default', ifile=''):
    """
    Reconstructs an input file with raw data using a caltag for calibration. 
    :@paths:      list of path objects
    :@caltag:     name of calibration tag (optional)
    :@ifile:      name of input file with raw data
    
    :author: benjamin.schwenker@phys.uni-goettinge.de  
    """  
    
    print ('[INFO] Starting to reconstruct file ' + ifile + ' ...')  
     
    self.import_caltag(caltag)
    self.run(paths)
    self.copy_rootfiles()
    
    print ('[INFO] Done processing file ' + ifile)  


class Calibration(Environment):
  """
  Class to run test beam calibration using calibration run
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self, steerfiles, name='cal'): 
    Environment.__init__(self, name=name, steerfiles=steerfiles)
  
  def calibrate(self, paths=[], caltag='default', ifile=''):
    """
    Calibrate beam telescope using a calibration run  
    :@paths:      list containing Marlin xml files that will be executed 
    :@caltag:     name of calibration tag (optional)
    :@ifile:      name of input file with raw data
    
    :author: benjamin.schwenker@phys.uni-goettinge.de  
    """   
      
    print ('[INFO] Starting to calibrate file ' + ifile + ' ...')    
    
    self.run(paths)
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




     
