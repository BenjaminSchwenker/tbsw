import os
import shutil
import subprocess
import sys
import glob
import math
import fileinput
import re

# Function which merges the result root files
def merge_rootfile(filename=None,RunList='',caltag=None):
  """
  Merges all root files from a list
    :@filename:    Name of merged output rootfile 
    :@RunList:     List of filenames , which will be merged 
    :@caltag:      caltag string needed to convert the list strings  
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """ 

  if filename == None:
    return None

  if caltag == None:
    return None

  flags='hadd '+filename+' '
  for run in RunList:
    name=os.path.splitext(os.path.basename(run))[0]
    if caltag == '':
      flags=flags+'root-files/X0-'+name+'.root '
    else:
      flags=flags+'root-files/X0-'+name+'-'+caltag+'.root '


  if os.path.isfile(filename):
    os.remove(filename)
  subprocess.call(flags, shell=True)

  return None

# Create List with filenames of root tree files from raw file list
def CreateRootFileList(rawlist,rootlist=[], caltag=''):
  """
  Creates a list with filenames of root tree files from a similar raw file list
    :@rawlist:    List of raw file names
    :@rootlist:   List of root file names to be created  
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """ 

  for rawfile in rawlist:
    name=os.path.splitext(os.path.basename(rawfile))[0]
    if caltag=='':
      rootlist.append('root-files/X0-'+name+'.root')
    else:
      rootlist.append('root-files/X0-'+name+'-'+caltag+'.root')


def x0calibration(rootfilelist=[],imagefile='',caltag='',steerfiles=''):
  """
  Performs radiation length calibration on well known material target
    :@filename:    List of input root file with calibration data 
    :@caltag:      Use x0 calibration parameters from cfg files in the caltag subdirectory and save results there    
    :@steerfiles:  Directory, where the image cfg file lies
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """   
  
  cfgfile=steerfiles+'x0.cfg'

  # Remember current working dir 
  fullpath = os.getcwd() 

  # Path to work directory
  workdir = 'tmp-runs/X0Calibration-'+caltag

  # remove old workdir if exists 
  if os.path.isdir(workdir):
    shutil.rmtree(workdir)
    
  # create workdir and change directory
  os.mkdir(workdir) 
  os.chdir(workdir)
  
  # Create Links of input root files
  for ifile,rootfile in enumerate(rootfilelist):
    this_x0filename = os.path.splitext(os.path.basename(rootfile))[0]
    
    runtags=re.findall('run\d+', this_x0filename )
    if len(runtags)>1:
      print('Input file ' + str(ifile+1) + ' contains two or more possible runnumbers! Only the first possibility is as link name.')
      if os.path.isfile(runtags[0]):
        print("Error link name already exists! ")
        sys.exit(2)
      else:
        os.symlink(fullpath+'/'+rootfile,runtags[0]) 

    elif len(runtags)==1:
      if os.path.isfile(runtags[0]):
        print("Error link name already exists! ")
      else:
        os.symlink(fullpath+'/'+rootfile,runtags[0]) 

    else:
      print('Input file ' + str(ifile+1) + ' contains no string with a runnumber! Creating link with default name.')
      runtag='run'+str(ifile+1)
      if os.path.isfile(runtag):
        print("Error link name already exists! ")
      else:
        os.symlink(fullpath+'/'+rootfile,runtag) 

  # Copy cfg file to current work dir
  cfgname="x0.cfg"
  shutil.copy(fullpath+'/'+cfgfile, cfgname)
    
  if os.path.isfile(fullpath+'/'+imagefile):
    #get X0 filename filename without path
    this_imagefilename = os.path.splitext(os.path.basename(imagefile))[0]
    
    # Create X0image root file link in the current work dir
    os.symlink(fullpath+'/'+imagefile,'X0image')
     
    subprocess.call( "$X0TOOLS/bin/DrawBoxes > x0-drawboxes.log  2>&1", shell=True)   
    print ('[INFO] Marking of measurement areas done... ')
  else:
    print ('[Print] No Image file found... Skip marking measurement areas! ')
  
  # Copy the results cfg file from previous calibrations, if it exists
  cfgfilename="x0cal_result.cfg"
  cfgfile=fullpath+'/localDB/'+caltag+'/'+cfgfilename
  if os.path.isfile(cfgfile):
    shutil.copy(cfgfile, cfgfilename)
  
  print ('[INFO] Start X0 Calibration fit ... ')  
  subprocess.call( "$X0TOOLS/bin/calibrationfit > x0cal.log  2>&1", shell=True)            
  print ('[INFO] X0 Calibration fit done! ')  
  
  caldir=fullpath+'/localDB/'+caltag
  
  if not os.path.isdir(caldir):
    os.mkdir(caldir)
    
  print ('[INFO] Copy cfg file ') 
  
  # save all interesting files to common folder   
  for cfgfile in glob.glob('*.cfg'): 
    shutil.copy(cfgfile, os.path.join(caldir,cfgfile))  
  
  # Remove files
  for file in glob.glob('*.cfg'):
    os.remove(file) 
  for file in glob.glob('*.cfg.bak'):
    os.remove(file) 
  
  # Go back to workspace
  os.chdir(fullpath) 

