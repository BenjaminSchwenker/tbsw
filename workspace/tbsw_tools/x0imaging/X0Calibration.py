#!/usr/bin/python3

import os
import shutil
import subprocess
import sys, getopt 
import glob
import math
import fileinput

if __name__ == '__main__':
  
  rootfile = ''
  cfgfile = ''
  imagefile = ''
  caltag=''

  try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:f:m:c:",["ifile=","cfgfile=","mfile=","caltag="])
  except getopt.GetoptError:
    print ('X0Calibration.py -i <inputfile> -f <cfgfile> -m <imagefile> -c <cal-tag>' )
    sys.exit(2)

  for opt, arg in opts:
    if opt == '-h':
      print ('X0Calibration.py -i <inputfile> -f <cfgfile> -m <imagefile> -c <cal-tag>')
      sys.exit()
    elif opt in ("-i", "--ifile"):
      rootfile = arg
    elif opt in ("-f", "--cfgfile"):
      cfgfile = arg
    elif opt in ("-m", "--mfile"):
      imagefile = arg
    elif opt in ("-c", "--caltag"):
      caltag = arg

  if rootfile == '':
    print ('missing option: -i path/to/inputfilename.root')
    sys.exit(2) 

  if cfgfile == '':
    print ('missing option: -i path/to/cfgfilename.cfg')
    sys.exit(2)   

  # remember current working dir 
  fullpath = os.getcwd() 

  # Directory where the results will be stored and where the image file should be
  resdir=os.path.splitext(os.path.dirname(rootfile))[0]

  # get X0 filename without path
  this_x0filename = os.path.splitext(os.path.basename(rootfile))[0]

  # Change to work directory
  workdir = 'tmp-runs/'+this_x0filename+'-X0Calibration'

  # remove old workdir if exists 
  if os.path.isdir(workdir):
    shutil.rmtree(workdir)
    
  # create workdir and change directory
  os.mkdir(workdir) 
  os.chdir(workdir)
  
  # Find directory with scripts
  scriptsfolder = fullpath+'/tbsw_tools/x0imaging'

  # Create X0 root file link in the current work dir
  os.symlink(fullpath+'/'+rootfile,'X0File')

  # Copy cfg file to current work dir
  cfgname="x0calibration.cfg"
  shutil.copy(fullpath+'/'+cfgfile, cfgname)
    
  if os.path.isfile(fullpath+'/'+imagefile):

    #get X0 filename filename without path
    this_imagefilename = os.path.splitext(os.path.basename(imagefile))[0]
    
    # Create X0image root file link in the current work dir
    os.symlink(fullpath+'/'+imagefile,'X0image')

    # Copy placeholder x0image and change it
    scriptname="DrawBoxes.C"
    shutil.copy(scriptsfolder+'/DrawBoxes.C', scriptname)
  
    subprocess.call('root -q -b '+scriptname, shell=True)
    print ('[Print] Marking of measurement areas done... ')

    # remove DrawBoxes.C script and image file
    os.remove(scriptname) 

  else:
    print ('[Print] No Image file found... Skip marking measurement areas! ')

  
  # Copy the calibration fit script
  scriptname="calibrationfit.C"
  shutil.copy(scriptsfolder+'/'+scriptname, scriptname)

  # Copy the results cfg file from previous calibrations, of it exists
  cfgfilename="x0cal_result.cfg"
  cfgfile=fullpath+'/cal-files/'+caltag+'/'+cfgfilename
  if os.path.isfile(cfgfile):
    shutil.copy(cfgfile, cfgfilename)
   
  subprocess.call('root -x -b '+scriptname+' | tee log-x0cal.txt', shell=True)            
  print ('[Print] Calibration done... ')  

  # remove calibrationfit.C script
  os.remove(scriptname) 

  # remove cfg file
  os.remove(cfgname) 

  caldir=fullpath+'/cal-files/'+caltag

  if not os.path.isdir(caldir):
    os.mkdir(caldir)
    
  print ('[Print] Copy cfg file ') 

  # save all interesting files to common folder   
  for cfgfile in glob.glob('*.cfg'): 
    shutil.copy(cfgfile, os.path.join(caldir,cfgfile))  

  # Remove files
  for file in glob.glob('*.cfg'):
    os.remove(file) 
  for file in glob.glob('*.cfg.bak'):
    os.remove(file) 

                   
                


     


