#!/usr/bin/python3

import sys, getopt 
import os
import shutil
import subprocess
import glob


if __name__ == '__main__':
  
  rawfile = ''
  xmlfile = ''
  caltag='dummy'
    
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:x:c:",["ifile=","xmlfile=","caltag="])
  except getopt.GetoptError:
    print ('reco.py -i <inputfiles>  -x <xmlfile> -c <caltag>')
    print ('-c is optional and defaults to: ' + caltag )
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print ('reco.py -i <inputfiles> -x <xmlfile> -c <caltag>')
      print ('-c is optional and defaults to: ' + caltag )
      sys.exit()
    elif opt in ("-i", "--ifile"):
      rawfile = arg
    elif opt in ("-x", "--xmlfile"):
      xmlfile = arg
    elif opt in ("-c", "--caltag"):
      caltag = arg
  
  if rawfile == '':
    print ('missing option: -i path/to/inputfilename')
    sys.exit(2)  
  
  if xmlfile == '':
    print ('missing option: -x path/to/steeringfilename.xml')
    sys.exit(2)  
  
  # remember current working dir 
  fullpath = os.getcwd() 
   
  store_lcio = False
  
  # get runtag
  runtag = os.path.splitext(os.path.basename(rawfile))[0]
                    
  print ('[Print] Starting to process file ' + runtag + ' ...')
  
  # create tmp dir
  if not os.path.isdir(fullpath+'/tmp-runs'):
    os.mkdir(fullpath+'/tmp-runs')
  tmpdir = os.path.join(fullpath+'/tmp-runs',runtag)  
  tmpdir = tmpdir + '-' + caltag + '-reco' 
  
  if os.path.isdir(tmpdir):
    shutil.rmtree(tmpdir)
       
  os.mkdir(tmpdir)
 
  # create symlink to raw file 
  os.symlink( os.path.join(fullpath, rawfile), tmpdir+'/inputfilename')
  # copy steering file 
  shutil.copy(xmlfile, tmpdir+'/reco.xml')
  
  # check that calibration files exist
  caldir = caldir = fullpath+'/cal-files/'+caltag   
  if not os.path.isdir(caldir):
    print ('[Print] warning: no calibration data found')
    return  
  else: 
    # copy calibration files 
    shutil.copytree(caldir,tmpdir+'/cal-files')
  
  # run reco in tmp dir 
  os.chdir(tmpdir)
                       
  subprocess.call('/$MARLIN/bin/Marlin reco.xml > log.txt 2>&1', shell=True)
  
  # clean up inputfile
  os.remove('inputfilename')

  # clean up tmp files 
  for tmpfile in glob.glob('tmp*'):
    os.remove(tmpfile)

  # store dqm files 
  if not os.path.isdir(fullpath+'/root-files'):
    os.mkdir(fullpath+'/root-files')
  
  for dqmfile in glob.glob('*.root'): 
    name = os.path.splitext(os.path.basename(dqmfile))[0]
    shutil.move(dqmfile, fullpath+'/root-files/'+name+'-'+runtag+'-'+caltag+'.root')  
  
  os.chdir(fullpath)

  print ('[Print] Processing file ' + runtag + ' done!')

 
