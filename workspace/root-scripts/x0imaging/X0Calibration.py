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
  imagefile = ''
  scriptsfolder = ''

  try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:m:s:",["ifile=","imagefile=","scriptsfolder="])
  except getopt.GetoptError:
    print ('GenerateImage.py -i <inputfile>  -m <imagefile> -s <scriptsfolder>' )
    sys.exit(2)

  for opt, arg in opts:
    if opt == '-h':
      print ('GenerateImage.py -i <inputfile> -m <imagefile> -s <scriptsfolder>')
      sys.exit()
    elif opt in ("-i", "--ifile"):
      rootfile = arg
    elif opt in ("-m", "--imagefile"):
      imagefile = arg
    elif opt in ("-s", "--scriptsfolder"):
      scriptsfolder = arg

  if rootfile == '':
    print ('missing option: -i path/to/inputfilename.root')
    sys.exit(2)  

  # remember current working dir 
  fullpath = os.getcwd() 

  # get X0 filename filename without path
  this_x0filename = os.path.splitext(os.path.basename(rootfile))[0]

  # copy rootfile to current work directory
  shutil.copy(rootfile, 'X0.root')

  # Copy cfg file
  cfgname="x0calibration.cfg"
  shutil.copy(scriptsfolder+'/'+cfgname, cfgname)

  if imagefile != '':

	#get X0 filename filename without path
  	this_imagefilename = os.path.splitext(os.path.basename(imagefile))[0]
	
	# Copy X0image file to work directory
  	shutil.copy(imagefile, 'X0image.root')

  	# Copy placeholder x0image and change it
  	scriptname="DrawBoxes.C"
  	shutil.copy(scriptsfolder+'/DrawBoxes.C', scriptname)
  
  	subprocess.call('root -q '+scriptname, shell=True)

  	# remove DrawBoxes.C script and image file
  	os.remove(scriptname) 
  	os.remove('X0image.root') 
  
  # Copy the calibration fit script
  scriptname="calibrationfit.C"
  shutil.copy(scriptsfolder+'/'+scriptname, scriptname)
   
  subprocess.call('root -x '+scriptname, shell=True)            
  print ('[Print] All partial images created and merged... ')  

  # remove calibrationfit.C script
  os.remove(scriptname) 

  # remove cfg file
  os.remove(cfgname) 
  os.remove('X0.root') 
               
		           
                


 	


