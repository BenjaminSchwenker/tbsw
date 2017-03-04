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
  caltag=''

  try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:c:",["ifile=","imagefile=","scriptsfolder=","caltag="])
  except getopt.GetoptError:
    print ('GenerateImage.py -i <inputfile> -c <cal-tag>' )
    sys.exit(2)

  for opt, arg in opts:
    if opt == '-h':
      print ('GenerateImage.py -i <inputfile> -c <cal-tag>')
      sys.exit()
    elif opt in ("-i", "--ifile"):
      rootfile = arg
    elif opt in ("-c", "--caltag"):
      caltag = arg

  if rootfile == '':
    print ('missing option: -i path/to/inputfilename.root')
    sys.exit(2)   

  # remember current working dir 
  fullpath = os.getcwd() 
  
  # Find directory with scripts
  scriptsfolder = fullpath+'/root-scripts/x0imaging'

  # get X0 filename filename without path
  this_x0filename = os.path.splitext(os.path.basename(rootfile))[0]

  # copy rootfile to current work directory
  shutil.copy(rootfile, 'X0.root')

  # Copy cfg file
  cfgname="x0calibration.cfg"
  shutil.copy(scriptsfolder+'/'+cfgname, cfgname)


  # Directory where the results will be stored and where the image file should be
  resdir=fullpath+'/root-files/'+caltag
  imagefile=resdir+'/X0-completeimage.root'

  if os.path.isfile(imagefile):

	#get X0 filename filename without path
  	this_imagefilename = os.path.splitext(os.path.basename(imagefile))[0]
	
	# Copy X0image file to work directory
  	shutil.copy(imagefile, 'X0image.root')

  	# Copy placeholder x0image and change it
  	scriptname="DrawBoxes.C"
  	shutil.copy(scriptsfolder+'/DrawBoxes.C', scriptname)
  
  	subprocess.call('root -q '+scriptname, shell=True)
  	print ('[Print] Marking of measurement areas done... ')

  	# remove DrawBoxes.C script and image file
  	os.remove(scriptname) 
  	os.remove('X0image.root') 

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
   
  subprocess.call('root -x '+scriptname+' | tee log-x0cal.txt', shell=True)            
  print ('[Print] Calibration done... ')  

  # remove calibrationfit.C script
  os.remove(scriptname) 

  # remove cfg file
  os.remove(cfgname) 
  os.remove('X0.root') 

  caldir=fullpath+'/cal-files/'+caltag

  if os.path.isdir(caldir):
	
    print ('[Print] Copy cfg file ') 

    # save all interesting files to common folder   
    for cfgfile in glob.glob('*.cfg'): 
      shutil.copy(cfgfile, os.path.join(caldir,cfgfile))  

  else :
    print ('[Print] Cal tag folder does not exist! cfg file will not be copied. ') 

  if os.path.isdir(resdir):
	
    # save all interesting files to common folder   
    for pdffile in glob.glob('*.pdf'): 
      shutil.copy(pdffile, os.path.join(resdir,pdffile))  

    for rootfile in glob.glob('*.root'): 
      shutil.copy(rootfile, os.path.join(resdir,rootfile))  

    for txtfile in glob.glob('*.txt'): 
      shutil.copy(txtfile, os.path.join(resdir,txtfile)) 

  else :
    print ('[Print] Results folder does not exist! cfg file will not be copied. ') 


  # Remove files
  for file in glob.glob('*.pdf'):
    os.remove(file) 
  for file in glob.glob('*.root'):
    os.remove(file) 
  for file in glob.glob('*.cfg'):
    os.remove(file) 
  for file in glob.glob('*.cfg.bak'):
    os.remove(file) 
  for file in glob.glob('*.txt'):
    os.remove(file) 
               
		           
                


 	


