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
  scriptsfolder = ''
  deletetag=`1`

  try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:s:d:",["ifile=","scriptsfolder=","deletetag="])
  except getopt.GetoptError:
    print ('GenerateImage.py -i <inputfile> -s <scriptsfolder> -d <deletetag>')
    print ('-d is optional and defaults to: ' + deletetag )
    print (' 1: delete partial rootfiles' )
    print (' 0: dont delete partial rootfiles' )
    sys.exit(2)

  for opt, arg in opts:
    if opt == '-h':
      print ('GenerateImage.py -i <inputfile> -d <deletetag>')
      print ('-d is optional and defaults to: ' + deletetag )
      sys.exit()
    elif opt in ("-i", "--ifile"):
      rootfile = arg
    elif opt in ("-s", "--scriptsfolder"):
      scriptsfolder = arg
    elif opt in ("-d", "--deletetag"):
      deletetag = arg

  if rootfile == '':
    print ('missing option: -i path/to/inputfilename.root')
    sys.exit(2)  

  # remember current working dir 
  fullpath = os.getcwd() 

  # get X0 filename without path
  this_x0filename = os.path.splitext(os.path.basename(rootfile))[0]

  # copy rootfile to current work directory
  shutil.copy(rootfile, this_x0filename+'.root')

  #-----------------------------------------------------------------
  #-------------------- Change Parameters here ---------------------
  #-----------------------------------------------------------------

  # u and v size of complete X0 image 
  u_length=20.0  #mm
  v_length=10.0  #mm

  # umin and vmax of the complete X0 image
  umin=-10.0
  vmax=+5.0

  # Pixel sizes of the image
  u_pixel_size=200.0 #mum
  v_pixel_size=200.0 #mum

  # number of pixels needed
  num_u_pixels=u_length * 1000.0 / u_pixel_size
  num_v_pixels=v_length * 1000.0 / v_pixel_size
 
  # Test
  print num_u_pixels
  print num_v_pixels

  # Max number of pixels per x0image (these numbers ensure a stable imaging)
  max_u_pixels=100
  max_v_pixels=50

  # How many x0image parts are required in each direction?
  u_splits=math.ceil(num_u_pixels/max_u_pixels)
  v_splits=math.ceil(num_v_pixels/max_v_pixels)

  # Calibration factor for this telescope setup
  calibrationfactor=1.0

  # Mean value of the momentum of the beam particles
  momentummean=4.0

  # Mean value of the momentum of the beam particles
  momentumslope=0.0

  # Name of the results file
  this_resultsfilename="X0-merge-completearea-image"

  #-----------------------------------------------------------------
  #---------------------- End of Parameters ------------------------
  #-----------------------------------------------------------------

  # Define the placeholder strings, which will be replaced in the x0image scripts
  placeholder_functionname="PARAMETER-FUNCTIONNAME"
  placeholder_x0filename="PARAMETER-X0FILENAME"
  placeholder_x0fileidentifier="PARAMETER-FILEINDENTIFIER"
  placeholder_umin="PARAMETER-UMIN"
  placeholder_vmax="PARAMETER-VMAX"
  placeholder_ulength="PARAMETER-ULENGTH"
  placeholder_vlength="PARAMETER-VLENGTH"
  placeholder_calibrationfactor="PARAMETER-CALIBRATIONFACTOR"
  placeholder_momentummean="PARAMETER-MOMENTUMMEAN"
  placeholder_momentumslope="PARAMETER-MOMENTUMSLOPE"


  # Define the placeholder strings, which will be replaced in the X0 merge script
  placeholder_usplits="PARAMETER-USPLITS"
  placeholder_vsplits="PARAMETER-VSPLITS"
  placeholder_upixelsize="PARAMETER-UPIXELSIZE"
  placeholder_vpixelsize="PARAMETER-VPIXELSIZE"
  placeholder_maxupixels="PARAMETER-MAXUPIXELS"
  placeholder_maxvpixels="PARAMETER-MAXVPIXELS"
  placeholder_resultsfilename="PARAMETER-RESULTSFILENAME"

  # Overall number of X0 image parts required
  num_parts=u_splits * v_splits

  # Test
  print num_parts

  # copy x0image files 
  for i in range(0,int(u_splits)):
    for j in range(0,int(v_splits)):
      scriptname="x0imaging_part_"+`i+1`+"_"+`j+1`+".C"
      # Function name for this x0image script
      this_functionname="x0imaging_part_"+`i+1`+"_"+`j+1`
      # fileidentifier for this x0image
      this_x0fileidentifier="-part_"+`i+1`+"_"+`j+1`
      # umin for this x0image
      this_umin=`umin + i * max_u_pixels * u_pixel_size / 1000.0`
      # vmax for this x0image
      this_vmax=`vmax-j * max_v_pixels * v_pixel_size / 1000.0`
      # ulength for this x0image
      this_ulength=`max_u_pixels * u_pixel_size / 1000.0`
      # vlength for this x0image
      this_vlength=`max_v_pixels * v_pixel_size / 1000.0`
      # calibrationfactor
      this_calibrationfactor=`calibrationfactor`
      # mean momentum
      this_momentummean=`momentummean`
      # momentum slope
      this_momentumslope=`momentumslope`
      
      # Copy placeholder x0image and change it
      shutil.copy(scriptsfolder+'/x0imaging_mask.C', scriptname)

      # Open the x0imaging script of this x0image part
      tempFile = None

      with open(scriptname, 'r') as file:
	tempFile = file.read()

      # Replace the placeholder strings in the x0image_mask.C script
      tempFile = tempFile.replace(placeholder_functionname, this_functionname)
      tempFile = tempFile.replace(placeholder_x0filename, this_x0filename)
      tempFile = tempFile.replace(placeholder_x0fileidentifier, this_x0fileidentifier)
      tempFile = tempFile.replace(placeholder_umin, this_umin)
      tempFile = tempFile.replace(placeholder_vmax, this_vmax)
      tempFile = tempFile.replace(placeholder_ulength, this_ulength)
      tempFile = tempFile.replace(placeholder_vlength, this_vlength)
      tempFile = tempFile.replace(placeholder_calibrationfactor, this_calibrationfactor)
      tempFile = tempFile.replace(placeholder_momentummean, this_momentummean)
      tempFile = tempFile.replace(placeholder_momentumslope, this_momentumslope)


      # Write the script out again
      with open(scriptname, 'w') as file:
  	file.write(tempFile)
  
      subprocess.call('root -q '+scriptname, shell=True)
      print ('[Print] One part done ...')

  # Modify the merging script
  scriptname="MergeImages.C"

  # Copy placeholder x0image and change it
  shutil.copy(scriptsfolder+'/MergeImages_mask.C', scriptname)

  # Open the x0imaging script of this x0image part
  tempFile = None

  with open(scriptname, 'r') as file:
    tempFile = file.read()

  # Create strings, which will replace the placeholder strings
  this_usplits=`u_splits`
  this_vsplits=`v_splits`
  this_upixelsize=`u_pixel_size`
  this_vpixelsize=`v_pixel_size`
  this_maxupixels=`max_u_pixels`
  this_maxvpixels=`max_v_pixels`
  this_umin=`umin`
  this_vmax=`vmax`
  this_ulength=`u_length`
  this_vlength=`v_length`

  # Replace the placeholder strings

  tempFile = tempFile.replace(placeholder_x0filename, this_x0filename)
  tempFile = tempFile.replace(placeholder_resultsfilename, this_resultsfilename)
  tempFile = tempFile.replace(placeholder_usplits, this_usplits)
  tempFile = tempFile.replace(placeholder_vsplits, this_vsplits)
  tempFile = tempFile.replace(placeholder_upixelsize, this_upixelsize)
  tempFile = tempFile.replace(placeholder_vpixelsize, this_vpixelsize)
  tempFile = tempFile.replace(placeholder_maxupixels, this_maxupixels)
  tempFile = tempFile.replace(placeholder_maxvpixels, this_maxvpixels)
  tempFile = tempFile.replace(placeholder_umin, this_umin)
  tempFile = tempFile.replace(placeholder_vmax, this_vmax)
  tempFile = tempFile.replace(placeholder_ulength, this_ulength)
  tempFile = tempFile.replace(placeholder_vlength, this_vlength)

  # Write the script out again
  with open(scriptname, 'w') as file:
    file.write(tempFile)
   
  subprocess.call('root -q '+scriptname, shell=True)            
  print ('[Print] All partial images created and merged... ')  

  # clean up partial image scripts 
  for tmpfile in glob.glob('*part_*.C'):
    os.remove(tmpfile) 

  if deletetag == `1`:
    # clean up root files
    for tmpfile in glob.glob('*part_*.root'):
      os.remove(tmpfile) 

  # remove MergeImage.C script
  os.remove('MergeImages.C') 
               
		           
                


 	


