#!/usr/bin/python3

import os
import shutil
import subprocess
import sys, getopt 
import glob
import math
import fileinput
from configparser import ConfigParser

if __name__ == '__main__':
  
  rootfile = ''
  cfgfile = ''
  caltag=''
  tag='Uncalibrated'
  deletetag='0'

  try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:f:c:d:t:",["ifile=","cfgfile=","caltag=","deletetag=","tag="])
  except getopt.GetoptError:
    print ('GenerateImage.py -i <inputfile> -f <cfgfile> -c <cal-tag> -d <deletetag>')
    print ('-d is optional and defaults to: ' + deletetag )
    print (' 1: delete partial rootfiles' )
    print (' 0: dont delete partial rootfiles' )
    print ('-t is optional and defaults to: ' + tag )
    sys.exit(2)

  for opt, arg in opts:
    if opt == '-h':
      print ('GenerateImage.py -i <inputfile> -f <cfgfile> -c <cal-tag> -d <deletetag> -t <tag>')
      print ('-d is optional and defaults to: ' + deletetag )
      sys.exit()
    elif opt in ("-i", "--ifile"):
      rootfile = arg
    elif opt in ("-f", "--cfgfile"):
      cfgfile = arg
    elif opt in ("-c", "--caltag"):
      caltag = arg
    elif opt in ("-d", "--deletetag"):
      deletetag = arg
    elif opt in ("-t", "--tag"):
      tag = arg

  if rootfile == '':
    print ('missing option: -i path/to/inputfilename.root')
    sys.exit(2)  

  if cfgfile == '':
    print ('missing option: -i path/to/cfgfile.cfg')
    sys.exit(2)  

  # remember current working dir 
  fullpath = os.getcwd() 

  # get X0 filename without path
  this_x0filename = os.path.splitext(os.path.basename(rootfile))[0]

  # Get the path without the X0 filename
  this_x0filepath = os.path.splitext(os.path.dirname(rootfile))[0]

  # Change to work directory
  workdir = 'tmp-runs/'+this_x0filename+'-'+ tag +'-X0image'

  # remove old workdir if exists 
  if os.path.isdir(workdir):
    shutil.rmtree(workdir)
    
  # create workdir and change directory
  os.mkdir(workdir) 
  os.chdir(workdir)

  scriptsfolder = fullpath+'/tbsw/x0imaging'

  # copy cfg text file to current work directory
  default_cfg_file_name='x0.cfg'
  shutil.copy(fullpath+'/'+cfgfile,default_cfg_file_name)

  # Read config txt file
  config = ConfigParser(delimiters=(':'))
  fp = open(default_cfg_file_name)
  config.readfp(fp)

  # Read out all the relevant parameters
  # calibration factor
  calibrationfactor = config.getfloat('general', 'lambda')

  # Mean value of the momentum of the beam particles
  momentumoffset = config.getfloat('general', 'momentumoffset')

  # u momentum slope of the beam particles
  momentumugradient = config.getfloat('general', 'momentumugradient')

  # u momentum slope of the beam particles
  momentumvgradient = config.getfloat('general', 'momentumvgradient')

  # Name of the results file
  this_resultsfilename = config.get('x0image', 'resultsfilename')

  # Max fit chi2ndof value
  this_maxchi2ndof = config.get('x0image', 'maxchi2ndof')

  # fit range parameter value
  this_fitrangeparameter = config.get('general', 'fitrange_parameter')

  # Vertex multiplicity cut parameters
  this_vertex_multiplicity_min = config.get('general', 'vertex_multiplicity_min')
  this_vertex_multiplicity_max = config.get('general', 'vertex_multiplicity_max')

  # u and v length of complete X0 image 
  u_length = config.getfloat('x0image', 'u_length')
  v_length = config.getfloat('x0image', 'v_length')

  # umin and vmax of the complete X0 image
  umin = config.getfloat('x0image', 'umin')
  vmax = config.getfloat('x0image', 'vmax')

  # Pixel sizes of the image
  u_pixel_size = config.getfloat('x0image', 'u_pixel_size')
  v_pixel_size = config.getfloat('x0image', 'v_pixel_size')

  # histosettings
  num_bins = config.getfloat('x0image', 'num_bins')
  histo_range = config.getfloat('x0image', 'histo_range')

  # Fit options (log likelihood or chi2 fit?)
  fit_options = config.get('x0image', 'fit_options')

  fp.close()

  # Create X0 root file link in the current work dir
  os.symlink(fullpath+'/'+rootfile,this_x0filename)

  # number of pixels needed
  num_u_pixels = u_length * 1000.0 / u_pixel_size
  num_v_pixels = v_length * 1000.0 / v_pixel_size

  # Max number of pixels per x0image
  max_u_pixels = 75
  max_v_pixels = 75

  # Convert the number to a string
  this_maxupixels=str(max_u_pixels)
  this_maxvpixels=str(max_v_pixels)

  # How many x0image parts are required in each direction?
  u_splits=math.ceil(num_u_pixels/max_u_pixels)
  v_splits=math.ceil(num_v_pixels/max_v_pixels)

  # Overall number of X0 image parts required
  num_parts=u_splits * v_splits

  # Copy the results cfg file from previous calibrations, if it exists
  calicfgfilename="x0cal_result.cfg"
  calicfgfile=fullpath+'/localDB/'+caltag+'/'+calicfgfilename
  if os.path.isfile(calicfgfile):
    shutil.copy(calicfgfile, calicfgfilename)

  scriptname="x0imaging.C"

  # copy x0image files 
  for i in range(0,int(u_splits)):
    for j in range(0,int(v_splits)):
      # Function name for this x0image script
      this_functionname="x0imaging-part-"+str(i+1)+"-"+str(j+1)
      # fileidentifier for this x0image
      this_x0fileidentifier="-part-"+str(i+1)+"-"+str(j+1)
      # umin for this x0image
      this_umin=str(umin + i * max_u_pixels * u_pixel_size / 1000.0)
      # vmax for this x0image
      this_vmax=str(vmax-j * max_v_pixels * v_pixel_size / 1000.0)
      # ulength for this x0image
      this_ulength=str(max_u_pixels * u_pixel_size / 1000.0)
      # vlength for this x0image
      this_vlength=str(max_v_pixels * v_pixel_size / 1000.0)
      # calibrationfactor
      this_calibrationfactor=str(calibrationfactor)
      # mean momentum
      this_momentumoffset=str(momentumoffset)
      # momentum u gradient
      this_momentumugradient=str(momentumugradient)
      # momentum v gradient
      this_momentumvgradient=str(momentumvgradient)

      # Angle histo range
      this_histo_range=str(histo_range)

      # fit options
      this_fit_options=str(fit_options)

      # Number of angle histo bins
      this_num_bins=str(num_bins)
      
      # Copy placeholder x0image and change it
      shutil.copy(scriptsfolder+'/'+scriptname, scriptname)

      # Open new cfg txt file
      open('x0image-partial.cfg', 'a').close()
      fp = open('x0image-partial.cfg','r+')
      config.readfp(fp)

      config.remove_section('image')
      config.remove_section('x0image')

      # fill a temporary cfg file with the parameters relevant for the current partial x0image
      config.add_section('image')
      config.set('image', 'x0filename', this_x0filename)
      config.set('image', 'x0fileidentifier', this_x0fileidentifier)
      config.set('image', 'maxchi2ndof', this_maxchi2ndof)
      config.set('image', 'fitrange_parameter', this_fitrangeparameter)
      config.set('image', 'umin', this_umin)
      config.set('image', 'vmax', this_vmax)
      config.set('image', 'maxupixels', this_maxupixels)
      config.set('image', 'maxvpixels', this_maxvpixels)
      config.set('image', 'ulength', this_ulength)
      config.set('image', 'vlength', this_vlength)
      config.set('image', 'lambda', this_calibrationfactor)
      config.set('image', 'momentumoffset', this_momentumoffset)
      config.set('image', 'momentumugradient', this_momentumugradient)
      config.set('image', 'momentumvgradient', this_momentumvgradient)
      config.set('image', 'vertexmultiplicitymin', this_vertex_multiplicity_min)
      config.set('image', 'vertexmultiplicitymax', this_vertex_multiplicity_max)
      config.set('image', 'num_bins', this_num_bins)
      config.set('image', 'histo_range', this_histo_range)
      config.set('image', 'fit_options', this_fit_options)

      # Writing the configuration file
      with open('x0image-partial.cfg', 'w') as configfile:
         config.write(configfile,space_around_delimiters=False)
  
      subprocess.call('root -q -b '+scriptname, shell=True)
      print ('[Print] One part done ...')

  # Modify the merging script
  scriptname="MergeImages.C"

  # Copy placeholder x0image and change it
  shutil.copy(scriptsfolder+'/MergeImages.C', scriptname)

  # Open the x0imaging script of this x0image part
  tempFile = None

  with open(scriptname, 'r') as file:
    tempFile = file.read()

  # Create strings, which will replace the placeholder strings
  this_usplits=str(u_splits)
  this_vsplits=str(v_splits)
  this_upixelsize=str(u_pixel_size)
  this_vpixelsize=str(v_pixel_size)
  this_umin=str(umin)
  this_vmax=str(vmax)
  this_ulength=str(u_length)
  this_vlength=str(v_length)

  # Open new cfg txt file
  open('x0merge.cfg', 'a').close()
  fp = open('x0merge.cfg','r+')
  config.readfp(fp)

  config.remove_section('image')
  config.remove_section('x0image')

  # fill a temporary cfg file with the parameters relevant for the current partial x0image
  config.add_section('merge')
  config.set('merge', 'x0filename', this_x0filename)
  config.set('merge', 'resultsfilename', this_resultsfilename)
  config.set('merge', 'usplits', this_usplits)
  config.set('merge', 'vsplits', this_vsplits)
  config.set('merge', 'upixelsize', this_upixelsize)
  config.set('merge', 'vpixelsize', this_vpixelsize)
  config.set('merge', 'maxupixels', this_maxupixels)
  config.set('merge', 'maxvpixels', this_maxvpixels)
  config.set('merge', 'umin', this_umin)
  config.set('merge', 'vmax', this_vmax)

  # Writing the configuration file
  with open('x0merge.cfg', 'w') as configfile:
    config.write(configfile,space_around_delimiters=False)
   
  subprocess.call('root -q -b '+scriptname, shell=True)            
  print ('[Print] All partial images created and merged... ')  

  # clean up partial image scripts 
  for tmpfile in glob.glob('*part_*.C'):
    os.remove(tmpfile) 

  if deletetag == '1':
    # clean up root files
    for tmpfile in glob.glob('*part*.root'):
      os.remove(tmpfile) 
    # also remove the input cfg file
    os.remove(default_cfg_file_name) 

  # remove MergeImage.C and config file script
  os.remove('x0merge.cfg') 
  os.remove('x0image-partial.cfg')  
  for file in glob.glob('*.C'):
    os.remove(file) 

  # Copy complete image file to root-files directory
  shutil.copy(this_resultsfilename+'.root', fullpath+'/root-files/'+this_x0filename+'-'+ tag +'X0image.root')
  		           
                


 	


