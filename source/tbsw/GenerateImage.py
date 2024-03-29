import os
import shutil
import subprocess
import sys
import glob
import math
import fileinput
import argparse
import configparser

class MyConfigParser(configparser.ConfigParser):

    def write(self, fp):
        """Write an .ini-format representation of the configuration state."""
        if self._defaults:
            fp.write("[%s]\n" % DEFAULTSECT)
            for (key, value) in self._defaults.items():
                fp.write("%s = %s\n" % (key, str(value).replace('\n', '\n\t')))
            fp.write("\n")
        for section in self._sections:
            fp.write("[%s]\n" % section)
            for (key, value) in self._sections[section].items():
                if key == "__name__":
                    continue
                if (value is not None) or (self._optcre == self.OPTCRE):
                    key = ": ".join((key, str(value).replace('\n', '\n\t')))
                fp.write("%s\n" % (key))
            fp.write("\n")

def x0imaging(rootfilelist=[],caltag='',steerfiles='',nametag=''):
  """
  Performs radiation length imaging from reconstructed root files in filelist
    :@rootlist:    List of input root files with TTree of scattering angles
    :@caltag:      Use x0 calibration parameters from cfg files in the caltag subdirectory    
    :@steerfiles:  Directory, where the image cfg file lies
    :@nametag:     Name tag, that is used to modify the output filename 
    :author: ulf.stolzenberg@phys.uni-goettingen.de  
  """   
  
  cfgfile=steerfiles+'x0.cfg'
   
  # Remember current working dir 
  fullpath = os.getcwd() 
  
  # Change to work directory
  workdir = 'tmp-runs/'+ nametag
  
  # Remove old workdir if exists 
  if os.path.isdir(workdir):
    shutil.rmtree(workdir)
    
  # Create workdir and change directory
  os.mkdir(workdir) 
  os.chdir(workdir)

  # Copy cfg text file to current work directory
  default_cfg_file_name='x0.cfg'
  shutil.copy(fullpath+'/'+cfgfile,default_cfg_file_name)

  # Read config txt file
  #config = ConfigParser.ConfigParser(delimiters=(':'))
  config = MyConfigParser()
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

  # Weight of energy loss calculation
  epsilon = config.get('general', 'epsilon')

  # Name of the results file
  resultsfilename = config.get('x0image', 'resultsfilename')

  # Max fit chi2ndof value
  maxchi2ndof = config.get('x0image', 'maxchi2ndof')

  # fit range parameter value
  fitrangeparameter = config.get('general', 'fitrange_parameter')

  # Vertex multiplicity cut parameters
  vertex_multiplicity_min = config.get('general', 'vertex_multiplicity_min')
  vertex_multiplicity_max = config.get('general', 'vertex_multiplicity_max')

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
  min_tracks = config.getfloat('x0image', 'min_tracks')

  # Fit options (log likelihood or chi2 fit?)
  fit_options = config.get('x0image', 'fit_options')

  fp.close()

  # Create Links of input root files
  for rootfile in rootfilelist:
    this_x0filename = os.path.splitext(os.path.basename(rootfile))[0]
    os.symlink(fullpath+'/'+rootfile,this_x0filename)  

  # number of pixels needed
  num_u_pixels = u_length * 1000.0 / u_pixel_size
  num_v_pixels = v_length * 1000.0 / v_pixel_size

  # Max number of pixels per x0image
  max_u_pixels = 50
  max_v_pixels = 50

  # How many x0image parts are required in each direction?
  u_splits=math.ceil(num_u_pixels/max_u_pixels)
  v_splits=math.ceil(num_v_pixels/max_v_pixels)

  # Copy the results cfg file from previous calibrations, if it exists
  calicfgfilename="x0cal_result.cfg"
  calicfgfile=fullpath+'/localDB/'+caltag+'/'+calicfgfilename
  if os.path.isfile(calicfgfile):
    shutil.copy(calicfgfile, calicfgfilename)
  
  # Copy x0image files 
  for i in range(0,int(u_splits)):
    for j in range(0,int(v_splits)):

      # Input files
      for it,rootfile in enumerate(rootfilelist):
        if it==0:
          this_inputfiles=os.path.splitext(os.path.basename(rootfile))[0]
        else:
          this_inputfiles=this_inputfiles+","+os.path.splitext(os.path.basename(rootfile))[0]
      
      # fileidentifier for this x0image
      this_x0imagename=nametag+"-part-"+str(i+1)+"-"+str(j+1)+".root"
      # umin for this x0image
      this_umin=str(umin + i * max_u_pixels * u_pixel_size / 1000.0)
      # vmax for this x0image
      this_vmax=str(vmax-j * max_v_pixels * v_pixel_size / 1000.0)
      # ulength for this x0image
      this_ulength=str(max_u_pixels * u_pixel_size / 1000.0)
      # vlength for this x0image
      this_vlength=str(max_v_pixels * v_pixel_size / 1000.0)
        
      # Open new cfg txt file
      open('x0image-partial.cfg', 'a').close()
      fp = open('x0image-partial.cfg','r+')
      config.readfp(fp)

      config.remove_section('image')
      config.remove_section('x0image')
      config.remove_section('x0calibration')

      # fill a temporary cfg file with the parameters relevant for the current partial x0image
      config.add_section('image')
      config.set('image', 'inputfile', this_inputfiles)
      config.set('image', 'x0imagename', this_x0imagename)
      config.set('image', 'maxchi2ndof', maxchi2ndof)
      config.set('image', 'fitrange_parameter', fitrangeparameter)
      config.set('image', 'umin', this_umin)
      config.set('image', 'vmax', this_vmax)
      config.set('image', 'maxupixels', str(max_u_pixels))
      config.set('image', 'maxvpixels', str(max_v_pixels))
      config.set('image', 'ulength', this_ulength)
      config.set('image', 'vlength', this_vlength)
      config.set('image', 'lambda', str(calibrationfactor))
      config.set('image', 'momentumoffset', str(momentumoffset))
      config.set('image', 'momentumugradient', str(momentumugradient))
      config.set('image', 'momentumvgradient', str(momentumvgradient))
      config.set('image', 'epsilon', str(epsilon))
      config.set('image', 'vertexmultiplicitymin', vertex_multiplicity_min)
      config.set('image', 'vertexmultiplicitymax', vertex_multiplicity_max)
      config.set('image', 'min_tracks', str(min_tracks))
      config.set('image', 'num_bins', str(num_bins))
      config.set('image', 'histo_range', str(histo_range))
      config.set('image', 'fit_options', str(fit_options))

      # Writing the configuration file
      with open('x0image-partial.cfg', 'w') as configfile:
         #config.write(configfile,space_around_delimiters=False)
         config.write(configfile)
  
      print ('[INFO] Create x0image for image part {} {}'.format(i,j))    
      subprocess.call( "$MARLIN/bin/x0imaging > x0-imaging_{}_{}.log 2>&1".format(i,j), shell=True)
      
  # Open new cfg txt file
  open('x0merge.cfg', 'a').close()
  fp = open('x0merge.cfg','r+')
  config.readfp(fp)

  config.remove_section('image')
  config.remove_section('x0image')

  # fill a temporary cfg file with the parameters relevant for the current partial x0image
  config.add_section('merge')
  config.set('merge', 'x0filename', nametag)
  config.set('merge', 'resultsfilename', resultsfilename)
  config.set('merge', 'usplits', str(u_splits))
  config.set('merge', 'vsplits', str(v_splits))
  config.set('merge', 'upixelsize', str(u_pixel_size))
  config.set('merge', 'vpixelsize', str(v_pixel_size))
  config.set('merge', 'maxupixels', str(max_u_pixels))
  config.set('merge', 'maxvpixels', str(max_v_pixels))
  config.set('merge', 'umin', str(umin))
  config.set('merge', 'vmax', str(vmax))

  # Writing the configuration file
  with open('x0merge.cfg', 'w') as configfile:
    #config.write(configfile,space_around_delimiters=False)
    config.write(configfile)
   
  subprocess.call( "$MARLIN/bin/MergeImages > merger.log 2>&1", shell=True)            
  print ('[Print] All partial images created and merged... ')  

  # clean up partial image scripts 
  for tmpfile in glob.glob('*part_*.C'):
    os.remove(tmpfile) 

  # remove MergeImage.C and config file script
  os.remove('x0merge.cfg') 
  os.remove('x0image-partial.cfg')  
  for file in glob.glob('*.C'):
    os.remove(file) 

  # Copy complete image file to root-files directory
  shutil.copy(resultsfilename+'.root', fullpath+'/root-files/'+nametag+'.root')
  		                      
  # Go back to workspace
  os.chdir(fullpath) 

 	


