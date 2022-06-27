#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:03:13 2022

@author: astrostudent
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 16:47:47 2022

@author: bezawitmekashakassaye
"""

#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Sat Thursday January 13 11:12:40 2022

@author: Bezawit Mekasha Kassaye

What the program does: 
- Reads in CSV file with the YBs informations 
- Plots Guassian fit, save it and show it to the user
- Takes an input of whether it was a good fit or not 
- Creates a CSV file output that contains the users input of yes/no
- It also let's the user start where they left off or begin again as well as revist any YB.

"""



import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting

from astropy.io import fits
from astropy.wcs import WCS
from astropy.io import ascii 
from operator import add

import os.path



import glob
import imageio
import matplotlib.image as mpimg

#checks if new.catalog.csv exists if does it asks if the user 
# wants to start again, visit a specific YB or start where they left off. 

# if the path doesn't exist it asks the user if they want to start 
#from the begining or visit a specific YB.
if os.path.exists('new.catalog.csv'): 
    begin = input("Enter 'restart' if you want to star over. Enter 'a' if you want to vist a specific YB. Anything else to start where you left off: ")
    if begin == 'Restart' or begin == 'restart' or begin == 'RESTART':  
        h_catalog = pd.read_csv('vel_with_gf.csv') #should be vel_with_gf.csv
        h_catalog["Good fit"] = ""
        h_catalog ["Multiple emission line"] = ""
        YB_inpt = h_catalog['YB'][0]
    elif begin =='a' or begin =='A':
        h_catalog = pd.read_csv('new.catalog.csv')
        YB_inpt = int(input("Enter the YB number you want to revist: "))         
    else: 
        h_catalog = pd.read_csv('new.catalog.csv')
        YB_inpt = h_catalog['YB'][0]
else: 
    begin = input("Enter 'a' if you want to vist a YB, anything else to start from the beginning: ")
    if begin == 'a' or begin == 'A':
        h_catalog = pd.read_csv('vel_with_gf.csv') # should be vel_with_gf.csv#quoting=csv.QUOTE_NONE, error_bad_lines=False, header = None, names = ['YB','GLON',  'GLAT', 'r', 'Vstrongest (km/s)','g amplitude', 'g velocity','g stddev','uncertamp','uncertgvel','uncertstd', 'Good fit']) #modified from 'rb' to 'r' because of this error Error: iterator should return strings, not bytes (did you open the file in text mode?)
        h_catalog["Good fit"] = ""
        h_catalog ["Multiple emission line"] = ""
        #h_catalog = pd.read_csv('new.catalog.csv')
        YB_inpt = int(input("Enter the YB number you want to revist: "))
    else:
        h_catalog = pd.read_csv('vel_with_gf.csv') # should be vel_with_gf.csv#quoting=csv.QUOTE_NONE, error_bad_lines=False, header = None, names = ['YB','GLON',  'GLAT', 'r', 'Vstrongest (km/s)','g amplitude', 'g velocity','g stddev','uncertamp','uncertgvel','uncertstd', 'Good fit']) #modified from 'rb' to 'r' because of this error Error: iterator should return strings, not bytes (did you open the file in text mode?)
        h_catalog["Good fit"] = ""
        h_catalog ["Multiple emission line"] = ""
        YB_inpt = h_catalog['YB'][0]




class YellowBall:
   """ YellowBall objects represent a source from the YB catalog with the following attributes:
      ----------- Attributes -----------
      cube: fits object containing the data cube corresponding to the yb
      glon: Galactic Longitude (in degrees)
      glat: Galactuc Latitude (in degrees)
      r: Radius as reported in the catalog, i.e. user specified (in degrees)
      wcs: World Coordinate System object for the yb
      data: 3D array containing the flux at each pixel for every velocity "slice" of the cube.

      ----------- Member Functions -----------
      get_glon: getter for YB Galactic Longitude
      get_glat: getter for YB Galactic Latitude
      get_r: getter for YB user-specified radius
      plot_image: plot the entire fits image containing the yb for a given velocity.
      get_spectrum: return an array containing the spectrum across all 
         velocities of the single pixel at the yb coordinates.
      avg_spectrum: return an array containing the average spectrum in a box 
         composed of all the pixels corresponding to the user-specified radius.
      get_main_peak: return a tuple where the first value corresponds to the
         velocity of the main peak and the second to it's flux value.
      sd_fun:return sd_value which is noise cutoff for each cube. 
   """
   def __init__(self, yb_lon, yb_lat, yb_r, cube):

             self.glon = float(yb_lon)
             self.glat = float(yb_lat) 
            
             self.r = float(yb_r) 
             
             self.cube = cube 
             self.wcs = WCS(cube)
             self.data = cube.data

   ######################## GETTERS ########################

   def get_glon(self):
      return self.glon

   def get_glat(self):
      return self.glat

   def get_r(self):
      return self.r

   ###################### END GETTERS ######################

   def plot_image(self, vel):
      empty1, empty2, pixvel = self.wcs.all_world2pix(0.,0.,vel*1000.,0)
      pixvel = int((pixvel))
      plt.figure(1)
      ax = plt.subplot(projection=self.wcs, slices=('x', 'y', 100))
      ax.imshow((self.data[pixvel,:,:]), origin='lower')
      
      if self.cube.header['object'] == '(31.000,-1.000,-5000) to (29.000,1.000,135000)':
         cube = '29-31'
      elif self.cube.header['object'] == '(33.000,-1.000,-5000) to (31.000,1.000,135000)':
         cube = '31-33'
      else:
         cube = '33-35'

      plt.savefig('13Co_images/13co_%s.png' % cube)
      plt.show()
      plt.close()
      #plt.figure(0)

# returns spectrum at single pixel centered at YB long, lat
   def get_spectrum(self):
      pixlon, pixlat, empty = self.wcs.all_world2pix(self.glon,self.glat,0.,0)
      pixlon = int(round(pixlon))
      pixlat = int(round(pixlat))
      return self.data[:,pixlat+1,pixlon+1]

#returns average spectrum of all pixels covered by YB extent, defined by central pixel +/- radius
   def avg_spectrum(self):
      b_pixlon, b_pixlat, empty = self.wcs.all_world2pix(self.glon + self.r, self.glat - self.r, 0., 0) 
      e_pixlon, e_pixlat, empty = self.wcs.all_world2pix(self.glon - self.r, self.glat + self.r, 0., 0)
      #converting the glon and glat value to pixle values for plotting. 
      #b_pixlon and b_pixlat are the beginnings while e_pixlon and e_pixlat are the endings for the plot. 
      
      b_pixlon = int((b_pixlon)) 
      e_pixlon = int((e_pixlon))
      b_pixlat = int((b_pixlat))
      e_pixlat = int((e_pixlat))
      
      box= []
      
      for i in range(b_pixlon, e_pixlon+1): 
         for j in range(b_pixlat, e_pixlat+1):
            box.append([i,j])
      sum_spectrum = [0.]*(self.data.shape[0]-1) # The number of elements in the velocity dimention - 1.

      for pixel in box:
         spectrum = [x for x in self.data[:, pixel[1],pixel[0]]] 
         sum_spectrum = map(add, spectrum, sum_spectrum)
      avg_spectrum = [x/len(box) for x in sum_spectrum]
      return avg_spectrum
  
   def get_main_peak(self):
      avg_spectrum = self.avg_spectrum() #where the previous fun is being called. 
      max_flux = max(avg_spectrum)
      pixvel = avg_spectrum.index(max_flux)
      empty1, empty2, vel = self.wcs.all_pix2world(0.,0.,pixvel,0)
      vel = vel / 1000
      return (vel, max_flux)
 
   
   def sd_fun (self,YB_number, bvel, evel,cube):
       #bvel is the begining velocity
       #evel is the ending velocity
       self.wcs = WCS(cube)
       self.data = cube.data

       #define step as the km/s per velocity channel; GRS cubes are varying vel chan length 

       step = cube.header['CDELT3']/1000.

       bchan = int(bvel/step) #coverting the bvel value to velocity channel 

       echan = int(evel/step)
       temparray = self.avg_spectrum()
       subarray=temparray[bchan:echan]
       sd_value = np.std(subarray)
#
       return sd_value


  

    
if __name__ == '__main__':
    #SEDIGISM Data server
    cube0 = fits.open("data_cubes/G000_13CO21_Tmb_DR1.fits")[0]
    
    cube02= fits.open("data_cubes/G002_13CO21_Tmb_DR1.fits")[0]
    
    cube04 = fits.open("data_cubes/G004_13CO21_Tmb_DR1.fits")[0]
    
    cube06 = fits.open("data_cubes/G006_13CO21_Tmb_DR1.fits")[0]
   
    cube08 = fits.open("data_cubes/G008_13CO21_Tmb_DR1.fits")[0]
   
    cube10 = fits.open("data_cubes/G010_13CO21_Tmb_DR1.fits")[0]

    cube12 = fits.open("data_cubes/G012_13CO21_Tmb_DR1.fits")[0]
    cube14 = fits.open("data_cubes/G014_13CO21_Tmb_DR1.fits")[0]
   
   
    #Boston University GRS DATA 
    cube16 = fits.open("data_cubes/13COdatacube_015-017.fits")[0]
    cube18 = fits.open("data_cubes/13COdatacube_017-019.fits")[0]
    cube20 = fits.open("data_cubes/13COdatacube_019-021.fits")[0]
    cube22 = fits.open("data_cubes/13COdatacube_021-023.fits")[0]
    cube24 = fits.open("data_cubes/13COdatacube_023-025.fits")[0]
    cube26 = fits.open("data_cubes/13COdatacube_025-027.fits")[0]
    cube28 = fits.open("data_cubes/13COdatacube_027-029.fits")[0]
    cube30 = fits.open("data_cubes/13COdatacube_029-031.fits")[0]
    cube32 = fits.open("data_cubes/13COdatacube_031-033.fits")[0]
    cube34 = fits.open("data_cubes/13COdatacube_033-035.fits")[0]
    cube36 = fits.open("data_cubes/13COdatacube_035-037.fits")[0]
    cube38 = fits.open("data_cubes/13COdatacube_037-039.fits")[0]
    cube40 = fits.open("data_cubes/13COdatacube_039-041.fits")[0]
    cube42 = fits.open("data_cubes/13COdatacube_041-043.fits")[0]
    cube44 = fits.open("data_cubes/13COdatacube_043-045.fits")[0]
    cube46 = fits.open("data_cubes/13COdatacube_045-047.fits")[0]
    cube48 = fits.open("data_cubes/13COdatacube_047-049.fits")[0]
    cube50 = fits.open("data_cubes/13COdatacube_049-051.fits")[0]
    cube52 = fits.open("data_cubes/13COdatacube_051-053.fits")[0]
    cube54 = fits.open("data_cubes/13COdatacube_053-055.fits")[0]
    cube56 = fits.open("data_cubes/13COdatacube_055-057.fits")[0]
   
   
   
   #SEDIGISM Data server

    cube302 = fits.open("data_cubes/G302_13CO21_Tmb_DR1.fits")[0]
    cube304 = fits.open("data_cubes/G304_13CO21_Tmb_DR1.fits")[0]
    cube306 = fits.open("data_cubes/G306_13CO21_Tmb_DR1.fits")[0]
    cube308 = fits.open("data_cubes/G308_13CO21_Tmb_DR1.fits")[0]
    cube310 = fits.open("data_cubes/G310_13CO21_Tmb_DR1.fits")[0]
    cube312 = fits.open("data_cubes/G312_13CO21_Tmb_DR1.fits")[0]
    cube314 = fits.open("data_cubes/G314_13CO21_Tmb_DR1.fits")[0]
    cube316 = fits.open("data_cubes/G316_13CO21_Tmb_DR1.fits")[0]
    cube318 = fits.open("data_cubes/G318_13CO21_Tmb_DR1.fits")[0]
    cube320 = fits.open("data_cubes/G320_13CO21_Tmb_DR1.fits")[0]
    cube322 = fits.open("data_cubes/G322_13CO21_Tmb_DR1.fits")[0]
    cube324 = fits.open("data_cubes/G324_13CO21_Tmb_DR1.fits")[0]
    cube326 = fits.open("data_cubes/G326_13CO21_Tmb_DR1.fits")[0]
    cube328 = fits.open("data_cubes/G328_13CO21_Tmb_DR1.fits")[0]
    cube330 = fits.open("data_cubes/G330_13CO21_Tmb_DR1.fits")[0]
    cube332 = fits.open("data_cubes/G332_13CO21_Tmb_DR1.fits")[0]
    cube334 = fits.open("data_cubes/G334_13CO21_Tmb_DR1.fits")[0]
    cube336 = fits.open("data_cubes/G336_13CO21_Tmb_DR1.fits")[0]
    cube338 = fits.open("data_cubes/G338_13CO21_Tmb_DR1.fits")[0]
    cube340 = fits.open("data_cubes/G340_13CO21_Tmb_DR1.fits")[0]
    cube342 = fits.open("data_cubes/G342_13CO21_Tmb_DR1.fits")[0]
    cube344 = fits.open("data_cubes/G344_13CO21_Tmb_DR1.fits")[0]
    cube346 = fits.open("data_cubes/G346_13CO21_Tmb_DR1.fits")[0]
    cube348 = fits.open("data_cubes/G348_13CO21_Tmb_DR1.fits")[0]
    cube350 = fits.open("data_cubes/G350_13CO21_Tmb_DR1.fits")[0]
    cube352 = fits.open("data_cubes/G352_13CO21_Tmb_DR1.fits")[0]
    cube354 = fits.open("data_cubes/G354_13CO21_Tmb_DR1.fits")[0]
    cube356 = fits.open("data_cubes/G356_13CO21_Tmb_DR1.fits")[0]
    cube358 = fits.open("data_cubes/G358_13CO21_Tmb_DR1.fits")[0]
    cube359 = fits.open("data_cubes/G359_13CO21_Tmb_DR1.fits")[0]
   





i=0

while i< len(h_catalog):       
    YB = int(h_catalog['YB'][i])   
    YB_long = float(h_catalog['GLON'][i])
    YB_lat = float(h_catalog['GLAT'][i])
    YB_rad = float(h_catalog['r'][i])
    vlsr = h_catalog['Vstrongest (km/s)'][i]
    if h_catalog['g velocity'][i] != 'Data_NA':
        #vlsr = float(h_catalog['Vstrongest (km/s)'][i])
        gvel = float(h_catalog['g velocity'][i])
        gamp = float(h_catalog['g amplitude'][i])
        gstddev = float(h_catalog['g stddev'][i])
    GF = h_catalog['Good fit'][i]

  
       # #SEDGISM dataset 
      
    if YB_long >= -1 and YB_long < 1 and YB_lat >= -0.5 and YB_lat <= 0.5: #based on the cube's longitude and lattitude range. 
          yb=YellowBall(YB_long, YB_lat, YB_rad, cube0) #Calling the class depending on which cube the YBs fit in
          nvel = cube0.header['naxis3'] -1 #nvel is the number of velocity channels in the cube, and the information is extracted from the header of the cube fits files. 
    elif YB_long >=1 and YB_long < 3 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube02)
          nvel = cube02.header['naxis3'] -1
    elif YB_long >=3 and YB_long < 5 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube04)
          nvel = cube04.header['naxis3'] -1
    elif YB_long >=5 and YB_long < 7  and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube06)
          nvel = cube06.header['naxis3'] -1
    elif YB_long >=7 and YB_long < 9 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube08)
          nvel = cube08.header['naxis3'] -1
    elif YB_long >=9 and YB_long < 11 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube10)
          nvel = cube10.header['naxis3'] -1
    elif YB_long >=11 and YB_long < 13 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube12)
          nvel = cube12.header['naxis3'] -1
    elif YB_long >=13 and YB_long< 15 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube14)
          nvel = cube14.header['naxis3'] -1
    
       
       #BU dataset     
    elif YB_long >=15 and YB_long < 17 and YB_lat >= -0.806 and YB_lat <0.818: 
         yb=YellowBall(YB_long,YB_lat,YB_rad,cube16)
         nvel = cube16.header['naxis3'] -1
    elif YB_long >=17 and YB_long < 19:
         yb=YellowBall(YB_long,YB_lat,YB_rad,cube18)
         nvel = cube18.header['naxis3'] -1
    elif YB_long >=19 and YB_long< 21:
         yb=YellowBall(YB_long,YB_lat,YB_rad,cube20)
         nvel = cube20.header['naxis3'] -1
    elif YB_long >=21 and YB_long < 23:
         yb=YellowBall(YB_long,YB_lat,YB_rad,cube22)
         nvel = cube22.header['naxis3'] -1
    elif YB_long >=23 and YB_long < 25:
         yb=YellowBall(YB_long,YB_lat,YB_rad,cube24)
         nvel = cube24.header['naxis3'] -1
    elif YB_long >=25 and YB_long < 27:
         yb=YellowBall(YB_long,YB_lat,YB_rad,cube26)
         nvel = cube26.header['naxis3'] -1
    elif YB_long >=27 and YB_long < 29:
         yb=YellowBall(YB_long,YB_lat,YB_rad,cube28)
         nvel = cube28.header['naxis3'] -1 
    elif YB_long >=29 and YB_long< 31:
        yb=YellowBall(YB_long,YB_lat,YB_rad,cube30)
        nvel = cube30.header['naxis3'] -1  
    elif YB_long >=31 and YB_long< 33:
        yb=YellowBall(YB_long,YB_lat,YB_rad,cube32)
        nvel = cube32.header['naxis3'] -1
    elif YB_long>=33 and YB_long <35:
        yb=YellowBall(YB_long,YB_lat,YB_rad,cube34)
        nvel = cube34.header['naxis3'] -1
    elif YB_long >=35 and YB_long <37:
        yb=YellowBall(YB_long,YB_lat,YB_rad, cube36)
        nvel = cube36.header['naxis3'] -1
    elif YB_long >=37 and YB_long<39:
        yb=YellowBall(YB_long,YB_lat,YB_rad, cube38)
        nvel = cube38.header['naxis3'] -1
    elif YB_long >=39 and YB_long <41:
        yb=YellowBall(YB_long,YB_lat,YB_rad,cube40)
        nvel = cube40.header['naxis3'] -1
    elif YB_long >=41 and YB_long <43:
        yb=YellowBall(YB_long,YB_lat,YB_rad,cube42)
        nvel = cube42.header['naxis3'] -1
    elif YB_long >=43 and YB_long <45:
        yb=YellowBall(YB_long,YB_lat,YB_rad,cube44)
        nvel = cube44.header['naxis3'] -1
    elif YB_long >=45 and YB_long<47:
        yb=YellowBall(YB_long,YB_lat,YB_rad,cube46)
        nvel = cube46.header['naxis3'] -1
    elif YB_long>=47 and YB_long<49:
        yb=YellowBall(YB_long,YB_lat,YB_rad,cube48)
        nvel = cube48.header['naxis3'] -1
       
    elif YB_long >=49 and YB_long< 51:
        yb=YellowBall(YB_long,YB_lat,YB_rad,cube50)
        nvel = cube50.header['naxis3'] -1
    elif YB_long >=51 and YB_long< 53:
        yb=YellowBall(YB_long,YB_lat,YB_rad,cube52)
        nvel = cube52.header['naxis3'] -1
    elif YB_long >=53 and YB_long < 55:
        yb=YellowBall(YB_long,YB_lat,YB_rad,cube54)
        nvel = cube54.header['naxis3'] -1
    elif YB_long >=55 and YB_long< 55.707 and YB_lat >= -1.095 and YB_lat<1.095:
       yb=YellowBall(YB_long,YB_lat,YB_rad,cube56)
       nvel = cube56.header['naxis3'] -1
       print(nvel)
       
      # #SEDISIM data set 
      
    elif YB_long >=301 and YB_long <303 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube302)
          nvel = cube302.header['naxis3'] -1
    elif YB_long >=303 and YB_long <305  and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube304)
          nvel = cube304.header['naxis3'] -1
    elif YB_long >=305 and YB_long < 307 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube306)
          nvel = cube306.header['naxis3'] -1
    elif YB_long >=307 and YB_long< 309 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube308)
          nvel = cube308.header['naxis3'] -1
    elif YB_long >=309 and YB_long < 311 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube310)
          nvel = cube310.header['naxis3'] -1
    elif YB_long >=311 and YB_long < 313 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube312)
          nvel = cube312.header['naxis3'] -1
    elif YB_long >=313 and YB_long < 315 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube314)
          nvel = cube314.header['naxis3'] -1
    elif YB_long >=315 and YB_long < 317 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube316)
          nvel = cube316.header['naxis3'] -1
    elif YB_long >=317 and YB_long < 319 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube318)
          nvel = cube318.header['naxis3'] -1
    elif YB_long >=319 and YB_long < 321 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube320)
          nvel = cube320.header['naxis3'] -1
    elif YB_long >=321 and YB_long < 323 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube322)
          nvel = cube322.header['naxis3'] -1
    elif YB_long >=323 and YB_long < 325 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube324)
          nvel = cube324.header['naxis3'] -1
    elif YB_long >=325 and YB_long < 327 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube326)
          nvel = cube326.header['naxis3'] -1
    elif YB_long >=327 and YB_long < 329 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube328)
          nvel = cube328.header['naxis3'] -1
    elif YB_long >=329 and YB_long < 331 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube330)
          nvel = cube330.header['naxis3'] -1    
    elif YB_long >=331 and YB_long < 333 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube332)
          nvel = cube332.header['naxis3'] -1
    elif YB_long >=333 and YB_long < 335 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube334)
          nvel = cube334.header['naxis3'] -1
    elif YB_long >=335 and YB_long < 337 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube336)
          nvel = cube336.header['naxis3'] -1
    elif YB_long >=337 and YB_long < 339 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube338)
          nvel = cube338.header['naxis3'] -1
    elif YB_long >=339 and YB_long < 341 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube340)
          nvel = cube340.header['naxis3'] -1
    elif YB_long >=341 and YB_long < 343 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube342)
          nvel = cube342.header['naxis3'] -1
    elif YB_long >=343 and YB_long < 345 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube344)
          nvel = cube344.header['naxis3'] -1
    elif YB_long >= 345 and YB_long < 347 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube346)
          nvel = cube346.header['naxis3'] -1
    elif YB_long >=347 and YB_long < 349 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube348)
          nvel = cube348.header['naxis3'] -1
    elif YB_long >=349 and YB_long < 351 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube350)
          nvel = cube350.header['naxis3'] -1
    elif YB_long >=351 and YB_long < 353 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube352)
          nvel = cube352.header['naxis3'] -1
    elif YB_long >=353 and YB_long < 355 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube354)
          nvel = cube354.header['naxis3'] -1
    elif YB_long >=355 and YB_long < 357 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube356)
          nvel = cube356.header['naxis3'] -1
    elif YB_long >=357 and YB_long < 359 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube358)
          nvel = cube358.header['naxis3'] -1    
    elif YB_long >=359 and YB_long < 360 and YB_lat >= -0.5 and YB_lat <= 0.5:
          yb=YellowBall(YB_long,YB_lat,YB_rad,cube359)
          nvel = cube359.header['naxis3'] -1


    box_spec = yb.avg_spectrum()
    step= yb.cube.header['CDELT3']/1000
    main_peak_vel, max_flux = yb.get_main_peak()
     

  
    if begin == 'a' or begin == 'A': 
        
        if YB== YB_inpt:
            if YB_long >=15 and YB_long <41:
                #while (True):
                    print('YellowBall ID: ', YB)
               
                   
                    z = [(x)*step - 5.0 for x in range(nvel)]
                    y = (box_spec) #declaring z and y and bounding their std,mean and amplitude.
          #Bounding is helpful expecially if we have a wide peak or multiple peaks
                    def gaussian(x, a, mu, sig):
                        return a* (np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))

                    x_values = np.linspace(-5, 135, nvel)
                    #for mu, sig in [(-1, 1), (0, 2), (2, 3)]:
                    plt.plot(z,y, color='k')

                    plt.plot(z, gaussian(x_values, gamp, gvel, gstddev), color ='darkorange')
                        
                    
                    plt.savefig(('TRIAL_amplitude is bounded and sd is 2/box_spec_%03i' % (YB)))
                    #plt.clf() 
                    
                    plt.show()                      
                    #coords = gf_fig.ginput(n=0, timeout=300, show_clicks=False,mouse_stop=2)
    
                    #check = input("Is this a good fit? Y/N: " )
                    check = input("hit Enter to skip if a good fit. Anything else for a bad fit: " )

                   #if check == 'y' or check =='Y' :
                    if check == "":
                        result = 'yes'
                    else:
                        result = 'no'
                    mel = input("Enter 'y' if there are multiple strong emission lines, if not hit 'enter' to coninue ")
                    if mel == 'y' or mel == 'Y':
                        mel_result = 'Yes'
                    else:
                        mel_result = 'No'
                  
                    plt.close()
            elif YB_long >=41 and YB_long <57.707:
                #while (True):
                    print('YellowBall ID: ', YB)
               
                   
                    z = [(x)*step - 5.0 for x in range(nvel)]
                    y = (box_spec) #declaring z and y and bounding their std,mean and amplitude.
          #Bounding is helpful expecially if we have a wide peak or multiple peaks
                    def gaussian(x, a, mu, sig):
                        return a* (np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))

                    x_values = np.linspace(-5, 85, nvel)
                    #for mu, sig in [(-1, 1), (0, 2), (2, 3)]:
                    plt.plot(z,y, color='k')

                    plt.plot(z, gaussian(x_values, gamp, gvel, gstddev), color ='darkorange')
                        
                    
                    plt.savefig(('TRIAL_amplitude is bounded and sd is 2/box_spec_%03i' % (YB)))
                    #plt.clf() 
                    
                    plt.show()                      
                    #coords = gf_fig.ginput(n=0, timeout=300, show_clicks=False,mouse_stop=2)
    
                    #check = input("Is this a good fit? Y/N: " )
                    check = input("hit Enter to skip if a good fit. Anything else for a bad fit: " )

                   #if check == 'y' or check =='Y' :
                    if check == "":
                        result = 'yes'
                    else:
                        result = 'no'
                    mel = input("Enter 'y' if there are multiple strong emission lines, if not hit 'enter' to coninue ")
                    if mel == 'y' or mel == 'Y':
                        mel_result = 'Yes'
                    else:
                        mel_result = 'No'
                  
                    plt.close()
        
            else: 
                   print('YellowBall ID: ', YB)
    
    
                  
                   z = [(x)*step - 200.0 for x in range(nvel)]
                   y = (box_spec) 
                   
                   def gaussian(x, a, mu, sig):
                       return a* (np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))
                    
                   x_values = np.linspace(-200, 200, nvel)
                   #for mu, sig in [(-1, 1), (0, 2), (2, 3)]:
                   plt.plot(z,y, color='k')
        
                   plt.plot(z, gaussian(x_values, gamp, gvel, gstddev), color ='darkorange')
                                    
                   plt.savefig(('TRIAL_amplitude is bounded and sd is 2/box_spec_%03i' % (YB)))
                   #plt.clf() 
                   #plt.ion()
                   plt.show()
                   
    
                   check = input("hit Enter to skip if a good fit. Anything else for a bad fit: " )
                   #if check == 'y' or check =='Y' :
                   if check =="":
                       result = 'yes'
                   else:
                       result = 'no'
                   mel = input("Enter 'y' if there are multiple strong emission lines, if not hit 'enter' to coninue ")
                   if mel == 'y' or mel == 'Y':
                       mel_result = 'Yes'
                   else:
                       mel_result = 'No'
                  
                   plt.close()
         

            h_catalog.loc[h_catalog["YB"] == YB_inpt, 'Good fit'] = result 
            h_catalog.loc[h_catalog["YB"] == YB_inpt, 'Multiple emission line'] =mel_result
            
            
            h_catalog.to_csv("new.catalog.csv", index = False)  #pushing forward the data from h_catalog with the updates

            
            inpt = input('would you like to keep visiting a specific YB? (Y/N): ')
            if inpt == 'n' or inpt =='N':
                begin = input("Would you like to continue working through the rest of the YBs? (Y/N): ")
                if begin == 'y' or begin =='Y':
                    i=0
                else:
                    break
            else: 
                YB_inpt = int(input("Enter the YB you want to visit: "))
                i = 0
                
        
    else: 
        if vlsr == 'Data_NA':
            h_catalog.loc[h_catalog["YB"] == YB, 'Good fit'] = 'Data_NA'
            h_catalog.loc[h_catalog["YB"] == YB, 'Multiple emission line'] = 'Data_NA'
            
            h_catalog.to_csv("new.catalog.csv", index = False) #header = False) #pushing forward the data from h_catalog with the updates

        elif (GF == "" or pd.isnull(GF)) and (vlsr != 'Data_NA'):
            if YB_long >=15 and YB_long <41:
                #while (True):
                    print('YellowBall ID: ', YB)
               
                    
                    z = [(x)*step - 5.0 for x in range(nvel)]
                    y = (box_spec) #declaring z and y and bounding their std,mean and amplitude.
          #Bounding is helpful expecially if we have a wide peak or multiple peaks
              
                    def gaussian(x, a, mu, sig):
                        return a* (np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))
  
                    x_values = np.linspace(-5, 135, nvel)
                    #for mu, sig in [(-1, 1), (0, 2), (2, 3)]:
                    plt.plot(z,y, color='k')
  
                    plt.plot(z, gaussian(x_values, gamp, gvel, gstddev), color ='darkorange')
                      

                    plt.savefig(('TRIAL_amplitude is bounded and sd is 2/box_spec_%03i' % (YB)))
                    #plt.clf() 
                    
                    plt.show()                      
    
                    #check = input("Is this a good fit? Y/N: " )
                    check = input("hit Enter to skip if a good fit. Anything else for a bad fit: " )

                    #if check == 'y' or check =='Y' :
                    if check == "":
                        result = 'yes'
                    else:
                        result = 'no'
                    mel = input("Enter 'y' if there are multiple strong emission lines, if not hit 'enter' to coninue ")
                    if mel == 'y' or mel == 'Y':
                        mel_result = 'Yes'
                    else:
                        mel_result = 'No'
                  
                    plt.close()
            elif YB_long >=41 and YB_long <57.707:
                #while (True):
                    print('YellowBall ID: ', YB)
               
                    
                    z = [(x)*step - 5.0 for x in range(nvel)]
                    y = (box_spec) #declaring z and y and bounding their std,mean and amplitude.
          #Bounding is helpful expecially if we have a wide peak or multiple peaks
              
                    def gaussian(x, a, mu, sig):
                        return a* (np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))
  
                    x_values = np.linspace(-5, 85, nvel)
                    #for mu, sig in [(-1, 1), (0, 2), (2, 3)]:
                    plt.plot(z,y, color='k')
  
                    plt.plot(z, gaussian(x_values, gamp, gvel, gstddev), color ='darkorange')
                      

                    plt.savefig(('TRIAL_amplitude is bounded and sd is 2/box_spec_%03i' % (YB)))
                    #plt.clf() 
                    
                    plt.show()                      
    
                    #check = input("Is this a good fit? Y/N: " )
                    check = input("hit Enter to skip if a good fit. Anything else for a bad fit: " )

                    #if check == 'y' or check =='Y' :
                    if check == "":
                        result = 'yes'
                    else:
                        result = 'no'
                    mel = input("Enter 'y' if there are multiple strong emission lines, if not hit 'enter' to coninue ")
                    if mel == 'y' or mel == 'Y':
                        mel_result = 'Yes'
                    else:
                        mel_result = 'No'
                  
                    plt.close()    
        
            else: 
                   print('YellowBall ID: ', YB)
    
    
                  
                   z = [(x)*step - 200.0 for x in range(nvel)]
                   y = (box_spec) 
                   
                   def makegaussian(x,a,b,c):
                       output= a* (np.exp(1))**(-((x-b)**2)/((2*c)**2))
                       return output
                   
                   def gaussian(x, a, mu, sig):
                       return a* (np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))

                   x_values = np.linspace(-200, 200, nvel)
                   #for mu, sig in [(-1, 1), (0, 2), (2, 3)]:
                   plt.plot(z,y, color='k')

                   plt.plot(z, gaussian(x_values, gamp, gvel, gstddev), color ='darkorange')
                       

                   plt.savefig(('TRIAL_amplitude is bounded and sd is 2/box_spec_%03i' % (YB)))
                   #plt.clf() 
                   #plt.ion()
                   plt.show()
                       
                   #check = input("Is this a good fit? Y/N: " )
                   check = input("hit Enter to skip if a good fit. Anything else for a bad fit: " )

                   #if check == 'y' or check =='Y' :
                   if check == "":
                       result = 'yes'
                   else:
                       result = 'no'
                   mel = input("Enter 'y' if there are multiple strong emission lines, if not hit 'enter' to coninue ")
                   if mel == 'y' or mel == 'Y':
                       mel_result = 'Yes'
                   else:
                       mel_result = 'No'
                  
                   plt.close()
         

            h_catalog.loc[h_catalog["YB"] == YB, 'Good fit'] = result 
            h_catalog.loc[h_catalog["YB"] == YB, 'Multiple emission line'] =mel_result
            
            h_catalog.to_csv("new.catalog.csv", index = False) #header = False) #pushing forward the data from h_catalog with the updates

            
            ans = input("hit 'enter' to coninute, Enter 'n' to stop, enter 'a' if you would like to visit specific YB: ")
            if ans == 'n' or ans == 'N':
                break
            elif ans == 'a':
                YB_inpt = int(input("Enter the YB number you want to visit: "))
                begin ='a'
                i = 0
                if YB_inpt == 1153:
                    i=-1
               
    i +=1        
        


    
    
