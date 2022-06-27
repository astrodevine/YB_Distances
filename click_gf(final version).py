# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 15:00:41 2021

@author: Bezawit Mekasha Kassaye

What the program does: 
- Reads in CSV file with the guassian fit of the YBs informations 
- Finds the ones that were flagged with a 'No' for their goodfit.
- Display the spectra of these YBs and fit overplot for the source,
    asking for input of Gaussian parameters(intensity, cent vel, width) with an interactive clicking display.
- Refits, replots, asks Y/N if fit is good
- If N, redo again with inputs
-  If Y, overwrites in CSV the new fit parameters and a good fit (Y) flag

"""


# Reads in Hritik's output CSV file
# Skip to the first "N" (or bad fit) line
# Display the spectra and fit overplot for the source, ask for input of Gaussian parameters 
#intensity, cent vel, width)
# Refits, replots, asks Y/N if fit is good
# If N, redo again with inputs
# If Y, overwrites in CSV the new fit parameters and a good fit (Y) flag

#import cv2
import csv
#import glob
import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting

from astropy.io import fits
from astropy.wcs import WCS
from operator import add

import os.path


if os.path.exists('revisted(Good_fit)catalog.csv'): 
    begin = input("Do you want to start where you left off? enter 'no' if you want to star over. 'a' to vist a specific YB. Anything else to pick up where you left off.: ")
    if begin == 'no' or begin == 'No' or begin == 'NO':  
        h_catalog = pd.read_csv('new.catalog.csv')
    elif begin =='a' or begin =='A':
        h_catalog = pd.read_csv('new.catalog.csv')
        for i, row in h_catalog.iterrows(): 
            if h_catalog["Good fit"][i] == 'no':
                print("Here's the list of bad fits", h_catalog["YB"][i])
        YB_inpt = int(input("Enter the YB number you want to revist: ")) 
    else:
        h_catalog = pd.read_csv('revisted(Good_fit)catalog.csv')

else: 

    h_catalog = pd.read_csv("new.catalog.csv") #quoting=csv.QUOTE_NONE, error_bad_lines=False, header = None, names = ['YB','GLON',  'GLAT', 'r', 'Vstrongest (km/s)','g amplitude', 'g velocity','g stddev','uncertamp','uncertgvel','uncertstd', 'Good fit']) #modified from 'rb' to 'r' because of this error Error: iterator should return strings, not bytes (did you open the file in text mode?)

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
   
 

n = 0       #creating an index so that we can export the appended values of gaussian to the csv file


for i, row in h_catalog.iterrows():    
#get the YB's location and radius

    YB = int(h_catalog['YB'][i])
    YB_long = float(h_catalog['GLON'][i])
    YB_lat = float(h_catalog['GLAT'][i])
    YB_rad =float(h_catalog['r'][i])
    sdv = float(h_catalog['g stddev'][i])
    GF = h_catalog["Good fit"][i]
    
    
    
    
    if begin == 'a' or begin == 'A': 
        #print(YB, YB_inpt)

        if YB == YB_inpt and GF == 'no':
            
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
          
         
          
         
             if YB_long >=15 and YB_long <57.707:
                 while (True):
                     print(YB)
                     print(YB_long)
                     #print(nvel)
                     fig = plt.figure("Yellow Ball: %03i"%(YB))
                     fig.set_size_inches(15,8)
                     plt.xlabel('Velocity (km/s)')
                     plt.ylabel('Flux (k)')
                     plt.grid()
                     plt.plot([(x)*step - 5.0 for x in range(nvel)], box_spec) 
                     plt.xticks([-5, 0, 5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100, 105, 110, 115, 120, 125, 130, 135, 140])
                     plt.ylim(min(box_spec), max_flux+ (max_flux*50/100))

                     
                 #interactive clicking to fill up coords
                     fig.suptitle('(*Left)click to select    ' 
                                  '*Right click to undo selection    ' 
                                  '*when done click enter \n'
                                 '1) the peak point     ' 
                                 '2) on the min value for the width\n' 
                                 '3) on the max value for the width     '
                                 '4) on the max value for the amplitude\n' 
                                 '5) on the min value for the amplitude' 
                                 ,fontsize=10, fontweight ="bold",ha='left', va='top')
                    # coords = clkfig.ginput(n=-1, timeout=300, show_clicks=True, mouse_stop=2)
                     #x = fig.ginput(1)
                     coords = fig.ginput(n=-1, timeout=300, show_clicks=True, mouse_stop=2)
                     #plt.show()
                     plt.close()
                     #print(coords)

                     while len(coords) !=5:
                         print("Please try again. Select all the 5 points on the plot.")
                         coords =fig.ginput(n=-1, timeout=300, show_clicks=True, mouse_stop=2)
                         plt.close()
                     plt.close()
                     val_fv = coords[0]
                     max_flux = float(val_fv[1])
                     main_peak_vel = float(val_fv[0])
                     val_min_width = coords[1]
                     val_max_width = coords[2]
                     val_max_flux = coords[3]
                     val_min_flux = coords[4]
                     
                     
                     z = [(x)*step - 5.0 for x in range(nvel)]
                     y = (box_spec) 
                     
                     compound_model_bounded = models.Gaussian1D(max_flux,main_peak_vel, stddev= 2)
                     compound_model_bounded.stddev.bounds = (float(val_min_width[0])-5 , (float(val_max_width[0]))+5)
                     
                     compound_model_bounded.mean.bounds = (float(val_min_width[0]), (float(val_max_width[0])))
                     compound_model_bounded.amplitude.bounds = (float(val_min_flux[1]), float(val_max_flux[1]))
                                       
                                             
                     fitter = fitting.LevMarLSQFitter() #fitter to fit the data
                     compound_fit_bounded = fitter(compound_model_bounded,z, y) 
                     
                                     
                     gf_fig = plt.figure("Yellow Ball: %03i"%(YB))
                     plt.xlabel('velocity')
                     plt.ylabel('flux')
                     plt.plot(z,y, color='k')
                     plt.plot(z, compound_fit_bounded(z), color='darkorange')
         
                     gf_fig.suptitle("click 'Enter' to continue", fontsize =10, fontweight= 'bold')                                        
                     coords = gf_fig.ginput(n=0, timeout=300, show_clicks=False,mouse_stop=2)

                     trial = input("""Enter Y if you're happy with your fit. 
                                   Enter P to pass because you can't salvage it. 
                                   Anything else to try again: """)
                     #plt.close()
                     if trial == 'Y' or trial == 'y' or trial =='p' or trial == 'P':  
                          plt.savefig(('click_gf_trial(new)/box_spec_%03i' % (YB)))
                          break
                     else:
                         plt.clf()

                 
                 amp = (compound_fit_bounded.parameters[0])        #amplitude of the gaussian
                 gvel = (compound_fit_bounded.parameters[1])      #velocity of the gaussian
                 std = (compound_fit_bounded.parameters[2])      #standard deviation of the gaussian
                     
                 if (fitter.fit_info['param_cov']) is not None:
                     cov_diag = np.diag(fitter.fit_info['para%03im_cov'])
                 
                     uncertamp = (np.sqrt(cov_diag[0])) #uncertanity of amplitude
                     uncertgvel = (np.sqrt(cov_diag[1]))  #uncertanity of velocity
                     uncertstd = (np.sqrt(cov_diag[2])) #uncertanity of standard deviation
                 else:
                     #print('This is not good fit')  
                     uncertamp = ('bad fit')
                     uncertgvel =('bad fit')
                     uncertstd = ('bad fit')
      
         
             else: 
                 while (True):
                      fig = plt.figure("Yellow Ball: %03i"% (YB))
                      fig.set_size_inches(15,8)
                      plt.xlabel('Velocity (km/s)')
                      plt.ylabel('Flux (k)')
                      plt.grid()
                      plt.plot([(x)*step - 200.0 for x in range(nvel)], box_spec) #removed the +1 from x because the header was already skipped 
                      plt.ylim(min(box_spec), max_flux+ (max_flux*50/100))
                     
                     #plt.xticks([-5, 0, 5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90])
                      #plt.show()
                     
                      #plt.savefig('salvage_trial_graphs/box_spec_%03i' % (i[0])) #changed i +1 to i[0] so that the graph would associate with the YB index number.  
                      #plt.clf()
                      fig.suptitle('(*Left)click to select    ' 
                                   '*Right click to undo selection    ' 
                                   '*when done click enter \n'
                                  '1) the peak point     ' 
                                  '2) on the min value for the width\n' 
                                  '3) on the max value for the width     '
                                  '4) on the max value for the amplitude\n' 
                                  '5) on the min value for the amplitude' 
                                  ,fontsize=10, fontweight ="bold",ha='left', va='top')
                  
                      coords = fig.ginput(n=-1, timeout=300, show_clicks=True, mouse_stop=2)
                      plt.close()
                      
                      while len(coords) !=5:
                          print("Please try again. Select all the 5 points on the plot.")
                          coords =fig.ginput(n=-1, timeout=300, show_clicks=True, mouse_stop=2)
                          plt.close()
                      #print(coords)
                      plt.close()
                      
                      val_fv = coords[0]
                      max_flux = float(val_fv[1])
                      main_peak_vel = float(val_fv[0])
                     
                      val_min_width = coords[1]
                      val_max_width = coords[2]
                      val_max_flux = coords[3]
                      val_min_flux = coords[4]
                     
                     
                      z = [(x)*step - 200.0 for x in range(nvel)]
                      y = (box_spec) 
                     
                      compound_model_bounded = models.Gaussian1D(max_flux,main_peak_vel, stddev= 2)
                      compound_model_bounded.stddev.bounds = (float(val_min_width[0]) -5, (float(val_max_width[0])) +5)
                     
                      compound_model_bounded.mean.bounds = (float(val_min_width[0]), (float(val_max_width[0])))
                      compound_model_bounded.amplitude.bounds = (float(val_min_flux[1]), float(val_max_flux[1]))
                                       
                                             
                      fitter = fitting.LevMarLSQFitter() #fitter to fit the data
                      compound_fit_bounded = fitter(compound_model_bounded,z, y) 
                      gf_fig = plt.figure("Yellow Ball: %03i"%(YB))
                      plt.xlabel('velocity')
                      plt.ylabel('flux')
                      plt.plot(z,y, color='k')
                      plt.plot(z, compound_fit_bounded(z), color='darkorange')
         
                      gf_fig.suptitle("click 'Enter' to continue", fontsize =10, fontweight= 'bold')                                        
                      coords = gf_fig.ginput(n=0, timeout=300, show_clicks=False,mouse_stop=2)

                      trial = input("""Enter Y if you're happy with your fit. 
                                   Enter P to pass because you can't salvage it. 
                                   Anything else to try again: """)
                     #plt.close()
                      if trial == 'Y' or trial == 'y' or trial =='p' or trial == 'P':  
                          plt.savefig(('click_gf_trial(new)/box_spec_%03i' % (YB)))
                          plt.close()
                          break
                      else:
                          plt.clf()
                            
                  
                 amp = (compound_fit_bounded.parameters[0])        #appending the amplitude of the gaussian
                 gvel = (compound_fit_bounded.parameters[1])      #appending the velocity of the gaussian
                 std = (compound_fit_bounded.parameters[2])      #appending the standard deviation of the gaussian
                 
           
             #Getting the uncertanities of amplitude, velocity and standard deviation
             #some value will be none because they are bad fit and only have straight line
                 if (fitter.fit_info['param_cov']) is not None:
                     #print((fitter.fit_info['param_cov']))
                     cov_diag = np.diag(fitter.fit_info['param_cov'])
                 
                     uncertamp = (np.sqrt(cov_diag[0])) #uncertanity of amplitude
                     uncertgvel = (np.sqrt(cov_diag[1]))  #uncertanity of velocity
                     uncertstd = (np.sqrt(cov_diag[2])) #uncertanity of standard deviation
                 else:
                     uncertamp  = ('bad fit')
                     uncertgvel = ('bad fit')
                     uncertstd= ('bad fit')

           
             #overwriting the csv file with the updated values
             h_catalog.loc[h_catalog["YB"] == YB, 'Good fit'] = trial
             h_catalog.loc[h_catalog["YB"] == YB, 'g amplitude'] = amp
             h_catalog.loc[h_catalog["YB"] == YB, 'g velocity'] = gvel
             h_catalog.loc[h_catalog["YB"] == YB, 'g stddev'] = std
             h_catalog.loc[h_catalog["YB"] == YB, 'uncertamp'] = uncertamp
             h_catalog.loc[h_catalog["YB"] == YB, 'uncertgvel'] = uncertgvel
             h_catalog.loc[h_catalog["YB"] == YB, 'uncertstd'] = uncertstd   
            
             h_catalog.to_csv("revisted(Good_fit)catalog.csv", index = False) #pushing forward the data from h_catalog with the updates
       
             inpt_= input("Do you want to vist another YB? Enter 'y' to conintue. Anything else to stop revisiting:  ")
             if inpt_ == 'y' or inpt_ == 'Y':
                 begin = 'a'
                 YB_inpt = int(input("Enter the YB number you want to revist: "))
             else:
                 second_inpt = input("Do you want to continue working on the rest of the YBs? Enter 'n' to stop. Anything else to continue: ")
                 if second_inpt == 'n' or second_inpt == 'N':
                     break
                 else:
                     begin = 'yes'       
            
    else:        
        if GF =='no':
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
         
        
         
        
            if YB_long >=15 and YB_long <57.707:
                while (True):
                    print(YB)
                    print(YB_long)
                    #print(nvel)
                    fig = plt.figure("Yellow Ball: %03i"%(YB))
                    fig.set_size_inches(15,8)
                    plt.xlabel('Velocity (km/s)')
                    plt.ylabel('Flux (k)')
                    plt.grid()
                    plt.plot([(x)*step - 5.0 for x in range(nvel)], box_spec) 
                    plt.xticks([-5, 0, 5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100, 105, 110, 115, 120, 125, 130, 135, 140])
                    plt.ylim(min(box_spec), max_flux+ (max_flux*50/100))
    
                    
                #interactive clicking to fill up coords
                    fig.suptitle('(*Left)click to select    ' 
                                 '*Right click to undo selection    ' 
                                 '*when done click enter \n'
                                '1) the peak point     ' 
                                '2) on the min value for the width\n' 
                                '3) on the max value for the width     '
                                '4) on the max value for the amplitude\n' 
                                '5) on the min value for the amplitude' 
                                ,fontsize=10, fontweight ="bold",ha='left', va='top')
                   # coords = clkfig.ginput(n=-1, timeout=300, show_clicks=True, mouse_stop=2)
                    #x = fig.ginput(1)
                    coords = fig.ginput(n=-1, timeout=300, show_clicks=True, mouse_stop=2)
                    #plt.show()
                    plt.close()
                    #print(coords)
    
                    while len(coords) !=5:
                        print("Please try again. Select all the 5 points on the plot.")
                        coords =fig.ginput(n=-1, timeout=300, show_clicks=True, mouse_stop=2)
                        plt.close()
                    plt.close()
                    val_fv = coords[0]
                    max_flux = float(val_fv[1])
                    main_peak_vel = float(val_fv[0])
                    val_min_width = coords[1]
                    val_max_width = coords[2]
                    val_max_flux = coords[3]
                    val_min_flux = coords[4]
                    
                    
                    z = [(x)*step - 5.0 for x in range(nvel)]
                    y = (box_spec) 
                    
                    compound_model_bounded = models.Gaussian1D(max_flux,main_peak_vel, stddev= 2)
                    compound_model_bounded.stddev.bounds = (float(val_min_width[0])-5 , (float(val_max_width[0]))+5)
                    
                    compound_model_bounded.mean.bounds = (float(val_min_width[0]), (float(val_max_width[0])))
                    compound_model_bounded.amplitude.bounds = (float(val_min_flux[1]), float(val_max_flux[1]))
                                      
                                            
                    fitter = fitting.LevMarLSQFitter() #fitter to fit the data
                    compound_fit_bounded = fitter(compound_model_bounded,z, y) 
                    
                                    
                    gf_fig = plt.figure("Yellow Ball: %03i"%(YB))
                    plt.xlabel('velocity')
                    plt.ylabel('flux')
                    plt.plot(z,y, color='k')
                    plt.plot(z, compound_fit_bounded(z), color='darkorange')
        
                    gf_fig.suptitle("click 'Enter' to continue", fontsize =10, fontweight= 'bold')                                        
                    coords = gf_fig.ginput(n=0, timeout=300, show_clicks=False,mouse_stop=2)
    
                    trial = input("""Enter Y if you're happy with your fit. 
                                  Enter P to pass because you can't salvage it. 
                                  Anything else to try again: """)
                    #plt.close()
                    if trial == 'Y' or trial == 'y' or trial =='p' or trial == 'P':  
                         plt.savefig(('click_gf_trial(new)/box_spec_%03i' % (YB)))
                         break
                    else:
                        plt.clf()
    
                
                amp = (compound_fit_bounded.parameters[0])        #amplitude of the gaussian
                gvel = (compound_fit_bounded.parameters[1])      #velocity of the gaussian
                std = (compound_fit_bounded.parameters[2])      #standard deviation of the gaussian
                    
                if (fitter.fit_info['param_cov']) is not None:
                    cov_diag = np.diag(fitter.fit_info['para%03im_cov'])
                
                    uncertamp = (np.sqrt(cov_diag[0])) #uncertanity of amplitude
                    uncertgvel = (np.sqrt(cov_diag[1]))  #uncertanity of velocity
                    uncertstd = (np.sqrt(cov_diag[2])) #uncertanity of standard deviation
                else:
                    #print('This is not good fit')  
                    uncertamp = ('bad fit')
                    uncertgvel =('bad fit')
                    uncertstd = ('bad fit')
     
        
            else: 
                while (True):
                     fig = plt.figure("Yellow Ball: %03i"% (YB))
                     fig.set_size_inches(15,8)
                     plt.xlabel('Velocity (km/s)')
                     plt.ylabel('Flux (k)')
                     plt.grid()
                     plt.plot([(x)*step - 200.0 for x in range(nvel)], box_spec) #removed the +1 from x because the header was already skipped 
                     plt.ylim(min(box_spec), max_flux+ (max_flux*50/100))
                    
                    #plt.xticks([-5, 0, 5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90])
                     #plt.show()
                    
                     #plt.savefig('salvage_trial_graphs/box_spec_%03i' % (i[0])) #changed i +1 to i[0] so that the graph would associate with the YB index number.  
                     #plt.clf()
                     fig.suptitle('(*Left)click to select    ' 
                                  '*Right click to undo selection    ' 
                                  '*when done click enter \n'
                                 '1) the peak point     ' 
                                 '2) on the min value for the width\n' 
                                 '3) on the max value for the width     '
                                 '4) on the max value for the amplitude\n' 
                                 '5) on the min value for the amplitude' 
                                 ,fontsize=10, fontweight ="bold",ha='left', va='top')
                 
                     coords = fig.ginput(n=-1, timeout=300, show_clicks=True, mouse_stop=2)
                     plt.close()
                     
                     while len(coords) !=5:
                         print("Please try again. Select all the 5 points on the plot.")
                         coords =fig.ginput(n=-1, timeout=300, show_clicks=True, mouse_stop=2)
                         plt.close()
                     #print(coords)
                     plt.close()
                     
                     val_fv = coords[0]
                     max_flux = float(val_fv[1])
                     main_peak_vel = float(val_fv[0])
                    
                     val_min_width = coords[1]
                     val_max_width = coords[2]
                     val_max_flux = coords[3]
                     val_min_flux = coords[4]
                    
                    
                     z = [(x)*step - 200.0 for x in range(nvel)]
                     y = (box_spec) 
                    
                     compound_model_bounded = models.Gaussian1D(max_flux,main_peak_vel, stddev= 2)
                     compound_model_bounded.stddev.bounds = (float(val_min_width[0]) -5, (float(val_max_width[0])) +5)
                    
                     compound_model_bounded.mean.bounds = (float(val_min_width[0]), (float(val_max_width[0])))
                     compound_model_bounded.amplitude.bounds = (float(val_min_flux[1]), float(val_max_flux[1]))
                                      
                                            
                     fitter = fitting.LevMarLSQFitter() #fitter to fit the data
                     compound_fit_bounded = fitter(compound_model_bounded,z, y) 
                     gf_fig = plt.figure("Yellow Ball: %03i"%(YB))
                     plt.xlabel('velocity')
                     plt.ylabel('flux')
                     plt.plot(z,y, color='k')
                     plt.plot(z, compound_fit_bounded(z), color='darkorange')
        
                     gf_fig.suptitle("click 'Enter' to continue", fontsize =10, fontweight= 'bold')                                        
                     coords = gf_fig.ginput(n=0, timeout=300, show_clicks=False,mouse_stop=2)
    
                     trial = input("""Enter Y if you're happy with your fit. 
                                  Enter P to pass because you can't salvage it. 
                                  Anything else to try again: """)
                    #plt.close()
                     if trial == 'Y' or trial == 'y' or trial =='p' or trial == 'P':  
                         plt.savefig(('click_gf_trial(new)/box_spec_%03i' % (YB)))
                         plt.close()
                         break
                     else:
                         plt.clf()
                           
                 
                amp = (compound_fit_bounded.parameters[0])        #appending the amplitude of the gaussian
                gvel = (compound_fit_bounded.parameters[1])      #appending the velocity of the gaussian
                std = (compound_fit_bounded.parameters[2])      #appending the standard deviation of the gaussian
                
          
            #Getting the uncertanities of amplitude, velocity and standard deviation
            #some value will be none because they are bad fit and only have straight line
                if (fitter.fit_info['param_cov']) is not None:
                    #print((fitter.fit_info['param_cov']))
                    cov_diag = np.diag(fitter.fit_info['param_cov'])
                
                    uncertamp = (np.sqrt(cov_diag[0])) #uncertanity of amplitude
                    uncertgvel = (np.sqrt(cov_diag[1]))  #uncertanity of velocity
                    uncertstd = (np.sqrt(cov_diag[2])) #uncertanity of standard deviation
                else:
                    uncertamp  = ('bad fit')
                    uncertgvel = ('bad fit')
                    uncertstd= ('bad fit')
    
          
            #overwriting the csv file with the updated values
            h_catalog.loc[h_catalog["YB"] == YB, 'Good fit'] = trial
            h_catalog.loc[h_catalog["YB"] == YB, 'g amplitude'] = amp
            h_catalog.loc[h_catalog["YB"] == YB, 'g velocity'] = gvel
            h_catalog.loc[h_catalog["YB"] == YB, 'g stddev'] = std
            h_catalog.loc[h_catalog["YB"] == YB, 'uncertamp'] = uncertamp
            h_catalog.loc[h_catalog["YB"] == YB, 'uncertgvel'] = uncertgvel
            h_catalog.loc[h_catalog["YB"] == YB, 'uncertstd'] = uncertstd


            
    
        h_catalog.to_csv("revisted(Good_fit)catalog.csv", index = False) #pushing forward the data from h_catalog with the updates
        
        if GF =='no':
            continue_ = input("Click 'ENTER' to coninue or 'n' to stop.")
            if continue_ =='n' or continue_ == 'N':
                break
            else:
                begin = 'yes'
        
        