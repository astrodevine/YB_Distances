#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Jun 14 12:30:15 2021

Originally written by Leonard, modified for use and expanded by Bezawit Mekasha Kassaye.
Gausian fit section contributed by Hritik Rawat and merged with this code for use by Bezawit Mekasha Kassaye.

This code was written in Python3.6.

This code requires cube fits files that can be downloaded from SEDIGISM database and the BU GRS database. 

Links to download the cubes from:
       https://sedigism.mpifr-bonn.mpg.de/cgi-bin-seg/SEDIGISM_DATABASE.cgi
       http://grunt.bu.edu/grs-stitch/download-all.php

"""


import csv
import numpy as np 
import matplotlib.pyplot as plt
import math
import astropy
import astropy.modeling

from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model, Gaussian1D
from astropy.modeling.fitting import LevMarLSQFitter

from astropy.io import fits
from astropy.wcs import WCS
from operator import add


catalog_csv = open ("USE_THIS_CATALOG_ybcat_MWP_with_ID.csv","r") #modified from 'rb' to 'r' because of this error Error: iterator should return strings, not bytes (did you open the file in text mode?)
catalog = csv.reader(catalog_csv, delimiter=',', quotechar='"')
next(catalog, None) #skips the header
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
      plt.figure(0)

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
       print("these are the sd vals", sd_value)
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
   


   for i in catalog:  
       i[1] = float(i[1])
       i[2]=float(i[2])
       i[3]=float(i[3])
       i[0]=int(i[0])
          
    
#      Hardcoding the standard deviations to grab standard deviation from noise free channels in specified YBs
#      3 selected values for each cube. 
       #NOTE: This section is hardcoded meaning, you need to do visiual inspection of your spectrum plots and change 
       #the YB ID and the bvel and evel values accordingly. 
      
      #CUBE 0
       if i[0] == 1:
           yb=YellowBall(i[1],i[2],i[3], cube0)
           sd0_1 = float(yb.sd_fun(1,0,145,cube0))
       elif i[0]== 8:
           yb=YellowBall(i[1],i[2],i[3], cube0)
           sd0_2 = float(yb.sd_fun(8,250,395,cube0))
       elif i[0] ==9 :
           yb=YellowBall(i[1],i[2],i[3], cube0)
           sd0_3 = float(yb.sd_fun(9,250,395,cube0))
      #CUBE 2
       elif i[0] == 42:
           yb=YellowBall(i[1],i[2],i[3], cube02)
           sd2_1 = float(yb.sd_fun(42,0,100,cube02))
       elif i[0]== 46:
           yb=YellowBall(i[1],i[2],i[3], cube02)
           sd2_2 = float(yb.sd_fun(46,0,125,cube02))
       elif i[0] == 65 :
           yb=YellowBall(i[1],i[2],i[3], cube02)
           sd2_3 = float(yb.sd_fun(65,0,150,cube02))
      
        #CUBE 4
       elif i[0] == 86:
           yb=YellowBall(i[1],i[2],i[3], cube04)
           sd4_1 = float(yb.sd_fun(86,200,395,cube04))
       elif i[0]== 96:
           yb=YellowBall(i[1],i[2],i[3], cube04)
           sd4_2 = float(yb.sd_fun(96,250,395,cube04))
       elif i[0] ==110 :
           yb=YellowBall(i[1],i[2],i[3], cube04)
           sd4_3 = float(yb.sd_fun(110,250,395,cube04))
           
       #CUBE 6
       elif i[0] == 146:
           yb=YellowBall(i[1],i[2],i[3], cube06)
           sd6_1 = float(yb.sd_fun(146,0,145,cube06))
       elif i[0]== 148:
           yb=YellowBall(i[1],i[2],i[3], cube06)
           sd6_2 = float(yb.sd_fun(148, 0,150,cube06))
       elif i[0] ==152:
           yb=YellowBall(i[1],i[2],i[3], cube06)
           sd6_3 = float(yb.sd_fun(152,250,395,cube06))
           
       #CUBE 8
       
       elif i[0] == 213:
           yb=YellowBall(i[1],i[2],i[3], cube08)
           sd8_1 = float(yb.sd_fun(213,0,145,cube08))
       elif i[0]== 222:
           yb=YellowBall(i[1],i[2],i[3], cube08)
           sd8_2 = float(yb.sd_fun(222,0, 200,cube08))
       elif i[0] ==225:
           yb=YellowBall(i[1],i[2],i[3], cube08)
           sd8_3 = float(yb.sd_fun(225,0,175,cube08))
           
       #CUBE 10
       
       elif i[0] == 248:
           yb=YellowBall(i[1],i[2],i[3], cube10)
           sd10_1 = float(yb.sd_fun(248,0,150,cube10))
       elif i[0]== 287:
           yb=YellowBall(i[1],i[2],i[3], cube10)
           sd10_2 = float(yb.sd_fun(287,0,150,cube10))
       elif i[0] ==291 :
           yb=YellowBall(i[1],i[2],i[3], cube10)
           sd10_3 = float(yb.sd_fun(291,25,150,cube10))
           
       #CUBE 12
       
       elif i[0] == 361:
           yb=YellowBall(i[1],i[2],i[3], cube12)
           sd12_1 = float(yb.sd_fun(361,0,200,cube12))
       elif i[0]==378:
           yb=YellowBall(i[1],i[2],i[3], cube12)
           sd12_2 = float(yb.sd_fun(378,5,200,cube12))
       elif i[0] ==380 :
           yb=YellowBall(i[1],i[2],i[3], cube12)
           sd12_3 = float(yb.sd_fun(380,5,200,cube12))
           
       #CUBE 14
       
       elif i[0] == 397:
           yb=YellowBall(i[1],i[2],i[3], cube14)
           sd14_1 = float(yb.sd_fun(397,0,150,cube14))
       elif i[0]== 468:
           yb=YellowBall(i[1],i[2],i[3], cube14)
           sd14_2 = float(yb.sd_fun(468,250,395,cube14))
       elif i[0] ==480:
           yb=YellowBall(i[1],i[2],i[3], cube14)
           sd14_3 = float(yb.sd_fun(480,0,200,cube14))
        
        #CUBE 16
      
       if i[0] == 523:
           yb=YellowBall(i[1],i[2],i[3], cube16)
           sd1 = float(yb.sd_fun(523,60,130,cube16))
       elif i[0] == 525:
           yb=YellowBall(i[1],i[2],i[3], cube16)
           sd2 = float(yb.sd_fun(525,60,130,cube16))
       elif i[0] == 551:
           yb=YellowBall(i[1],i[2],i[3], cube16)
           sd3= float(yb.sd_fun(551,60,130,cube16))
     
        #CUBE 18
        
       elif i[0] == 590:
           yb=YellowBall(i[1],i[2],i[3], cube18)
           sd4 = float(yb.sd_fun(590,60,120,cube18))
       elif i[0]== 604:
           yb=YellowBall(i[1],i[2],i[3], cube18)
           sd5 = float(yb.sd_fun(604,80,130,cube18))
       elif i[0] == 608:
           yb=YellowBall(i[1],i[2],i[3], cube18)
           sd6 = float(yb.sd_fun(608,55,120,cube18))
       
        #CUBE 20
        
       elif i[0]== 678:
           yb=YellowBall(i[1],i[2],i[3], cube20)
           sd7 = float(yb.sd_fun(678,40,110,cube20))
       elif i[0] == 707:
           yb=YellowBall(i[1],i[2],i[3], cube20)
           sd8 = float(yb.sd_fun(707,60,130,cube20))
       elif i[0] == 713:
           yb=YellowBall(i[1],i[2],i[3], cube20)
           sd9 = float(yb.sd_fun(713,75,130,cube20))
       
       #CUBE22
       
       elif i[0]== 755:
           yb=YellowBall(i[1],i[2],i[3], cube22)
           sd10 = float(yb.sd_fun(755,90,130,cube22))
       elif i[0] == 785:
           yb=YellowBall(i[1],i[2],i[3], cube22)
           sd11 = float(yb.sd_fun(785,0,45,cube22))
       elif i[0]==805:
           yb=YellowBall(i[1],i[2],i[3], cube22)
           sd12 = float(yb.sd_fun(805,50,130,cube22))
           
      #CUBE 24           
      
       elif i[0]== 830:
           yb=YellowBall(i[1],i[2],i[3], cube24)
           sd13 = float(yb.sd_fun(830,0,45,cube24))
       elif i[0]== 893:
           yb=YellowBall(i[1],i[2],i[3], cube24)
           sd14 = float(yb.sd_fun(893,0,80,cube24))
       elif i[0]==895:
           yb=YellowBall(i[1],i[2],i[3], cube24)
           sd15 = float(yb.sd_fun(895,0,80,cube24))
           
       #CUBE 26   
       
       elif i[0]== 955:
           yb=YellowBall(i[1],i[2],i[3], cube26)
           sd16 = float(yb.sd_fun(955,0,65,cube26))
       elif i[0] == 976:
           yb=YellowBall(i[1],i[2],i[3], cube26)
           sd17 = float(yb.sd_fun(976,60,130,cube26))
       elif i[0] == 980:
           yb=YellowBall(i[1],i[2],i[3], cube26)
           sd18 = float(yb.sd_fun(980,0,80,cube26))
         
      #CUBE 28
      
       elif i[0] == 1040:
           yb=YellowBall(i[1],i[2],i[3], cube28)
           sd19 = float(yb.sd_fun(1040,60,130,cube28))
       elif i[0] == 1049:
           yb=YellowBall(i[1],i[2],i[3], cube28)
           sd20 = float(yb.sd_fun(1049,90,130,cube28))
       elif i[0]==1088:
           yb=YellowBall(i[1],i[2],i[3], cube28)
           sd21 = float(yb.sd_fun(1088,0,60,cube28))
           
           
           
       #CUBE 30  
       
       elif i[0] == 1166:
           yb=YellowBall(i[1],i[2],i[3], cube30)
           sd22 = float(yb.sd_fun(1166,0,60,cube30))
       elif i[0] == 1186:
           yb=YellowBall(i[1],i[2],i[3], cube30)
           sd23 = float(yb.sd_fun(1186,100,130,cube30))
       elif i[0]==1211:
            yb=YellowBall(i[1],i[2],i[3], cube30)
            sd24 = float(yb.sd_fun(1211,20,70,cube30))
       
       
       #CUBE 32
       
       elif i[0] == 1253:
           yb=YellowBall(i[1],i[2],i[3], cube32)
           sd25 = float(yb.sd_fun(1253,5,25,cube32))
       elif i[0] == 1267:
           yb=YellowBall(i[1],i[2],i[3], cube32)
           sd26= float(yb.sd_fun(1267, 0,35,cube32))
       elif i[0] ==1286:
            yb=YellowBall(i[1],i[2],i[3], cube32)
            sd27 = float(yb.sd_fun(1286,20,60,cube32))
            
       #CUBE 34     
       
       elif i[0] == 1358:
           yb=YellowBall(i[1],i[2],i[3], cube34)
           sd28= float(yb.sd_fun(1358, 20,50, cube34))
       elif i[0] == 1401:
           yb=YellowBall(i[1],i[2],i[3], cube34)
           sd29= float(yb.sd_fun(1401,65,135,cube34))
       elif i[0] == 1431:
            yb=YellowBall(i[1],i[2],i[3], cube34)
            sd30 = float(yb.sd_fun(1431,60,105,cube34))
            
            
       #CUBE 36   
       
       elif i[0] == 1458:
           yb=YellowBall(i[1],i[2],i[3], cube36)
           sd31= float(yb.sd_fun(1458,50,135,cube36))
       elif i[0] == 1498:
           yb=YellowBall(i[1],i[2],i[3], cube36)
           sd32= float(yb.sd_fun(1498,80,135,cube36))
       elif i[0]== 1486:
            yb=YellowBall(i[1],i[2],i[3], cube36)
            sd33 = float(yb.sd_fun(1486,100,130,cube36))
            
            
       #CUBE 38     
       
       elif i[0] == 1577:
           yb=YellowBall(i[1],i[2],i[3], cube38)
           sd34= float(yb.sd_fun(1577, 95, 135,cube38))
       elif i[0] == 1603:
           yb=YellowBall(i[1],i[2],i[3], cube38)
           sd35= float(yb.sd_fun(1603,35,70,cube38))
       elif i[0] == 1569:
           yb=YellowBall(i[1],i[2],i[3], cube38)
           sd36 = float(yb.sd_fun(1569,60,130,cube38))
           
           
       #CUBE 40   
       
       elif i[0] == 1618:  
           yb=YellowBall(i[1],i[2],i[3], cube40)
           sd37 = float(yb.sd_fun(1618,50,135,cube40))
       elif i[0] == 1663:
           yb=YellowBall(i[1],i[2],i[3], cube40)
           sd38 = float(yb.sd_fun(1663,75,135,cube40))   
       elif i[0] == 1688:
           yb=YellowBall(i[1],i[2],i[3], cube40)
           sd39 = float(yb.sd_fun(1688,60,130,cube40))
           
           
       #CUBE 42
       
       elif i[0] == 1753:
           yb=YellowBall(i[1],i[2],i[3], cube42)
           sd40=float(yb.sd_fun(1753, 0, 40, cube42))
       elif i[0] ==1740:
           yb=YellowBall(i[1],i[2],i[3], cube42)
           sd41=float(yb.sd_fun(1740,40,70, cube42))
       elif i[0]==1742:
           yb=YellowBall(i[1],i[2],i[3], cube42)
           sd42 = float(yb.sd_fun(1742,35,70,cube42))
           
       #CUBE 44
       
       elif i[0]== 1833:
           yb=YellowBall(i[1],i[2],i[3], cube44)
           sd43= float(yb.sd_fun(1833,0,35,cube44))
       elif i[0] == 1893:
           yb=YellowBall(i[1],i[2],i[3], cube44)
           sd44=float(yb.sd_fun(1893,0,50,cube44))
       elif i[0]==1854:
           yb=YellowBall(i[1],i[2],i[3], cube44)
           sd45 = float(yb.sd_fun(1854,0,45,cube44))
           
           
       #CUBE 46
       
       elif i[0]== 1895:
           yb=YellowBall(i[1],i[2],i[3], cube46)
           sd46= float(yb.sd_fun(1895, 0, 45, cube46))
       elif i[0]== 1942:
           yb=YellowBall(i[1],i[2],i[3], cube46)
           sd47= float(yb.sd_fun(1942, 0,40, cube46))
       elif i[0] == 1896:
           yb=YellowBall(i[1],i[2],i[3], cube46)
           sd48 = float(yb.sd_fun(1896,0,50,cube46))
           
           
       #CUBE 48
       
       elif i[0]== 2042:
           yb=YellowBall(i[1],i[2],i[3], cube48)
           sd49=float(yb.sd_fun(2042, 0, 40, cube48))
       elif i[0]==2028:
           yb=YellowBall(i[1],i[2],i[3], cube48)
           sd50=float(yb.sd_fun(2028, 60,80,cube48))
       elif i[0] == 2032:
           yb=YellowBall(i[1],i[2],i[3], cube48)
           sd51 = float(yb.sd_fun(2032,45,80,cube48))
           
           
       #CUBE 50
       
       elif i[0]==2050:
           yb=YellowBall(i[1],i[2],i[3], cube50)
           sd52 = float(yb.sd_fun(2050, 0, 40,cube50))
       elif i[0] == 2052:
           yb=YellowBall(i[1],i[2],i[3], cube50)
           sd53= float(yb.sd_fun(2052, 0,40,cube50))
       elif i[0] == 2091:
           yb=YellowBall(i[1],i[2],i[3], cube50)
           sd54 = float(yb.sd_fun(2091,0,35,cube50))
           
        #CUBE 52
        
       elif i[0] == 2188:
           yb=YellowBall(i[1],i[2],i[3], cube52)
           sd55 = float(yb.sd_fun(2188,50,80,cube52))
       elif i[0]==2202:
           yb=YellowBall(i[1],i[2],i[3], cube52)
           sd56 = float(yb.sd_fun(2202,0,40,cube52))
       elif i[0] == 2210:
           yb=YellowBall(i[1],i[2],i[3], cube52)
           sd57 = float(yb.sd_fun(2210,0,50,cube52))
           
           
       #CUBE 54
       
       elif i[0] ==2252:
           yb=YellowBall(i[1],i[2],i[3], cube54)
           sd58 = float(yb.sd_fun(2252,55,80,cube54))
       elif i[0]==2255:
           yb=YellowBall(i[1],i[2],i[3], cube54)
           sd59 = float(yb.sd_fun(2255,35,80,cube54))
       elif i[0] == 2282:
           yb=YellowBall(i[1],i[2],i[3], cube54)
           sd60 = float(yb.sd_fun(2282,50,80,cube54))
           
        
       #CUBE 56
           
       elif i[0] == 2309:
           yb=YellowBall(i[1],i[2],i[3], cube56)
           sd61 = float(yb.sd_fun(2309,55,80,cube56))
       elif i[0] == 2310:
           yb=YellowBall(i[1],i[2],i[3], cube56)
           sd62 = float(yb.sd_fun(2310,50,80,cube56))
       elif i[0] == 2313:
           yb=YellowBall(i[1],i[2],i[3], cube56)
           sd63 = float(yb.sd_fun(2313,10,80,cube56))
           
       #CUBE 302
       
       elif i[0] == 3525:
           yb=YellowBall(i[1],i[2],i[3], cube302)
           sd302_1 = float(yb.sd_fun(3525,200,395,cube302))
       elif i[0]== 3547:
           yb=YellowBall(i[1],i[2],i[3], cube302)
           sd302_2 = float(yb.sd_fun(3547,200,395,cube302))
       elif i[0] ==3549 :
           yb=YellowBall(i[1],i[2],i[3], cube302)
           sd302_3 = float(yb.sd_fun(3549,200,395,cube302))
           
       #CUBE 304
       
       elif i[0] == 3595:
           yb=YellowBall(i[1],i[2],i[3], cube304)
           sd304_1 = float(yb.sd_fun(3595,0,150,cube304))
       elif i[0]== 3599:
           yb=YellowBall(i[1],i[2],i[3], cube304)
           sd304_2 = float(yb.sd_fun(3599,5,200,cube304))
       elif i[0] ==3601 :
           yb=YellowBall(i[1],i[2],i[3], cube304)
           sd304_3 = float(yb.sd_fun(3601,250,395,cube304))
           
       #CUBE 306
       
       elif i[0] == 3678:
           yb=YellowBall(i[1],i[2],i[3], cube306)
           sd306_1 = float(yb.sd_fun(3678,200,395,cube306))
       elif i[0]== 3705:
           yb=YellowBall(i[1],i[2],i[3], cube306)
           sd306_2 = float(yb.sd_fun(3705,200,395,cube306))
       elif i[0] ==3729 :
           yb=YellowBall(i[1],i[2],i[3], cube306)
           sd306_3 = float(yb.sd_fun(3729,250,395,cube306))
       
       #CUBE 308
           
       elif i[0] == 3770:
           yb=YellowBall(i[1],i[2],i[3], cube308)
           sd308_1 = float(yb.sd_fun(3770,200,395,cube308))
       elif i[0]== 3778:
           yb=YellowBall(i[1],i[2],i[3], cube308)
           sd308_2 = float(yb.sd_fun(3778,200,395,cube308))
       elif i[0] ==3789 :
           yb=YellowBall(i[1],i[2],i[3], cube308)
           sd308_3 = float(yb.sd_fun(3789,250,395,cube308))
           
       #CUBE 310
       
       elif i[0] == 3834:
           yb=YellowBall(i[1],i[2],i[3], cube310)
           sd310_1 = float(yb.sd_fun(3834,200,395,cube310))
       elif i[0]== 3862:
           yb=YellowBall(i[1],i[2],i[3], cube310)
           sd310_2 = float(yb.sd_fun(3862,200,395,cube310))
       elif i[0] ==3885 :
           yb=YellowBall(i[1],i[2],i[3], cube310)
           sd310_3 = float(yb.sd_fun(3885,250,395,cube310))
       
        #CUBE 312
        
       elif i[0] == 3943:
           yb=YellowBall(i[1],i[2],i[3], cube312)
           sd312_1 = float(yb.sd_fun(3943,250,395,cube312))
       elif i[0]== 4005:
           yb=YellowBall(i[1],i[2],i[3], cube312)
           sd312_2 = float(yb.sd_fun(4005,200,395,cube312))
       elif i[0] ==4019 :
           yb=YellowBall(i[1],i[2],i[3], cube312)
           sd312_3 = float(yb.sd_fun(4019,0,125,cube312))
           
       #CUBE 314
       elif i[0] == 4086:
           yb=YellowBall(i[1],i[2],i[3], cube314)
           sd314_1 = float(yb.sd_fun(4086,200,395,cube314))
       elif i[0]== 4095:
           yb=YellowBall(i[1],i[2],i[3], cube314)
           sd314_2 = float(yb.sd_fun(4095,250,395,cube314))
       elif i[0] ==4104 :
           yb=YellowBall(i[1],i[2],i[3], cube314)
           sd314_3 = float(yb.sd_fun(4104,200,395,cube314))
       
      #CUBE 316
      
       elif i[0] == 4185:
           yb=YellowBall(i[1],i[2],i[3], cube316)
           sd316_1 = float(yb.sd_fun(4185,0,145,cube316))
       elif i[0]== 4199:
           yb=YellowBall(i[1],i[2],i[3], cube316)
           sd316_2 = float(yb.sd_fun(4199,200,395,cube316))
       elif i[0] ==4211:
           yb=YellowBall(i[1],i[2],i[3], cube316)
           sd316_3 = float(yb.sd_fun(4211,200,395,cube316))
           
       #CUBE 318
       
       elif i[0] == 4234:
           yb=YellowBall(i[1],i[2],i[3], cube318)
           sd318_1 = float(yb.sd_fun(4234,0,125,cube318))
       elif i[0]== 4275:
           yb=YellowBall(i[1],i[2],i[3], cube318)
           sd318_2 = float(yb.sd_fun(4275,250,395,cube318))
       elif i[0] ==4315:
           yb=YellowBall(i[1],i[2],i[3], cube318)
           sd318_3 = float(yb.sd_fun(4315,200,320,cube318))
           
       #CUBE 320
       
       elif i[0] == 4316:
           yb=YellowBall(i[1],i[2],i[3], cube320)
           sd320_1 = float(yb.sd_fun(4316,200,300,cube320))
       elif i[0]== 4327:
           yb=YellowBall(i[1],i[2],i[3], cube320)
           sd320_2 = float(yb.sd_fun(4327,250,395,cube320))
       elif i[0] ==4374:
           yb=YellowBall(i[1],i[2],i[3], cube320)
           sd320_3 = float(yb.sd_fun(4374,200,395,cube320))
           
       #CUBE 322
       
       elif i[0] == 4412:
           yb=YellowBall(i[1],i[2],i[3], cube322)
           sd322_1 = float(yb.sd_fun(4412,200,395,cube322))
       elif i[0]== 4443:
           yb=YellowBall(i[1],i[2],i[3], cube322)
           sd322_2 = float(yb.sd_fun(4443,200,395,cube322))
       elif i[0] ==4444:
           yb=YellowBall(i[1],i[2],i[3], cube322)
           sd322_3 = float(yb.sd_fun(4444,25,100,cube322))
           
       #CUBE 324
       
       elif i[0] == 4480:
           yb=YellowBall(i[1],i[2],i[3], cube324)
           sd324_1 = float(yb.sd_fun(4480,200,395,cube324))
       elif i[0]== 4487:
           yb=YellowBall(i[1],i[2],i[3], cube324)
           sd324_2 = float(yb.sd_fun(4487,200,395,cube324))
       elif i[0] ==4512:
           yb=YellowBall(i[1],i[2],i[3], cube324)
           sd324_3 = float(yb.sd_fun(4512,200,395,cube324))
           
       #CUBE 326
       
       elif i[0] == 4542:
           yb=YellowBall(i[1],i[2],i[3], cube326)
           sd326_1 = float(yb.sd_fun(4542,200,395,cube326))
       elif i[0]== 4575:
           yb=YellowBall(i[1],i[2],i[3], cube326)
           sd326_2 = float(yb.sd_fun(4575,200,395,cube326))
       elif i[0] ==4576:
           yb=YellowBall(i[1],i[2],i[3], cube326)
           sd326_3 = float(yb.sd_fun(4576,200,395,cube326))
           
       #CUBE 328
       
       elif i[0] == 4630:
           yb=YellowBall(i[1],i[2],i[3], cube328)
           sd328_1 = float(yb.sd_fun(4630,200,395,cube328))
       elif i[0]== 4672:
           yb=YellowBall(i[1],i[2],i[3], cube328)
           sd328_2 = float(yb.sd_fun(4672,200,395,cube328))
       elif i[0] ==4679:
           yb=YellowBall(i[1],i[2],i[3], cube328)
           sd328_3 = float(yb.sd_fun(4679,150,395,cube328))
           
       #CUBE 330
       
       elif i[0] == 4748:
           yb=YellowBall(i[1],i[2],i[3], cube330)
           sd330_1 = float(yb.sd_fun(4748,150,395,cube330))
       elif i[0]== 4755:
           yb=YellowBall(i[1],i[2],i[3], cube330)
           sd330_2 = float(yb.sd_fun(4755,200,395,cube330))
       elif i[0] ==4793:
           yb=YellowBall(i[1],i[2],i[3], cube330)
           sd330_3 = float(yb.sd_fun(4793,150,320,cube330))
           
       #CUBE 332
       
       elif i[0] == 4849:
           yb=YellowBall(i[1],i[2],i[3], cube332)
           sd332_1 = float(yb.sd_fun(4849,200,395,cube332))
       elif i[0]== 4860:
           yb=YellowBall(i[1],i[2],i[3], cube332)
           sd332_2 = float(yb.sd_fun(4860,200,395,cube332))
       elif i[0] ==4925:
           yb=YellowBall(i[1],i[2],i[3], cube332)
           sd332_3 = float(yb.sd_fun(4925,200,395,cube332))
           
       #CUBE 334
       
       elif i[0] == 4963:
           yb=YellowBall(i[1],i[2],i[3], cube334)
           sd334_1 = float(yb.sd_fun(4963,200,395,cube334))
       elif i[0]== 4982:
           yb=YellowBall(i[1],i[2],i[3], cube334)
           sd334_2 = float(yb.sd_fun(4982,200,395,cube334))
       elif i[0] ==5063:
           yb=YellowBall(i[1],i[2],i[3], cube334)
           sd334_3 = float(yb.sd_fun(5063,0,125,cube334))
           
       #CUBE 336
       
       elif i[0] == 5084:
           yb=YellowBall(i[1],i[2],i[3], cube336)
           sd336_1 = float(yb.sd_fun(5084,250,395,cube336))
       elif i[0]== 5156:
           yb=YellowBall(i[1],i[2],i[3], cube336)
           sd336_2 = float(yb.sd_fun(5156,150,395,cube336))
       elif i[0] ==5166:
           yb=YellowBall(i[1],i[2],i[3], cube336)
           sd336_3 = float(yb.sd_fun(5166,200,395,cube336))
           
       #CUBE 338
       
       elif i[0] == 5254:
           yb=YellowBall(i[1],i[2],i[3], cube338)
           sd338_1 = float(yb.sd_fun(5254,200,395,cube338))
       elif i[0]== 5302:
           yb=YellowBall(i[1],i[2],i[3], cube338)
           sd338_2 = float(yb.sd_fun(5302,150,395,cube338))
       elif i[0] ==5324:
           yb=YellowBall(i[1],i[2],i[3], cube338)
           sd338_3 = float(yb.sd_fun(5324,250,395,cube338))
           
       #CUBE 340
       
       elif i[0] == 5421:
           yb=YellowBall(i[1],i[2],i[3], cube340)
           sd340_1 = float(yb.sd_fun(5421,200,395,cube340))
       elif i[0]== 5422:
           yb=YellowBall(i[1],i[2],i[3], cube340)
           sd340_2 = float(yb.sd_fun(5422,200,395,cube340))
       elif i[0] ==5438 :
           yb=YellowBall(i[1],i[2],i[3], cube340)
           sd340_3 = float(yb.sd_fun(5438,200,395,cube340))
           
       #CUBE 342
       
       elif i[0] == 5480:
           yb=YellowBall(i[1],i[2],i[3], cube342)
           sd342_1 = float(yb.sd_fun(5480,200,395,cube342))
       elif i[0]== 5516:
           yb=YellowBall(i[1],i[2],i[3], cube342)
           sd342_2 = float(yb.sd_fun(5516,200,395,cube342))
       elif i[0] ==5528:
           yb=YellowBall(i[1],i[2],i[3], cube342)
           sd342_3 = float(yb.sd_fun(5528,200,395,cube342))
           
       #CUBE 344
       
       elif i[0] == 5563:
           yb=YellowBall(i[1],i[2],i[3], cube344)
           sd344_1 = float(yb.sd_fun(5563,250,395,cube344))
       elif i[0]== 5566:
           yb=YellowBall(i[1],i[2],i[3], cube344)
           sd344_2 = float(yb.sd_fun(5566,250,395,cube344))
       elif i[0] ==5604:
           yb=YellowBall(i[1],i[2],i[3], cube344)
           sd344_3 = float(yb.sd_fun(5604,200,395,cube344))
           
       #CUBE 346
       
       elif i[0] == 5652:
           yb=YellowBall(i[1],i[2],i[3], cube346)
           sd346_1 = float(yb.sd_fun(5652,250,395,cube346))
       elif i[0]== 5668:
           yb=YellowBall(i[1],i[2],i[3], cube346)
           sd346_2 = float(yb.sd_fun(5668,200,395,cube346))
       elif i[0] ==5676:
           yb=YellowBall(i[1],i[2],i[3], cube346)
           sd346_3 = float(yb.sd_fun(5676,250,395,cube346))
           
       #CUBE 348
       
       elif i[0] == 5757:
           yb=YellowBall(i[1],i[2],i[3], cube348)
           sd348_1 = float(yb.sd_fun(5757,250,395,cube348))
       elif i[0]== 5758:
           yb=YellowBall(i[1],i[2],i[3], cube348)
           sd348_2 = float(yb.sd_fun(5758,225,395,cube348))
       elif i[0] ==5783:
           yb=YellowBall(i[1],i[2],i[3], cube348)
           sd348_3 = float(yb.sd_fun(5783,0,75,cube348))
           
       #CUBE 350
       
       elif i[0] == 5823:
           yb=YellowBall(i[1],i[2],i[3], cube350)
           sd350_1 = float(yb.sd_fun(5823,200,395,cube350))
       elif i[0]== 5843:
           yb=YellowBall(i[1],i[2],i[3], cube350)
           sd350_2 = float(yb.sd_fun(5843,225,320,cube350))
       elif i[0] ==5839:
           yb=YellowBall(i[1],i[2],i[3], cube350)
           sd350_3 = float(yb.sd_fun(5839,0,125,cube350))
           
       #CUBE 352
       
       elif i[0] == 5863:
           yb=YellowBall(i[1],i[2],i[3], cube352)
           sd352_1 = float(yb.sd_fun(5863,0,150,cube352))
       elif i[0]== 5864:
           yb=YellowBall(i[1],i[2],i[3], cube352)
           sd352_2 = float(yb.sd_fun(5864,215,395,cube352))
       elif i[0] ==5918:
           yb=YellowBall(i[1],i[2],i[3], cube352)
           sd352_3 = float(yb.sd_fun(5918,150,395,cube352))
           
       #CUBE 354
       
       elif i[0] == 5946:
           yb=YellowBall(i[1],i[2],i[3], cube354)
           sd354_1 = float(yb.sd_fun(5946,0,150,cube354))
       elif i[0]== 5961:
           yb=YellowBall(i[1],i[2],i[3], cube352)
           sd354_2 = float(yb.sd_fun(5961,225,375,cube354))
       elif i[0] ==5955:
           yb=YellowBall(i[1],i[2],i[3], cube354)
           sd354_3 = float(yb.sd_fun(5955,200,395,cube354))
           
       #CUBE 356
       
       elif i[0] == 6008:
           yb=YellowBall(i[1],i[2],i[3], cube356)
           sd356_1 = float(yb.sd_fun(6008,0,150,cube356))
       elif i[0]== 6045:
           yb=YellowBall(i[1],i[2],i[3], cube356)
           sd356_2 = float(yb.sd_fun(6045,225,395,cube356))
       elif i[0] ==6046:
           yb=YellowBall(i[1],i[2],i[3], cube356)
           sd356_3 = float(yb.sd_fun(6046,0,150,cube356))
           
       #CUBE 358
       
       elif i[0] == 6093:
           yb=YellowBall(i[1],i[2],i[3], cube358)
           sd358_1 = float(yb.sd_fun(6093,0,150,cube358))
       elif i[0]== 6101:
           yb=YellowBall(i[1],i[2],i[3], cube358)
           sd358_2 = float(yb.sd_fun(6101,0,150,cube358))
       elif i[0] ==6104:
           yb=YellowBall(i[1],i[2],i[3], cube358)
           sd358_3 = float(yb.sd_fun(6104,230,395,cube358))
           
       #CUBE 359
       
       elif i[0] == 6147:
           yb=YellowBall(i[1],i[2],i[3], cube359)
           sd359_1 = float(yb.sd_fun(6147,0,150,cube359))
       elif i[0]== 6159:
           yb=YellowBall(i[1],i[2],i[3], cube359)
           sd359_2 = float(yb.sd_fun(6159,0,175,cube359))
       elif i[0] ==6165:
           yb=YellowBall(i[1],i[2],i[3], cube359)
           sd359_3 = float(yb.sd_fun(6165,0,175,cube359))
           
           

           
   #Noise calculation for inclusion in output table for each cube separately
   sd_c0 = np.mean([sd0_1, sd0_2, sd0_3])
   sd_c2 = np.mean([sd2_1, sd2_2, sd2_3])
   sd_c4 = np.mean([sd4_1, sd4_2, sd4_3])
   sd_c6 = np.mean([sd6_1, sd6_2, sd6_3])
   sd_c8 = np.mean([sd8_1, sd8_2, sd8_3])
   sd_c10 = np.mean([sd10_1, sd10_2, sd10_3])
   sd_c12 = np.mean([sd12_1, sd12_2, sd12_3])
   sd_c14 = np.mean([sd14_1, sd14_2, sd14_3])
   
   sd_c16 = np.mean([sd1, sd2, sd3])
   sd_c18 = np.mean([sd4, sd5, sd6])
   sd_c20 = np.mean([sd7, sd8, sd9])
   sd_c22 = np.mean([sd10, sd11, sd12])
   sd_c24 = np.mean([sd13, sd14, sd15])
   sd_c26 = np.mean([sd16, sd17, sd18])
   sd_c28 = np.mean([sd19, sd20, sd21])
   sd_c30 = np.mean([sd22, sd23, sd24])
   sd_c32 = np.mean([sd25, sd26, sd27])
   sd_c34 = np.mean([sd28, sd29, sd30])
   sd_c36 = np.mean([sd31, sd32, sd33])
   sd_c38 = np.mean([sd34, sd35, sd36])
   sd_c40 = np.mean([sd37, sd38, sd39])
   sd_c42 = np.mean([sd40, sd41, sd42])
   sd_c44 = np.mean([sd43, sd44, sd45])
   sd_c46 = np.mean([sd46, sd47, sd48])
   sd_c48 = np.mean([sd49, sd50, sd51])
   sd_c50 = np.mean([sd52, sd53, sd54])
   sd_c52 = np.mean([sd55, sd56, sd57])
   sd_c54 = np.mean([sd58, sd59, sd60])
   sd_c56 = np.mean([sd61, sd62, sd63])
   
   
   sd_c302 = np.mean([sd302_1, sd302_2, sd302_3])
   sd_c304 = np.mean([sd304_1, sd304_2,sd304_3])
   sd_c306 = np.mean([sd306_1, sd306_2, sd306_3])
   sd_c308 = np.mean([sd308_1, sd308_2, sd306_3])
   sd_c310 = np.mean([sd310_1, sd310_2, sd310_3])
   sd_c312 = np.mean([sd312_1, sd312_2, sd312_3])
   sd_c314 = np.mean([sd314_1, sd314_2, sd314_3])
   sd_c316 = np.mean([sd316_1, sd316_2, sd316_3])
   sd_c318 = np.mean([sd318_1, sd318_2, sd318_3])
   sd_c320 = np.mean([sd320_1, sd320_2, sd320_3])
   sd_c322 = np.mean([sd322_1, sd322_2, sd322_3])
   sd_c324 = np.mean([sd324_1, sd324_2, sd324_3])
   sd_c326 = np.mean([sd326_1, sd326_2, sd326_3])
   sd_c328 = np.mean([sd328_1, sd328_2, sd328_3])
   sd_c330 = np.mean([sd330_1, sd330_2, sd330_3])
   sd_c332 = np.mean([sd332_1, sd332_2, sd332_3])
   sd_c334 = np.mean([sd334_1, sd334_2, sd334_3])
   sd_c336 = np.mean([sd336_1, sd336_2, sd336_3])
   sd_c338 = np.mean([sd338_1, sd338_2, sd338_3])
   sd_c340 = np.mean([sd340_1, sd340_2, sd340_3])
   sd_c342 = np.mean([sd342_1, sd342_2, sd342_3])
   sd_c344 = np.mean([sd344_1, sd344_2, sd344_3])
   sd_c346 = np.mean([sd346_1, sd346_2, sd346_3])
   sd_c348 = np.mean([sd348_1, sd348_2, sd348_3])
   sd_c350 = np.mean([sd350_1, sd350_2, sd350_3])
   sd_c352 = np.mean([sd352_1, sd352_2, sd352_3])
   sd_c354 = np.mean([sd354_1, sd354_2, sd354_3])
   sd_c356 = np.mean([sd356_1, sd356_2, sd356_3])
   sd_c358 = np.mean([sd358_1, sd358_2, sd358_3])
   sd_c359 = np.mean([sd359_1, sd359_2, sd359_3])
   
   
   catalog_csv = open ("USE_THIS_CATALOG_ybcat_MWP_with_ID.csv","r") #modified from 'rb' to 'r' because of this error Error: iterator should return strings, not bytes (did you open the file in text mode?)
   catalog = csv.reader(catalog_csv, delimiter=',', quotechar='"')
   next(catalog, None)
   
   vel = open('vel_with_gf.csv', 'w', newline='')
   field_names = ['YB','GLON', 'GLAT', 'r', 'Vstrongest (km/s)','g amplitude','g velocity','g stddev','uncertamp','uncertgvel','uncertstd',  '+/-',  'P(far)' , 'Extra']
   writer = csv.DictWriter(vel, fieldnames=field_names)
   writer.writeheader() 
   
   index = 0       #creating an index so that we can export the appended values of gaussian to the csv file
   amp = []        #stores an array of amplitude of the gaussian fit
   gvel = []       #stores an array of gaussian velocity
   std = []        #stores an array of standard deviation for gaussian fit 
   uncertamp = []  #stores an array of uncertanity of amplitude
   uncertgvel = [] #stores an array of uncertanity of gaussian velocity
   uncertstd = []  #stores an array of uncertanity of standard deviation
   
   for i in catalog: 
      i[1] = float(i[1]) #glon value
      i[2]=float(i[2]) #glat value
      i[3]=float(i[3]) #radius value
      i[0]=int(i[0]) #YB ID 
      
      
      #SEDGISM dataset 
      
      if i[1]>= -1 and i[1]< 1 and i[2] >= -0.5 and i[2] <= 0.5: #based on the cube's longitude and lattitude range. 
          yb=YellowBall(i[1], i[2], i[3], cube0) #Calling the class depending on which cube the YBs fit in
          nvel = cube0.header['naxis3'] -1 #nvel is the number of velocity channels in the cube, and the information is extracted from the header of the cube fits files. 
      elif i[1] >=1 and i[1]< 3 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube02)
           nvel = cube02.header['naxis3'] -1
      elif i[1] >=3 and i[1]< 5 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube04)
           nvel = cube04.header['naxis3'] -1
      elif i[1] >=5 and i[1]< 7  and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube06)
           nvel = cube06.header['naxis3'] -1
      elif i[1] >=7 and i[1]< 9 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube08)
           nvel = cube08.header['naxis3'] -1
      elif i[1] >=9 and i[1]< 11 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube10)
           nvel = cube10.header['naxis3'] -1
      elif i[1] >=11 and i[1]< 13 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube12)
           nvel = cube12.header['naxis3'] -1
      elif i[1] >=13 and i[1]< 15 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube14)
           nvel = cube14.header['naxis3'] -1
           
           
      #BU dataset     
      elif i[1] >=15 and i[1]< 17 and i[2]>= -0.806 and i[2] <0.818: 
           yb=YellowBall(i[1],i[2],i[3],cube16)
           nvel = cube16.header['naxis3'] -1
      elif i[1] >=17 and i[1]< 19:
           yb=YellowBall(i[1],i[2],i[3],cube18)
           nvel = cube18.header['naxis3'] -1
      elif i[1] >=19 and i[1]< 21:
           yb=YellowBall(i[1],i[2],i[3],cube20)
           nvel = cube20.header['naxis3'] -1
      elif i[1] >=21 and i[1]< 23:
           yb=YellowBall(i[1],i[2],i[3],cube22)
           nvel = cube22.header['naxis3'] -1
      elif i[1] >=23 and i[1]< 25:
           yb=YellowBall(i[1],i[2],i[3],cube24)
           nvel = cube24.header['naxis3'] -1
      elif i[1] >=25 and i[1]< 27:
           yb=YellowBall(i[1],i[2],i[3],cube26)
           nvel = cube26.header['naxis3'] -1
      elif i[1] >=27 and i[1]< 29:
           yb=YellowBall(i[1],i[2],i[3],cube28)
           nvel = cube28.header['naxis3'] -1 
      elif i[1] >=29 and i[1]< 31:
           yb=YellowBall(i[1],i[2],i[3],cube30)
           nvel = cube30.header['naxis3'] -1          
      elif i[1] >=31 and i[1]< 33:
           yb=YellowBall(i[1],i[2],i[3],cube32)
           nvel = cube32.header['naxis3'] -1
      elif i[1]>=33 and i[1] <35:
           yb=YellowBall(i[1],i[2],i[3],cube34)
           nvel = cube34.header['naxis3'] -1
      elif i[1]>=35 and i[1] <37:
           yb=YellowBall(i[1],i[2],i[3], cube36)
           nvel = cube36.header['naxis3'] -1
      elif i[1]>=37 and i[1] <39:
           yb=YellowBall(i[1],i[2],i[3], cube38)
           nvel = cube38.header['naxis3'] -1
      elif i[1] >=39 and i[1] <41:
           yb=YellowBall(i[1],i[2],i[3],cube40)
           nvel = cube40.header['naxis3'] -1
      elif i[1] >=41 and i[1] <43:
           yb=YellowBall(i[1],i[2],i[3],cube42)
           nvel = cube42.header['naxis3'] -1
      elif i[1] >=43 and i[1] <45:
           yb=YellowBall(i[1],i[2],i[3],cube44)
           nvel = cube44.header['naxis3'] -1
      elif i[1] >=45 and i[1] <47:
           yb=YellowBall(i[1],i[2],i[3],cube46)
           nvel = cube46.header['naxis3'] -1
      elif i[1] >=47 and i[1] <49:
           yb=YellowBall(i[1],i[2],i[3],cube48)
           nvel = cube48.header['naxis3'] -1
   
      elif i[1] >=49 and i[1]< 51:
           yb=YellowBall(i[1],i[2],i[3],cube50)
           nvel = cube50.header['naxis3'] -1
      elif i[1] >=51 and i[1]< 53:
           yb=YellowBall(i[1],i[2],i[3],cube52)
           nvel = cube52.header['naxis3'] -1
      elif i[1] >=53 and i[1]< 55:
           yb=YellowBall(i[1],i[2],i[3],cube54)
           nvel = cube54.header['naxis3'] -1
      elif i[1] >=55 and i[1]< 55.707 and i[2] >= -1.095 and i[2]<1.095:
          yb=YellowBall(i[1],i[2],i[3],cube56)
          nvel = cube56.header['naxis3'] -1
         
 #SEDISIM data set 
      elif i[1] >=301 and i[1]<303 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube302)
           nvel = cube302.header['naxis3'] -1
      elif i[1] >=303 and i[1]<305  and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube304)
           nvel = cube304.header['naxis3'] -1
      elif i[1] >=305 and i[1]< 307 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube306)
           nvel = cube306.header['naxis3'] -1
      elif i[1] >=307 and i[1]< 309 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube308)
           nvel = cube308.header['naxis3'] -1
      elif i[1] >=309 and i[1]< 311 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube310)
           nvel = cube310.header['naxis3'] -1
      elif i[1] >=311 and i[1]< 313 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube312)
           nvel = cube312.header['naxis3'] -1
      elif i[1] >=313 and i[1]< 315 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube314)
           nvel = cube314.header['naxis3'] -1
      elif i[1] >=315 and i[1]< 317 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube316)
           nvel = cube316.header['naxis3'] -1
      elif i[1] >=317 and i[1]< 319 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube318)
           nvel = cube318.header['naxis3'] -1
      elif i[1] >=319 and i[1]< 321 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube320)
           nvel = cube320.header['naxis3'] -1
      elif i[1] >=321 and i[1]< 323 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube322)
           nvel = cube322.header['naxis3'] -1
      elif i[1] >=323 and i[1]< 325 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube324)
           nvel = cube324.header['naxis3'] -1
      elif i[1] >=325 and i[1]< 327 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube326)
           nvel = cube326.header['naxis3'] -1
      elif i[1] >=327 and i[1]< 329 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube328)
           nvel = cube328.header['naxis3'] -1
      elif i[1] >=329 and i[1]< 331 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube330)
           nvel = cube330.header['naxis3'] -1    
      elif i[1] >=331 and i[1]< 333 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube332)
           nvel = cube332.header['naxis3'] -1
      elif i[1] >=333 and i[1]< 335 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube334)
           nvel = cube334.header['naxis3'] -1
      elif i[1] >=335 and i[1]< 337 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube336)
           nvel = cube336.header['naxis3'] -1
      elif i[1] >=337 and i[1]< 339 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube338)
           nvel = cube338.header['naxis3'] -1
      elif i[1] >=339 and i[1]< 341 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube340)
           nvel = cube340.header['naxis3'] -1
      elif i[1] >=341 and i[1]< 343 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube342)
           nvel = cube342.header['naxis3'] -1
      elif i[1] >=343 and i[1]< 345 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube344)
           nvel = cube344.header['naxis3'] -1
      elif i[1] >= 345 and i[1]< 347 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube346)
           nvel = cube346.header['naxis3'] -1
      elif i[1] >=347 and i[1]< 349 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube348)
           nvel = cube348.header['naxis3'] -1
      elif i[1] >=349 and i[1]< 351 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube350)
           nvel = cube350.header['naxis3'] -1
      elif i[1] >=351 and i[1]< 353 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube352)
           nvel = cube352.header['naxis3'] -1
      elif i[1] >=353 and i[1]< 355 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube354)
           nvel = cube354.header['naxis3'] -1
      elif i[1] >=355 and i[1]< 357 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube356)
           nvel = cube356.header['naxis3'] -1
      elif i[1] >=357 and i[1]< 359 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube358)
           nvel = cube358.header['naxis3'] -1    
      elif i[1] >=359 and i[1]< 360 and i[2] >= -0.5 and i[2] <= 0.5:
           yb=YellowBall(i[1],i[2],i[3],cube359)
           nvel = cube359.header['naxis3'] -1
           


      plt.figure(0)
      plt.xlabel('Velocity (km/s)')
      plt.ylabel('Flux (k)')
      plt.grid()
          
      box_spec = yb.avg_spectrum()

      step= yb.cube.header['CDELT3']/1000
      main_peak_vel, max_flux = yb.get_main_peak()#where get_main_peak fun is being called.
      
      
      if i[1] >=15 and i[1] <57.707:
          plt.plot([(x)*step - 5.0 for x in range(nvel)], box_spec) #removed the +1 from x because the header was already skipped 
          plt.savefig('66percent_YB_graphs/box_spec_%03i' % (i[0])) #changed i +1 to i[0] so that the graph would associate with the YB index number.  
          plt.clf()
          
          z = [(x)*step - 5.0 for x in range(nvel)]
          y = (box_spec) #declaring z and y and bounding their std,mean and amplitude.
      #Bounding is helpful expecially if we have a wide peak or multiple peaks
          
          compound_model_bounded = models.Gaussian1D(max_flux,main_peak_vel, stddev= 2)
          compound_model_bounded.stddev.bounds = (0,2)
          compound_model_bounded.mean.bounds = (main_peak_vel-5, main_peak_vel+5)
          compound_model_bounded.amplitude.bounds = (0.9*max_flux, 1.20*max_flux)
  
    #After bounding the parameters doing the gaussian fitting and exporting to a folder       
   
          fitter = fitting.LevMarLSQFitter() #fitter to fit the data
          compound_fit_bounded = fitter(compound_model_bounded,z, y) 
          plt.plot(z,y, color='k')
          plt.plot(z, compound_fit_bounded(z), color='darkorange')   
          plt.xlabel('velocity')
          plt.ylabel('flux')
          plt.savefig(('amplitude is bounded and sd is 2(new)/box_spec_%03i' % (i[0])))
          plt.clf()   
          
          
          
      else:
          plt.plot([(x)*step - 200.0 for x in range(nvel)], box_spec) #removed the +1 from x because the header was already skipped 
          plt.savefig('66percent_YB_graphs/box_spec_%03i' % (i[0])) #changed i +1 to i[0] so that the graph would associate with the YB index number.  
          plt.clf()
          z = [(x)*step - 200.0 for x in range(nvel)]
          y = (box_spec)
          
          compound_model_bounded = models.Gaussian1D(max_flux,main_peak_vel, stddev= 2)
          compound_model_bounded.stddev.bounds = (0,2)
          compound_model_bounded.mean.bounds = (main_peak_vel-5, main_peak_vel+5)
          compound_model_bounded.amplitude.bounds = (0.9*max_flux, 1.20*max_flux)
      
                    
          fitter = fitting.LevMarLSQFitter()
          compound_fit_bounded = fitter(compound_model_bounded,z, y)
          plt.plot(z,y, color='k')
          plt.plot(z, compound_fit_bounded(z), color='darkorange')   
          plt.xlabel('velocity')
          plt.ylabel('flux')
          plt.savefig(('amplitude is bounded and sd is 2(new)/box_spec_%03i' % (i[0])))
          plt.clf()
   
     
      amp.append(compound_fit_bounded.parameters[0])        #appending the amplitude of the gaussian
      gvel.append(compound_fit_bounded.parameters[1])      #appending the velocity of the gaussian
      std.append(compound_fit_bounded.parameters[2])      #appending the standard deviation of the gaussian
      
      #Getting the uncertanities of amplitude, velocity and standard deviation
      #some value will be none because they are bad fit and only have straight line
      if (fitter.fit_info['param_cov']) is not None:
         print((fitter.fit_info['param_cov']))
         cov_diag = np.diag(fitter.fit_info['param_cov'])

         uncertamp.append(np.sqrt(cov_diag[0])) #uncertanity of amplitude
         uncertgvel.append(np.sqrt(cov_diag[1]))  #uncertanity of velocity
         uncertstd.append(np.sqrt(cov_diag[2])) #uncertanity of standard deviation
       
      else:
          print('This is not good fit')  
          uncertamp.append('bad fit')
          uncertgvel.append('bad fit')
          uncertstd.append('bad fit')
          
      #calculating the residual and exporting their plot to a separate folder
      resid = []
      for j in range(0, len(z)):
         resid.append( abs(compound_fit_bounded(z)[j] - y[j]))
      plt.plot(z, resid, color='k')
      plt.savefig(('residuals when sd is 2/box_spec_%03i' % (i[0])))
      plt.clf()

      
      
      
      
      if i[1] >= -1 and i[1]< 1  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c0:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel,'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],  '+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
        
      elif i[1] >=1 and i[1]< 3  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c2:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none'})
           
      elif i[1] >= 3 and i[1]< 5  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c4:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 5 and i[1]< 7  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c6:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 7 and i[1]< 9  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c8:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 9 and i[1]< 11  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c10:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 11 and i[1]< 13  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c12:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 13 and i[1]< 15  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c14:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
        
        
      elif i[1] >=15 and i[1]< 17  and i[2]>= -0.806 and i[2] <0.818 and max_flux > 4 * sd_c16:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      

      elif i[1] >=17 and i[1]< 19 and max_flux > 4*sd_c18:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel , 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none'})

      elif i[1] >=19 and i[1]< 21 and max_flux > 4*sd_c20:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel , 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none'})

          
      elif i[1] >=21 and i[1]< 23 and max_flux > 4*sd_c22:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index] ,'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none'})

          
      elif i[1] >=23 and i[1]< 25 and max_flux> 4*sd_c24:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

          
      elif i[1] >=25 and i[1]< 27 and max_flux > 4*sd_c26:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

           
      elif i[1] >=27 and i[1]< 29 and max_flux > 4*sd_c28:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

           
      elif i[1] >=29 and i[1]< 31 and max_flux > 4*sd_c30:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

      elif i[1] >=31 and i[1]< 33 and max_flux > 4*sd_c32:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

      elif i[1]>=33 and i[1] <35 and max_flux > 4*sd_c34:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

          
      elif i[1]>=35 and i[1] <37 and max_flux > 4*sd_c36:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

      
      elif i[1]>=37 and i[1] <39 and max_flux > 4*sd_c38:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

           
      elif i[1] >=39 and i[1] <41 and max_flux > 4*sd_c40:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

        
      elif i[1] >=41 and i[1] <43 and max_flux > 4*sd_c42:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

      elif i[1] >=43 and i[1] <45 and max_flux > 4*sd_c44:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

          
      elif i[1] >=45 and i[1] <47 and max_flux > 4*sd_c46:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

           
      elif i[1] >=47 and i[1] <49 and max_flux > 4*sd_c48:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

   
      elif i[1] >=49 and i[1]< 51 and max_flux > 4*sd_c50:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })


      elif i[1] >=51 and i[1]< 53 and max_flux > 4*sd_c52:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

      elif i[1] >=53 and i[1]< 55 and max_flux > 4*sd_c54:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

          
      elif i[1] >=55 and i[1]< 55.707 and max_flux > 4*sd_c56:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })

     
        
      elif i[1] >= 301 and i[1]< 303  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c302:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 303 and i[1]< 305  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c304:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 305 and i[1]< 307  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c306:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 307 and i[1]< 309  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c308:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 309 and i[1]< 311  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c310:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 311 and i[1]< 313  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c312:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 313 and i[1]< 315  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c314:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 315 and i[1]< 317  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c316:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 317 and i[1]< 319  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c318:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 319 and i[1]< 321  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c320:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 321 and i[1]< 323  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c322:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 323 and i[1]< 325  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c324:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 325 and i[1]< 327  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c326:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 327 and i[1]< 329  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c328:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 329 and i[1]< 331  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c330:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 331 and i[1]< 333  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c332:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 333 and i[1]< 335  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c334:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 335 and i[1]< 337  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c336:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 337 and i[1]< 339  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c338:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 339 and i[1]< 341  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c340:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 341 and i[1]< 343  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c342:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 343 and i[1]< 345  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c344:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 345 and i[1]< 347  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c346:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 347 and i[1]< 349  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c348:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 349 and i[1]< 351  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c350:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 351 and i[1]< 353  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c352:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 353 and i[1]< 355  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c354:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 355 and i[1]< 357  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c356:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 357 and i[1]< 359  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c358:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
      elif i[1] >= 359 and i[1]< 360  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * sd_c359:        
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel, 'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
      
  
      elif i[1] >= -1 and i[1] < 15 and i[2] >0.5 or i[2] <-0.5:
           writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': 'Data_NA','g amplitude': amp[index],'g velocity':'Data_NA','g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none'}) 
           
      elif i[1] >= 301 and i[1] <360 and i[2] >0.5 or i[2] <-0.5:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': 'Data_NA','g amplitude': amp[index],'g velocity':'Data_NA','g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none'})   
          
      elif i[1] >=15 and i[1]< 17  and i[2]< -0.806 or i[2] >= 0.818:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': 'Data_NA','g amplitude': amp[index],'g velocity':'Data_NA','g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none'})    
          
      else:
          writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': 'NA','g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],'+/-' :5 ,'P(far)': 0.5, 'Extra': 'none'}) 
          
      index = index + 1

   vel.close()
#   file = open('parameters from gaussian .csv', 'w', newline ='')
#   with file:
#      header = ['amplitude','velocity','standard deviation','uncertamp','uncertgvel','uncertstd']
#      writer = csv.DictWriter(file, fieldnames = header)
#      writer.writeheader()
#      # print(compound_model_bounded.parameters[1])
#      for i in range (0, len(amp)):
#           writer.writerow({'amplitude' :(amp[i]),'velocity': gvel[i],'standard deviation': std[i],'uncertamp': uncertamp[i],'uncertgvel': uncertgvel[i], 'uncertstd': uncertstd[i]})
# 
catalog_csv.close()
