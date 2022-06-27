#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 08:24:41 2022

@author: bezawitmekashakassaye


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
      empty1, empty2, pixvel = self.wcs.all_world2pix(0.,vel*1000.,0)
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
      pixlon, pixlat, empty, empty= self.wcs.all_world2pix(self.glon,self.glat,0.,0, 0.,0,0)
      pixlon = int(round(pixlon))
      pixlat = int(round(pixlat))
      return self.data[:,:, pixlat+1,pixlon+1] 

#returns average spectrum of all pixels covered by YB extent, defined by central pixel +/- radius
   def avg_spectrum(self):
      b_pixlon, b_pixlat, empty, empty= self.wcs.all_world2pix(self.glon + self.r , self.glat - self.r, 0.,0.,0)
      e_pixlon, e_pixlat, empty, empty= self.wcs.all_world2pix(self.glon - self.r, self.glat + self.r, 0.,0,0)
      
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
      sum_spectrum = [0.]*(self.data.shape[1]-1) # The number of elements in the velocity dimention - 1.

      for pixel in box:
         spectrum = [x for x in self.data[:, :, pixel[1],pixel[0]]] 
         sum_spectrum = map(add, spectrum, sum_spectrum)

      avg_spectrum = [x/len(box) for x in sum_spectrum]
      
      return avg_spectrum
  
   def get_main_peak(self):
       
       
      avg_spectrum = self.avg_spectrum() #where the previous fun is being called. 
      #print("Here's avg_spectrum AGAIN", avg_spectrum)
     # max_flux = max(avg_spectrum)
     # pixvel = avg_spectrum.index(max_flux)
      A = avg_spectrum[0]
      #print(A)
      A[np.isnan(A)] = 0
      A = list(A)
      max_flux = max(A)
      # #max_flux = max(avg_spectrum)
      # print(max_flux)
      pixvel = A.index(max_flux)
      empty1, empty2, vel, empty3 = self.wcs.all_pix2world(0.,0.,pixvel,0.,0)
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
   cubeMF1 = fits.open("data_cubes/CGPS_MF1_CO_line_image.fits")[0] 
   cubeMF2 = fits.open("data_cubes/CGPS_MF2_CO_line_image.fits")[0]
   cubeMG1 = fits.open("data_cubes/CGPS_MG1_CO_line_image.fits")[0]
   cubeMG2 = fits.open('data_cubes/CGPS_MG2_CO_line_image.fits')[0]
    


   # for i in catalog:  
   #     i[1] = float(i[1])
   #     i[2]=float(i[2])
   #     i[3]=float(i[3])
   #     i[0]=int(i[0])
          
    
#      Hardcoding the standard deviations to grab standard deviation from noise free channels in specified YBs
#      3 selected values for each cube. 
       #NOTE: This section is hardcoded meaning, you need to do visiual inspection of your spectrum plots and change 
       #the YB ID and the bvel and evel values accordingly. 
      
      #CUBE 0
      
   
   
   # catalog_csv = open ("USE_THIS_CATALOG_ybcat_MWP_with_ID.csv","r") #modified from 'rb' to 'r' because of this error Error: iterator should return strings, not bytes (did you open the file in text mode?)
   # catalog = csv.reader(catalog_csv, delimiter=',', quotechar='"')
   # next(catalog, None)
   
   vel = open('outergalaxy_vel_with_gf.csv', 'w', newline='')
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
      
      
    
     
      if i[1] >=106.2 and i[1]< 111.35 and i[2] >= -3 and i[2] <= 1.55:
           yb=YellowBall(i[1],i[2],i[3],cubeMF1)
           nvel = cubeMF1.header['naxis3'] -1
           
      elif i[1] >=106.2 and i[1]< 111.35 and i[2] >= 0.4 and i[2] <= 5.4:
           yb =YellowBall(i[1],i[2],i[3],cubeMF2)
           nvel = cubeMF2.header['naxis3'] -1
            
      elif i[1] >=103 and i[1]< 107.3 and i[2] >= -3 and i[2] <= 1.55:
           yb =YellowBall(i[1],i[2],i[3],cubeMG1)
           nvel = cubeMG1.header['naxis3'] -1
            
      elif i[1] >=103 and i[1]< 107.3 and i[2] >= 0.4 and i[2] <= 5.4:
           yb =YellowBall(i[1],i[2],i[3],cubeMG2)
           nvel = cubeMG2.header['naxis3'] -1
      else:
          continue
         

          
      box_spec = yb.avg_spectrum()

      step= yb.cube.header['CDELT3']/1000
      main_peak_vel, max_flux = yb.get_main_peak()#where get_main_peak fun is being called.
      
      print(i[0], "vlsr", main_peak_vel)
          
                   
    
      plt.plot([(x)*step +59 for x in range(nvel +1 )], box_spec[0]) 
      plt.savefig('new_data_graph/box_spec_%03i' % (i[0])) #changed i +1 to i[0] so that the graph would associate with the YB index number.  
      plt.clf()
      z = [(x)*step +59 for x in range(nvel+1) ]
      y = (box_spec[0]) #declaring z and y and bounding their std,mean and amplitude.
      #   #Bounding is helpful expecially if we have a wide peak or multiple peaks
            
        
      compound_model_bounded = models.Gaussian1D(max_flux,main_peak_vel, stddev= 2) # try 0.2
      compound_model_bounded.stddev.bounds = (0,2)
      compound_model_bounded.mean.bounds = (main_peak_vel-5, main_peak_vel+5)
      compound_model_bounded.amplitude.bounds = (0.9*max_flux, 1.20*max_flux)
      
     #After bounding the parameters doing the gaussian fitting and exporting to a folder       

                
      fitter = fitting.LevMarLSQFitter()
      compound_fit_bounded = fitter(compound_model_bounded,z, y)
     
      plt.plot(z,y, color='k')
      plt.plot(z, compound_fit_bounded(z), color='darkorange')   
      plt.xlabel('velocity')
      plt.ylabel('flux')
      plt.savefig(('new_data_gf/box_spec_%03i' % (i[0])))
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
      # resid = []
      # for j in range(0, len(z)):
      #    resid.append( abs(compound_fit_bounded(z)[j] - y[j]))
      # plt.plot(z, resid, color='k')
      # plt.savefig(('new_data_residuals/box_spec_%03i' % (i[0])))
      # plt.clf()

    
      
      
      
      #if i[1] >= -1 and i[1]< 1  and i[2]>= -0.5 and i[2] <= 0.5 and max_flux > 4 * 0.2:        
      writer.writerow({'YB': i[0],'GLON': i[1], 'GLAT': i[2], 'r': i[3], 'Vstrongest (km/s)': main_peak_vel,'g amplitude': amp[index],'g velocity':gvel[index],'g stddev':std[index],'uncertamp':uncertamp[index],'uncertgvel':uncertgvel[index],'uncertstd':uncertstd[index],  '+/-' :5 ,'P(far)': 0.5, 'Extra': 'none' })
        
        
      index = index + 1

   vel.close()
   file = open('parameters from gaussian .csv', 'w', newline ='')
   with file:
       header = ['amplitude','velocity','standard deviation','uncertamp','uncertgvel','uncertstd']
       writer = csv.DictWriter(file, fieldnames = header)
       writer.writeheader()
       # print(compound_model_bounded.parameters[1])
       for i in range (0, len(amp)):
           writer.writerow({'amplitude' :(amp[i]),'velocity': gvel[i],'standard deviation': std[i],'uncertamp': uncertamp[i],'uncertgvel': uncertgvel[i], 'uncertstd': uncertstd[i]})

catalog_csv.close()