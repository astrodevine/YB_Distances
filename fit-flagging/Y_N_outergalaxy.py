#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 19:50:38 2022


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
#from the begining ornew.c visit a specific YB.
if os.path.exists('flagged(outer_galaxy).catalog.csv'): 
    begin = input("Enter 'restart' if you want to star over. Enter 'a' if you want to vist a specific YB. Anything else to start where you left off: ")
    if begin == 'Restart' or begin == 'restart' or begin == 'RESTART':  
        h_catalog = pd.read_csv('outergalaxy_vel_with_gf.csv') #should be vel_with_gf.csv
        h_catalog["Good fit"] = ""
        h_catalog ["Multiple emission line"] = ""
        YB_inpt = h_catalog['YB'][0]
    elif begin =='a' or begin =='A':
        h_catalog = pd.read_csv('lagged(outer_galaxy).catalog.csv')
        YB_inpt = int(input("Enter the YB number you want to revist: "))         
    else: 
        h_catalog = pd.read_csv('lagged(outer_galaxy).catalog.csv')
        YB_inpt = h_catalog['YB'][0]
else: 
    begin = input("Enter 'a' if you want to vist a YB, anything else to start from the beginning: ")
    if begin == 'a' or begin == 'A':
        h_catalog = pd.read_csv('outergalaxy_vel_with_gf.csv') # should be vel_with_gf.csv#quoting=csv.QUOTE_NONE, error_bad_lines=False, header = None, names = ['YB','GLON',  'GLAT', 'r', 'Vstrongest (km/s)','g amplitude', 'g velocity','g stddev','uncertamp','uncertgvel','uncertstd', 'Good fit']) #modified from 'rb' to 'r' because of this error Error: iterator should return strings, not bytes (did you open the file in text mode?)
        h_catalog["Good fit"] = ""
        h_catalog ["Multiple emission line"] = ""
        #h_catalog = pd.read_csv('new.catalog.csv')
        YB_inpt = int(input("Enter the YB number you want to revist: "))
    else:
        h_catalog = pd.read_csv('outergalaxy_vel_with_gf.csv') # should be vel_with_gf.csv#quoting=csv.QUOTE_NONE, error_bad_lines=False, header = None, names = ['YB','GLON',  'GLAT', 'r', 'Vstrongest (km/s)','g amplitude', 'g velocity','g stddev','uncertamp','uncertgvel','uncertstd', 'Good fit']) #modified from 'rb' to 'r' because of this error Error: iterator should return strings, not bytes (did you open the file in text mode?)
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
    
    #SOutergalaxy Data server
  cubeMF1 = fits.open("data_cubes/CGPS_MF1_CO_line_image.fits")[0] 
  cubeMF2 = fits.open("data_cubes/CGPS_MF2_CO_line_image.fits")[0]
  cubeMG1 = fits.open("data_cubes/CGPS_MG1_CO_line_image.fits")[0]
  cubeMG2 = fits.open('data_cubes/CGPS_MG2_CO_line_image.fits')[0]



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

  
    if YB_long >=106.2 and YB_long< 111.35 and YB_lat >= -3 and  YB_lat <= 1.55:
           yb=YellowBall(YB_long, YB_lat,YB_rad,cubeMF1)
           nvel = cubeMF1.header['naxis3'] -1
           
    elif YB_long >=106.2 and YB_long< 111.35 and  YB_lat >= 0.4 and  YB_lat <= 5.4:
           yb =YellowBall(YB_long, YB_lat,YB_rad,cubeMF2)
           nvel = cubeMF2.header['naxis3'] -1
            
    elif YB_long >=103 and YB_long< 107.3 and  YB_lat >= -3 and  YB_lat <= 1.55:
           yb =YellowBall(YB_long, YB_lat,YB_rad,cubeMG1)
           nvel = cubeMG1.header['naxis3'] -1
            
    elif YB_long >=103 and YB_long < 107.3 and  YB_lat >= 0.4 and  YB_lat <= 5.4:
           yb =YellowBall(YB_long, YB_lat,YB_rad,cubeMG2)
           nvel = cubeMG2.header['naxis3'] -1
    else:
          continue


    box_spec = yb.avg_spectrum()
    step= yb.cube.header['CDELT3']/1000
    main_peak_vel, max_flux = yb.get_main_peak()
     

  
    if begin == 'a' or begin == 'A': 
        
        if YB== YB_inpt:
            
                #while (True):
              print('YellowBall ID: ', YB)
             
             
              z = [(x)*step + 59.0 for x in range(nvel+1 )]
              y = (box_spec[0]) #declaring z and y and bounding their std,mean and amplitude.
            #Bounding is helpful expecially if we have a wide peak or multiple peaks
              def gaussian(x, a, mu, sig):
                  return a* (np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))
            
              x_values = np.linspace(59, -164, nvel+1)
              #for mu, sig in [(-1, 1), (0, 2), (2, 3)]:
              plt.plot(z,y, color='k')
            
              plt.plot(z, gaussian(x_values, gamp, gvel, gstddev), color ='darkorange')
              
              
              plt.savefig(('Outergalaxy_amplitude is bounded and sd is 2/box_spec_%03i' % (YB)))
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
           
         

              h_catalog.loc[h_catalog["YB"] == YB_inpt, 'Good fit'] = result 
              h_catalog.loc[h_catalog["YB"] == YB_inpt, 'Multiple emission line'] =mel_result
            
            
              h_catalog.to_csv("flagged(outer_galaxy).catalog.csv", index = False)  #pushing forward the data from h_catalog with the updates

            
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
            
            h_catalog.to_csv("flagged(outer_galaxy).catalog.csv", index = False) #header = False) #pushing forward the data from h_catalog with the updates

        elif (GF == "" or pd.isnull(GF)) and (vlsr != 'Data_NA'):
     
      #while (True):
          print('YellowBall ID: ', YB)
     
          
          z = [(x)*step +59 for x in range(nvel +1)]
          y = (box_spec[0]) #declaring z and y and bounding their std,mean and amplitude.
    #Bounding is helpful expecially if we have a wide peak or multiple peaks
    
          def gaussian(x, a, mu, sig):
              return a* (np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))
    
          x_values = np.linspace(59,-164,nvel+1)
          #for mu, sig in [(-1, 1), (0, 2), (2, 3)]:
          plt.plot(z,y, color='k')
    
          plt.plot(z, gaussian(x_values, gamp, gvel, gstddev), color ='darkorange')
            
    
          plt.savefig(('Outergalaxy_amplitude is bounded and sd is 2/box_spec_%03i' % (YB)))
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

      
     
        
           
        
        h_catalog.loc[h_catalog["YB"] == YB, 'Good fit'] = result 
        h_catalog.loc[h_catalog["YB"] == YB, 'Multiple emission line'] =mel_result
        
        h_catalog.to_csv("flagged(outer_galaxy).catalog.csv", index = False) #header = False) #pushing forward the data from h_catalog with the updates
        
        
        ans = input("hit 'enter' to coninute, Enter 'n' to stop, enter 'a' if you would like to visit specific YB: ")
        if ans == 'n' or ans == 'N':
            break
        elif ans == 'a':
            YB_inpt = int(input("Enter the YB number you want to visit: "))
            begin ='a'
            i = 0
            if YB_inpt == 3296:
                i=-1
           
    i +=1        
        


    
    