


""" 
Written by Leonardo Trujillo. 08-03-2017
Modified for use by Bezawit Mekasha Kassaye and Hritik Rawat. 06-09-2021, 
Expanded to 66% of the catalogue by Bezawit Mekasha Kassaye. 07-27-2021.

This program runs the Bayesian Distance Calculator written by Reid et al. on a
list of sources provided on a csv file.

To see instructions on how to use this program, refer to the README.txt file.

This program was modified to run the second version of the Bayseian distance calculator code written by Reid et al in 2019. 
"""
"""
NOTE: the csv table input obtained from the velocity code contains two different sets of velocity, one with 
      the guassian fit and another with the peak value. To choose which one of these velocities you want to use
      you need to change the value of 'Vcent' either to True or False. 
      Use True, if you would like to use the gaussian fit velocity 
      Use False, if you would like to use the peak velocity. 
"""
#Make sure that you have all your files in the same directory including this python code. This includes the velocity csv
# and the baysian_distance calculator code.  
import csv
import os
csvfile = open('outergalaxy_vel_with_gf.csv', 'r') # csv table with the sources
inpfile = open('sources_info.inp', 'w') # text file used by the calculator
spreadsheet = csv.reader(csvfile, delimiter=',', quotechar='"')
next(spreadsheet, None) # Skip the header
#next(spreadsheet, None) # Skip the header

#writing in the inpfile we opened for writing.
inpfile.write('!  Source          Long    Lat    Vlsr      +/-    P(far)  Extra\n') #heading1
inpfile.write('!                  (deg)   (deg)  (km/s)                        \n') #heading 2

Vcent=True
for row in spreadsheet:
   if Vcent ==True:
       if float(row[2]) >= 0:
           source_name = 'G%06.2f+%05.2f' % (float(row[1]), float(row[2]))
       else :
           source_name = 'G%06.2f%06.2f' % (float(row[1]), float(row[2]))
  
       if row[6] != 'bad fit' and row[6] != 'Data_NA':       
           inpfile.write('%s %8.2f %6.2f %3.1f %3.1f %s %s \n' % (source_name, float(row[1]), float(row[2]),  float(row[6]),  float(row[11]), row[12].center(8), row[13].center(8) ))
   else: 
      if float(row[2]) >= 0:
           source_name = 'G%06.2f+%05.2f' % (float(row[1]), float(row[2]))
      else :
           source_name = 'G%06.2f%06.2f' % (float(row[1]), float(row[2]))
  
      if row[4] != 'NA' and row[4] != 'Data_NA':       
           inpfile.write('%s %8.2f %6.2f %3.1f %3.1f %s %s \n' % (source_name, float(row[1]), float(row[2]),  float(row[4]),  float(row[11]), row[12].center(8), row[13].center(8) ))
 
#Make sure that '+/-' and P(far) values are float type so that the fortran code can compile them. 
inpfile.close()

# Reset the csv file to the beginning 
csvfile.seek(0)

os.system("gfortran-7 Bayesian_distance_2019_fromlist_v2.4.f") #gfortran is a command to make it the f file run and give output in a.out file.
os.system("./a.out > outergalaxy_YB_results.txt") # then we pump the reuslt into a results.txt file which will be created automatically in the same direcotry. 

csvfile.close()

csvfile = open('outergalaxy_vel_with_gf.csv', 'r') # csv table with the sources (reffered below as sources_info_excel.csv, as the name can change depening on your input file name)
spreadsheet = csv.reader(csvfile, delimiter=',', quotechar='"')
next(spreadsheet, None) # Skip the header

results = open('outergalaxy_YB_results.txt', 'r')
csv_results = open('Outergalaxy_YB_Distance_results.csv', 'w')

for i in range(116): #put the number in the range depending on the number of lines the fortran code producees as headers. In this case The calculator produces a header that's 116 lines long so we skip them. 
   next(results, None)
  
#Depending on the result.txt file variables we create columns in our csv file. 
#The following is the heading:
fieldnames = [ 'YB', 'Lon', 'Lat', 'Vlsr','+/-', 'Dist','+/-', 'Integrated Probability','Arm', 'Distance','+/-', 'Integrated Probability', 'Arm'] 
writer = csv.DictWriter(csv_results, fieldnames=fieldnames)
writer.writeheader()

# Go through every row in both results.txt and sources_info_excel.csv at the same time
for row_results, row_csv in zip(results, spreadsheet):
   words = row_results.split()
   if row_csv[4] != 'NA' and row_csv[4] != 'Data_NA' and words[0] != 'No': #skips the lines that don't have data velocity data in them from both the files
       
       writer.writerow({'YB': row_csv[0],
                    'Lon':words[0],
                    'Lat':words[1], 
                    'Vlsr':words[2], 
                    '+/-':words[3],
                    'Dist':words[4],
                    '+/-': words[5],
                    'Integrated Probability':words[6], 
                    'Arm':words[7], 
                    'Distance':words[8], 
                    '+/-': words[9], 
                    'Integrated Probability': words[10],
                    'Arm': words[11]})

results.close() 
csv_results.close()
csvfile.close()