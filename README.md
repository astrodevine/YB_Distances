# YB_Distances
Codes for calculating distances to YBs using the Reid et al. 2019 Bayesian distance calculator. Below readme notes are from Bezawit Kasaye, June 2022:


# GENERAL READ ME 

The first thing done is calculating velocities. After getting the output velocities, move on to flagging for bad and good fit. Following that to refitting the bad fits and finally to the Distance Calculator. For the outergalaxy refitting isn’t important at the moment as it doesn’t have any proper gaussian fit regardless due to the way the spectrum is. 


# How to use the velocity code:
To use the velocity code for the inner galaxy, you will first need to download the cubes you’ll be needing from the SEDIGISM database and the BU GRS database. 
Link:  http://grunt.bu.edu/grs-stitch/download-all.php 
https://sedigism.mpifr-bonn.mpg.de/cgi-bin-seg/SEDIGISM_DATABASE.cgi 
For the outer galaxy, the cubes are already present in the folder for you. 
After downloading the cubes, create a folder for your velocity code and a sub directory named ‘data_cubes’, where you’ll put your data cubes in. Rename all your data cubes downloaded from the BU-GRS dataset with the format 13COdatacubes_015-017.fits, 13COdatacubes_017-019 (the numbers are the degrees each datacube covers, so it should be different for all the codes). Leave the ones downloaded from SEDIDISM as they are. 
You should also create other subdirectories in the velocity directory named the following: ‘66percent_YB_graphs’, ‘amplitude is bounded and sd is 2(new)’, ‘residuals when sd is 2’ where the code will store the different plots in, and lastly the catalogue input csv file.
After making sure all these are in one folder you can go ahead and run the code. 
I have made sure to add the graphs I selected to calculate sigma value (SD) under folders BU_selected_graphs_for_sd and SEDGISM_selected_graphs_for_sd. You can either use the same ones or choose another one. But if you decide to change them, make sure to change the hardcoded part of the SD function calling section in the code, which tells you what to do with the comments on top of it. 
 
# How to use the Distance calculator code: 
To use this code that incorporates the Bayesian distance calculator you need to be using a machine that can run Fortran code (Linux is preferable). 
You will want to have a folder that contains the distance code, the Baysian_distance_2019_fromlist_v2.4.f , and the output csv from the velocity code. 
The velocity code outputs two different velocities; hence you need to choose which velocity you would like to use to calculate the distance with. To do that, you can change the value of the variable ‘Vcel’ to True or False where if True, it will use the gaussian fit velocity, and if false it will use the peak velocity. 
The results of this code are saved in a csv file,(distance_results.csv) which you can change the name of with what you want to name it. 

# How to use the Guassian Fit checking code: 
name of code: 
Inner Galaxy:  GF_y_N(final version).py
Outer Galaxy: Y_N_outergalaxyy.py
To use this code you need to have the code in the same directory as the velocity code, that way it will have access to the cube files and the csv output file the velocity code generates as well. On top of that, you will need to create a folder in this directory and name it “TRIAL_amplitude is bounded and sd is 2”
IMPORTANT: You will only need this directory if you want to save the plots that you have already worked through. Remember, these plots are the replica of the ones that are created by the velocity code, hence saving them is not really necessary. If you do choose to not save them make sure to go to the two lines that have “ plt.savefig(TRIAL_amplitude is bounded and sd is 2/box_spec_%03i)” commented out, so you won’t run into an error. 

When answering ‘y/n’ what you will want to check is if the orange plot is fitting properly over the highest peak of the spectrum that is plotted in black below the orange. Below are some examples of bad fit and good fits. 
					Bad fits 


					Good fits



# How to use the Gaussian fit replotting clicking code:
Name of code: click_gf(final version)
This code also needs to be in the same directory as the velocity code so it will have access to the data cubes as well as the output csv file from the Guassian fit checking code that has information of which YBs were bad fit and need to be manually refitted.
In this directory you will need to create a new folder and name it click_gf_trial(new)
This is where the new plots you make will be saved. 

When working with this code, if it’s the first time you’re using it, it will automatically open an interactive plot where you get to click on the figure five main points from where the program will extract informations. follow these steps when clicking to make sure you have all the right informations clicked on. These informations will also be displayed on the interactive plot. 
(*Left)click to select   
*Right click to undo selection  
*when done click enter 
 the peak point
 the min value for the width
 the max value for the width     
 the max value for the amplitude 
 the min value for the amplitude

After clicking on these five points and clicking enter a plot with the new fit made based on the points you clicked will appear. Click enter again and you will be prompted to say if this was a good fit or not and you will enter ‘y/n/p’. Refere to the figure examples on how to use the gaussian fit checking code to see which ones are good fit and which ones aren’t. If you feel like the fit cannot be improved or salvaged anymore you can simply enter p and move on to the next YB that had a bad fit. 

