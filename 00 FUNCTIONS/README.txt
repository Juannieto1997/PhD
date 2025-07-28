README
# Last update
2025 07 28
# Functions Folder 

This folder contains all the functions necesary to run the above mentioned codes. 
It contains both public functions and custom functions

#Function Directory

- abfload = function to read abf files from clampfit (obtained from matworks Copyright (c) 2004, Harald Hentschke)
- dxf2binary_mask = Create binary mask from geometry, (by TOM AUBIER)
- f_AppendCells = Append 2 arrays of cells using a for loop (Custom Function)
- f_boxplot = Fast plot of a box plot, used to generate the comparison between different duration of FUS (Custom Function)
	TODO:   Add documentation 
		Generalize for more cases 
		Validate with more data  
-f_cell2matwExclusion = Convert a cell into a matriz with the option to eliminate specific cells of the original cell.  
- f_cut = Cut a single channel of data acording to given parameters. (Custom Function) 
	TODO: 	Modify for efficiency, generate a new vector with all the parameters to cut and reshape the final vector
- f_CutData = Implementation of f_cut, Identify artifacts for cutting the data or uses given parameters to identify points to cut the data. (Custom Function) 
- f_detectPeaks = Fully automatic peak detection, does not detect FUS artifact main use for electrical stimulation. (Custom Function)
	TODO: 	Test for espontaneus activity
		Add documentation for tested parameters 
- f_FeatExtract = Extract the features of a sigal 
	TODO: improve to be consistent with data 
- f_FFTfilt = FFT filter, good for notch and visualization
- f_FindParams = extract the parameters (amplitude, duration, etc..) from the EFP (Custom Function)
	TODO: 	Validate points
- f_FindParamsMat = implementation of f_FindParams for several channels. (Custom Function)
-f_HeatMapMEA12x12 = NOT SURE WHAT DOES IT DOES 
- f_ImageMatrix = Marios version of imagesc (Mario Valderrama)
- f_loadFeatures = ...
	TODO:	Check implentation to understand function 
		Add Documentation
- f_MeanCompleteMatriz = Modified from Davids original function, plots the average of the data in the MEA 120 electrode chip format. (Custom Function)
	TODO: 	Correct to include the raw data current state displays only the average
		Export with transparent background
- f_mkVid = Make the video of the displacement (based on Ivan Suarez Code) 
- f_MorseAWTransformMatlab = Perform Moars transformation in a vector (Mario Valderrama)
- f_organizeChannels = Function to organize the channels in a growing order a cording to the nomenclature in the MEA system 120 grid. (Letter-ChannelNo)
- f_plotBigMat = Plot for the big matrix of the MEA 
- f_Promsg = Display a progress message in loops (Custom function)
- f_rasterplot = VALIDATE
- f_ReadRawData = Fuction to open RAW files from the MCS data tool (David Henao)
		TODO: Add documentation
- f_Searchfiles = Recursive search for all the files with a given extension in a given folder and its subfolders (Custom Function)
- f_TimeFreq = Simplified implementation fo Moarse function above (Custom function)
- f_ValidatePeaks = validation of predetected peaks. 
- f_ValWindow = 
- f_WinDetect = 
- f_ZeroCross = 
- kICA = independent component analisis (Brian Moore)
- pvpmod = evaluate parameter/value pairs (obtained from matworks Copyright (c) 2004, Harald Hentschke)
- violin = Violin plot (Hoffmann H, 2015: violin.m - Simple violin plot using matlab default kernel density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany. hhoffmann@uni-bonn.de)