This code is from Cappell et al., 2016, Cell 166, 167-180. June 30, 2016. 

Running the scripts:
Run Script_1 first. This will track cells and extract the fluoresence itensity for each channel imaged. 
Run Script_2 second. This will merge the data from daughter cells with mother cells. 

File Naming: 
Each image file must be named in the following format: 
Row_Column_Site_Channel_Frame. 
For example, an image taken in a 96 well plate in well B7, in the RFP channel, and at
timepoint number 25 will have the name:  "2_7_1_RFP_25.tif"
If you take multiple images per well you can designate different well sites. (eg "2_7_2_RFP_25.tif")

Any deviations from this naming format will require adjustments to the tracking and analysis codes. 

File Hierarchy: 
Image files from an experiment are stored in 1 parent folder, with subfolders for each Well and site. 
For example, if the experiment includes timelapse images from well B7 but 3 different sites, the folder hierarchy will be: 
ParentFolder>2_7_1 and 2_7_2 and 2_7_3
Each subfolder contains all the images taken from all the timepoints and channels for that well and site. 
Deviations from this hierarchy will require some small adjustments to the tracking code so that the correct image files are found.

Data Structure:
The data from Script_1 is stored in a file named 'Tracedata_well_column_row.mat'
It includes three data matraces: 'tracedata', 'genealogy', and 'jitters'
Tracedata is a 3 dimensional matrix. Rows are each cell tracked, Columns are each Timepoint, the 3rd dimension encodes various parameters extracted. In order they are:  X coordinates of centroid of nucleus, Y coordinates of centroid of nucleus, Nuclear area, Nuclear mass, median nuclear intensity of Channel 1, median nuclear intensity of Channel 2, 75th percentile of ring intensity of Channel 2, Median ring intensity of Channel 2, median nuclear intensity of Channel 3, 75th percentile of ring intensity of Channel 3, Median ring intensity of Channel 3.
Genealogy includes data on when a cell was born and who its mother cell was
Jitters includes the X and Y pixel shift between frames. Correcting for the jitters greatly improves the tracking. Script_1 measures the jitters and corrects for them automatically, but it also stores the data in Jitters in case it could be used later.
Script_2 extracts all the data from Tracedata and Genealogy and combines data from mothers and daughters. It also filters out poorly tracked cells and cells with signals that are too noisy. 

