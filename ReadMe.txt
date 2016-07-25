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

