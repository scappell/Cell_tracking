clear all;close all;clc;
%%% Designate wells to analyze %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rows = 1:4;     
cols = 1:12;   
sites= 1:4;     
%%% Manual Well Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
manualcontrol = 1; %1 or 0
manualwells = [
    2 1 1;
];
%%% DeBug Mode Setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debug_mode=0; %1 or 0   %%This will run the first well and stop half-way through the analysis so you can plot images, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Paths
settings.experiment_name='2016-06-07-MCF10A-CDK12i';  %Make up a new folder name that the data for this experiment will get stored
settings.projectpath=['C:\Users\MeyerLab\Cell_tracking\SampleData\',settings.experiment_name,'\Data\']; %Path where this folder will get made
settings.imagepath='C:\Users\MeyerLab\Cell_tracking\Sample_images\';  % Path where the images are stored
settings.shadingpath='C:\Users\MeyerLab\Cell_tracking\shadingpath\';  %path were a background noise image is located %%Optional
settings.biaspath=[settings.imagepath,'Bias\'];  %Optional path where a shading correction image is located
%%% Options
settings.signal3_option=1; %1=3rd fluorescent protein present; 0=Only 2 fluorescent proteins present
settings.separatedirectories_option=1; %1=Folder for each well,col,site; 0=all files in 1 folder
settings.maskwrite_option=0; %1=save an image with the nucmask; 0=dont save an image with the nucmask
settings.ringcalc_sig2=1; %1=calcuclate ring pixels for signal 2; 0=only calculate nuc signal for signal 2
settings.ringcalc_sig3=0; %1=calcuclate ring pixels for signal 3; 0=only calculate nuc signal for signal 3
settings.bias_sig1_option=0; %1=correct for bleedthrough; 0=dont correct for bleedthrough;
settings.bias_sig2_option=0; %1=correct for bleedthrough; 0=dont correct for bleedthrough;
settings.bias_sig3_option=0; %1=correct for bleedthrough; 0=dont correct for bleedthrough;
settings.shadingcorrection=0; %1=correct for shading; 0=dont correct for shading;
%%% Experiment parameters
settings.StartFrame=1;   %Frame to start analyzing from
settings.EndFrame=10;  %Frame to stop analyzing  %must be greater than StartFrame
settings.magnification=10; %10=10x or 20=20x
settings.binsize=1; %1=bin 1 or 2=bin 2
settings.nucleus_name='CFP_'; %File prefix for nuclear marker
settings.signal2='YFP_'; %File prefix for signal 2
settings.signal3='RFP_'; %File prefix for signal 3
settings.firstsegmethod='single'; %'log' or 'single' or 'double'
settings.secondsegmethod='apriori'; %'log' or 'double' or 'apriori'
%%% Segmentation parameters
settings.nucr=12; %10x bin1=12 and 20x bin2=12
settings.debrisarea=100; %parameter for throwing out small objects
settings.boulderarea=1500; %parameter for throwing out large objects
settings.blobthreshold=-0.03; %parameter for throwing out large objects
%%% Tracking parameters
settings.maxjump=settings.nucr*4; %parameter for throwing out objects that move too fast (due to bad tracking)
settings.masschangethreshold=0.60; %parameter for throwing out objects that have a large change in mass from frame to frame
settings.areachangethreshold=0.60; %parameter for throwing out objects that have a large change in area from frame to frame
settings.daughtervariance=0.10;
settings.blurradius=3; 
settings.blocksize=10000;
settings.eccentricitythresh=0.85;  %threshold for separating two nuclei in contact with each other
settings.compression=4;
settings.ringthresh=100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DeBug Mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if debug_mode==1 & manualcontrol==0;  %enter Debug mode to have the program stop running in the middle of the script
    Timelapse_Tracking_v1_PC(rows(1),cols(1),sites(1),settings,debug_mode)
elseif debug_mode==1 & manualcontrol==1;
    Timelapse_Tracking_v1_PC(manualwells(1,1),manualwells(1,2),manualwells(1,3),settings,debug_mode)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows=length(rows);
numcols=length(cols); 
numsites=length(sites);
shots=numrows*numcols*numsites;
if manualcontrol==1
    shots=size(manualwells,1);
end
time1=tic;
parfor shot=1:shots
    if manualcontrol==1  
        row=manualwells(shot,1); %creates single vector for the rows
        col=manualwells(shot,2); %creates single vector for the columns
        site=manualwells(shot,3); %creates single vector for the sites
    else    %generates a random list of the well, cols, sites so they dont repeat in the parrallel processing
        siteidx=mod(shot,numsites);
        if siteidx==0
            siteidx=numsites;
        end
        site=sites(siteidx);
        colidx=mod(ceil(shot/numsites),numcols);
        if colidx==0
            colidx=numcols;
        end
        col=cols(colidx);
        rowidx=ceil(shot/(numcols*numsites));
        row=rows(rowidx);
    end
    fprintf([num2str(row),'_',num2str(col),'_',num2str(site),'\n']); %display each well being analyzed
    try 
    %%% Timelapse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Timelapse_Tracking_v1_PC(row,col,site,settings,debug_mode)  %This is the engine of the analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    catch
        disp(['Error: ',num2str(row),'_',num2str(col),'_',num2str(site)]);   %Display an error and move to the next well for analysis
    end

end
save([settings.projectpath,settings.experiment_name,'_settings.mat'],'settings'); %saves all the settings to the Data folder. will overwrite with the most recent
toc(time1)
