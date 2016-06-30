clear all;close all;clc;
%%
pth='/Users/Steve/Documents/Meyer_Lab/Data/Fucci Sensors/2016-02-03-MCF10A-MLN-siRNA/Data/'; %Path where data is stored
%%% Experiment Paremeters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conditions={ 
    'siControl',2,1:4,1:2,[0 0 0];  %condition1 Name, Rows, Columns, Sites, Color (eg [R G B])
    'sip21',6,1:4,1:2,[0.8 0 0];    %condition2 Name, Rows, Columns, Sites, Color (eg [R G B])
    'sip27',7,1:4,1:2,[0 0.8 0];    %condition3 Name, Rows, Columns, Sites, Color (eg [R G B])
    'sip21/sip27',8,1:4,1:2,[1 0.5 0];   %condition4 Name, Rows, Columns, Sites, Color (eg [R G B])
};

frames_per_hour=5;  %Frame Rate
frame_drug_added=0;  %Frame treatment was added. Set to Zero if no treatment was added
%%% Analysis Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ring2option=1;       %0:Nuclear Levels; 1:Cytoplasmic to Nuclear Ratio
ring3option=0;       %0:Nuclear Levels; 1:Cytoplasmic to Nuclear Ratio
signal3option=1;     %0: onl1 1 sensor expressed;  1: Two sensors Expressed
IFoption=0;          %0:No IF 1:IF
motheroption=0;      %0:no gating 1:Must have been a mother 2: Was not a Mother
daughteroption=1;    %0:no gating 1:Must have been seen born 2:Was not seen born
quiescentanalysis=0;   %0:cycling cells 1:Cells coming out of serum starvation
POIoption=1;           %0:dont find points of interest 1:find points of interest (eg, G1/S time)
POIcdkoption=0;
POIgemininoption=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sensor]=combine_wells_analysis_v1(conditions,pth,signal3option,ring2option,ring3option,IFoption,motheroption,daughteroption,quiescentanalysis,POIoption,POIcdkoption,POIgemininoption); %combines all the data and traces from various wells into 1 condition structure
numframes=size(sensor(1).signal2,2);  %extracts number of frames in the movie
xtime=1:numframes;  %creates a vector for each frame in the movie. Useful as an X-axis. 
%%
condnum=length(sensor); %Extracts number of conditions
for cond=1:condnum;     %Empty "For" Loop for doing analysis on each condition.  
   sensor(cond).generic_data=sensor(cond).signal2;   %example
end
