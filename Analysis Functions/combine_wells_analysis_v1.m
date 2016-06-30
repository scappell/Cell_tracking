function [sensor]=combine_wells_mingyu_analysis_rev10(conditions,datadir,signal3option,ring2option,ring3option,IFoption,motheroption,daughteroption,quiescentanalysis,POIoption,POIsignal2option,POIsignal3option)
Drugoption=0;
% motheroption=0; %0:no gating 1:mothers 2:no mothers
% daughteroption=0; %0:no gating 1:daughters 2:no daughters
if quiescentanalysis
    motheroption=2; daughteroption=2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allnames=conditions(:,1);
[~,uidx]=unique(allnames,'first');
uniquenames=allnames(sort(uidx));
uniquecondnum=numel(uniquenames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condnum=size(conditions,1);
for ii=1:uniquecondnum
    condrow=find(ismember(conditions(:,1),uniquenames{ii}));
    tracedata=[];
    tracestats=[];
    motherstats=[];
    IFdata=[];
    wellindex=[];
    cellID=[];
    cc=0;
    for c=condrow'
        rowmat=cell2mat(conditions(c,2));
        colmat=cell2mat(conditions(c,3));
        sitemat=cell2mat(conditions(c,4));
        for row=rowmat
            for col=colmat
                for site=sitemat
                    cc=cc+1;
                    shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                    if exist([datadir,'tracedata_',shot,'.mat']);
                        if IFoption & exist([datadir,'IF_',shot,'.mat']);
                            [tracedatatemp,tracestatstemp,motherstatstemp,IFdatatemp,jitters,samplecellsID]=gathertracedata_1_rev05(datadir,shot,motheroption,daughteroption,IFoption);
                            tracedata=[tracedata;tracedatatemp];
                            tracestats=[tracestats;tracestatstemp];
                            motherstats=[motherstats;motherstatstemp];
                            IFdata=[IFdata;IFdatatemp];
                            cellID=[cellID;samplecellsID];
                            wellindextemp=ones(size(tracedatatemp,1),3);
                            wellindextemp(:,1)=wellindextemp(:,1)*row;wellindextemp(:,2)=wellindextemp(:,2)*col;wellindextemp(:,3)=wellindextemp(:,3)*site;
                            wellindex=[wellindex;wellindextemp];
                            
                        elseif ~IFoption;
                            [tracedatatemp,tracestatstemp,motherstatstemp,IFdatatemp,jitters,samplecellsID]=gathertracedata_1_rev05(datadir,shot,motheroption,daughteroption,IFoption);
                            tracedata=[tracedata;tracedatatemp];
                            tracestats=[tracestats;tracestatstemp];
                            motherstats=[motherstats;motherstatstemp];
                            IFdata=[IFdata;IFdatatemp];
                            cellID=[cellID;samplecellsID];
                            wellindextemp=ones(size(tracedatatemp,1),3);
                            wellindextemp(:,1)=wellindextemp(:,1)*row;wellindextemp(:,2)=wellindextemp(:,2)*col;wellindextemp(:,3)=wellindextemp(:,3)*site;
                            wellindex=[wellindex;wellindextemp];
                        end
                    end
                end
            end
        end
    end
    numframes=size(tracedata,2);
    minlengthtrace=numframes-ceil((numframes*0.4));
    if quiescentanalysis
        minlengthtrace=numframes-10;
    end
    minlengthmother=5;
%%% gate on length %%%%%%%%
badlengths=tracestats(:,2)-motherstats(:,1)<minlengthtrace | tracestats(:,3)<1;
if ~quiescentanalysis
    badmitosis=tracestats(:,1)<=5;
else
    badmitosis=badlengths;
end
badcells=badlengths | badmitosis;

tracedata=tracedata(~badcells,:,:);
tracestats=tracestats(~badcells,:);
motherstats=motherstats(~badcells,:);
IFdata=IFdata(~badcells,:);
cellID=cellID(~badcells,:);
wellindex=wellindex(~badcells,:);
%%% smooth and gate signal 2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ring2option
    nucchannel_2=6; cytochannel_2=8;
    maxthresh=5; %threshold above which max of each trace must be  %150
    noisethresh=0.25; %threshold rate of DHBratio change (not absolute) above which trace is considered noisy
    [traces1,traces1_cyto,badtraces1]=gate_Cdk2_1_rev05(tracedata,nucchannel_2,cytochannel_2,tracestats,minlengthtrace,maxthresh,noisethresh,quiescentanalysis);
else
    %%% smooth and gate nuclear signal 2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nucchannel_2=6; %6
    maxthresh=100;  %threshold above which max of each trace must be   %500
    minthresh=100; %threshold below which min of each trace must be     %20
    if quiescentanalysis
        maxthresh=1;  %threshold above which max of each trace must be   %500
        minthresh=500; %threshold below which min of each trace must be     %20
    end
    maxnoisethresh=2000; %threshold below which max derivative must be
    minnoisethresh=-10000; %threshold above which min derivative must be
    [traces1,badtraces1]=gate_lengthandrange_noise_rev05(tracedata,tracestats,nucchannel_2,minlengthtrace,maxthresh,minthresh,maxnoisethresh,minnoisethresh);
    traces1_cyto=traces1;
end

%%% smooth and gate signal 3 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ring3option & signal3option
    nucchannel_3=9; cytochannel_3=11;
    maxthresh=5; %threshold above which max of each trace must be  %150
    noisethresh=0.25; %threshold rate of DHBratio change (not absolute) above which trace is considered noisy
    [traces2,traces2_cyto,badtraces2]=gate_Cdk2_1_rev05(tracedata,nucchannel_3,cytochannel_3,tracestats,minlengthtrace,maxthresh,noisethresh,quiescentanalysis);
elseif signal3option
    %%% smooth and gate nuclear signal 3 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nucchannel_3=9; %9
    maxthresh=100;  %threshold above which max of each trace must be  %50  %500    %500
    minthresh=100; %threshold below which min of each trace must be  %1000   %60   %200
    if quiescentanalysis
        maxthresh=1;  %threshold above which max of each trace must be   %500
        minthresh=50; %threshold below which min of each trace must be     %20
    end
    maxnoisethresh=2000; %threshold below which max derivative must be %200
    minnoisethresh=-10000; %threshold above which min derivative must be %-400
    [traces2,badtraces2]=gate_lengthandrange_noise_rev05(tracedata,tracestats,nucchannel_3,minlengthtrace,maxthresh,minthresh,maxnoisethresh,minnoisethresh);
    traces2_cyto=traces2;
elseif ~signal3option
    traces2=traces1;
    traces2_cyto=traces1;
    badtraces2=badtraces1;
end

%%% Extract Nuclear Area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_area_channel=3;
nuclear_area=tracedata(:,:,nuc_area_channel);

%%% transform and gate IF data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IFoption
    channelIF=9;
    IFvals=IFdata(:,channelIF); 
    badtracesIF=IFvals<1 | isnan(traces2(:,end)); 
    IFvals=log2(IFvals);
else
    IFvals=ones(size(badtraces1))*NaN;
    badtracesIF=badtraces1;
end
%%% combine gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badtraces=badtraces1 |badtraces2 | badtracesIF; %badtraces

%%% remove gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces1=traces1(~badtraces,:);
traces1_cyto=traces1_cyto(~badtraces,:);
traces2=traces2(~badtraces,:);
traces2_cyto=traces2_cyto(~badtraces,:);
nuclear_area=nuclear_area(~badtraces,:);
tracestats=tracestats(~badtraces,:);
motherstats=motherstats(~badtraces,:);
IFvals=IFvals(~badtraces);
wellindex=wellindex(~badtraces,:);
cellID=cellID(~badtraces);

%%% Quiescence Gate for Geminin data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if quiescentanalysis
%     if ring3option
%         [badtraces_nuclear_sig3_quiescence]=getbadtraces_cdk2_quiescence(traces1);
%     else
%         [badtraces_nuclear_sig3_quiescence]=getbadtraces_quiescence(traces2);
%     end
%     if ring2option
%         [badtraces_ring_sig2_quiescence]=getbadtraces_cdk2_quiescence(traces1);
%     else badtraces_ring_sig2_quiescence=badtraces_nuclear_sig3_quiescence;
%     end
%     badtraces_quiescence=badtraces_nuclear_sig3_quiescence | badtraces_ring_sig2_quiescence;
%     traces1=traces1(~badtraces_quiescence,:);
%     traces2=traces2(~badtraces_quiescence,:);
%     nuclear_area=nuclear_area(~badtraces_quiescence,:);
%     tracestats=tracestats(~badtraces_quiescence,:);
%     motherstats=motherstats(~badtraces_quiescence,:);
%     IFvals=IFvals(~badtraces_quiescence);
%     wellindex=wellindex(~badtraces_quiescence,:);
%     cellID=cellID(~badtraces_quiescence);
% end
%%% gate POI calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if POIoption;
    if quiescentanalysis & ~ring3option
            traces2_norm=normalizeMyTracesGeminin_serum_addition_ming_rev05(traces2,0.01);
    elseif ~quiescentanalysis & ~ring3option
        traces2_norm=normalizetraces(traces2);
    end
    
    if ring3option & POIsignal3option
        [POI_signal3,badtraces_signal3]=getCdk2features_steve_drug(traces2,tracestats,minlengthtrace,Drugoption,quiescentanalysis);
    elseif POIsignal3option
        %[POI_signal3,badtraces_signal3]=getGemininfeatures_steve_drug(traces2_norm,tracestats,quiescentanalysis,Drugoption(ii));
        [POI_signal3,badtraces_signal3]=getGemininfeatures_steve_drug_frameadjust(traces2_norm,tracestats,quiescentanalysis,Drugoption);
    end
    
    if ring2option & POIsignal2option & POIsignal3option & ~ring3option
        [POI_signal2,badtraces_signal2]=getCdk2features_steve_drug(traces1,tracestats,minlengthtrace,Drugoption,quiescentanalysis);
        POI_G1_length=POI_signal3-POI_signal2;
        badtraces_POI_logic=POI_G1_length<0;
        badtraces_POI_logic2=isnan(POI_signal2) & ~isnan(POI_signal3);
    elseif ring2option & POIsignal2option & POIsignal3option & ring3option
        [POI_signal2,badtraces_signal2]=getCdk2features_steve_drug(traces1,tracestats,minlengthtrace,Drugoption,quiescentanalysis);
        badtraces_POI_logic=badtraces_signal2;
        badtraces_POI_logic2=isnan(POI_signal2) & ~isnan(POI_signal3);
    elseif ring2option & POIsignal2option & ~POIsignal3option
        [POI_signal2,badtraces_signal2]=getCdk2features_steve_drug(traces1,tracestats,minlengthtrace,Drugoption,quiescentanalysis);
        badtraces_POI_logic=badtraces_signal2;
        badtraces_POI_logic2=badtraces_signal2;
        badtraces_signal3=badtraces_signal2;
        POI_signal3=NaN(size(POI_signal2));
    elseif ~ring2option | (POIsignal3option & ~POIsignal2option)
        badtraces_signal2=badtraces_signal3;
        badtraces_POI_logic=badtraces_signal3;
        badtraces_POI_logic2=badtraces_signal3;
        POI_signal2=NaN(size(POI_signal3));
    end
    
   badtracesallPOI=badtraces_signal2 | badtraces_signal3 | badtraces_POI_logic | badtraces_POI_logic2;
   traces1=traces1(~badtracesallPOI,:);
   traces1_cyto=traces1_cyto(~badtracesallPOI,:);
   traces2=traces2(~badtracesallPOI,:);
   traces2_cyto=traces2_cyto(~badtracesallPOI,:);
   nuclear_area=nuclear_area(~badtracesallPOI,:);
   tracestats=tracestats(~badtracesallPOI,:);
   motherstats=motherstats(~badtracesallPOI,:);
   IFvals=IFvals(~badtracesallPOI);
   POI_signal3=POI_signal3(~badtracesallPOI);
   if ring2option
       POI_signal2=POI_signal2(~badtracesallPOI);
   else POI_signal2=NaN(length(POI_signal3),1);
   end
   POI_matrix=[tracestats(:,1),POI_signal2,POI_signal3];
else POI_matrix=tracestats(:,1);
end
sensor(ii).signal2=traces1;
sensor(ii).signal2_cyto=traces1_cyto;
sensor(ii).signal3=traces2;
sensor(ii).signal3_cyto=traces2_cyto;
sensor(ii).tracestats=tracestats;
sensor(ii).motherstats=motherstats;
sensor(ii).IF=IFvals;
sensor(ii).nucleus_area=nuclear_area;
sensor(ii).POI=POI_matrix;
sensor(ii).wellindex=wellindex;
sensor(ii).cellID=cellID;
end
end