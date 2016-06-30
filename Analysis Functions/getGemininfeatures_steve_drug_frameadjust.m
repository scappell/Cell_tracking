function [risetime,badtraces]=getGemininfeatures_steve_drug_frameadjust(sampletraces,samplestats,quiescentanalysis,Drugoption)

cutoff_value_for_signal=0.015;  %0.025   %value cutoff for the maximum difference between points. 
size_of_window=10; %10 
totalframes=size(sampletraces,2);

[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
risetime=ones(samplesize,1)*NaN;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
% keyboard;
altstore=sigstore;
for i=1:samplesize
    signal_total=sampletraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2); %added this back in;
    if firstframe(i)<totalframes*0.8; %65;
        lastframe(i)=totalframes;
    if quiescentanalysis & firstframe(i)<10;
       new_firstframe=find(signal_total(20:end)<0.03,1,'first');
       if isempty(new_firstframe)
           new_firstframe=-18;
       end
       firstframe(i)=new_firstframe+19;
       
       if lastframe(i)-firstframe(i)<10
           firstframe(i)=lastframe(i)-10;
       end
    end
    signal=signal_total(firstframe(i):lastframe(i));
    signal_smooth=smooth(signal,5)'; %incoming signal always raw (unsmoothened)
    sigstore_smooth(i,firstframe(i):lastframe(i))=signal_smooth;
    sigstore(i,firstframe(i):lastframe(i))=signal;
    %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prebuffer=5;
%     if size(signal_smooth,2)<=prebuffer
% %         prebuffer=size(signal_smooth,2)-1;
%         extension_of_NaNs_for_short_traces=NaN(1,20);
%         signal=[signal extension_of_NaNs_for_short_traces];
%         signal_smooth=[signal_smooth extension_of_NaNs_for_short_traces];
%         lastframe(i)=size(signal_smooth,2)+20;
%     end
    badtraces_1(i)=signal_smooth(prebuffer+1)>0.03;
    badtraces_noise(i)=max(diff(signal(1:end-5),1,2))>0.26;
    badtraces(i)=badtraces_1(i)| badtraces_noise(i);
    
    sig_fwdslope=10*getslope_forward_avg(signal_smooth,1:5);
    altstore(i,firstframe(i):lastframe(i))=sig_fwdslope;
    tempsearch=find(sig_fwdslope>0.05 & abs(signal_smooth)<0.032,1,'last');
    if isempty(tempsearch) || signal_smooth(end)<0.05
        risetime(i)=NaN;
        risetime_temp(i)=NaN;
    else
        risetime_temp(i)=tempsearch;
        risetime(i)=risetime_temp(i)+firstframe(i)-1; %return absolute POI rather than relative to mitosis
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
moving_value=NaN(size(signal,2),size_of_window); 
signal_matrix_first_frame=NaN(size(signal,2),1);    %make empty matrix the same size as moving_value. 

signal_matrix_first_frame=signal';    %make a 3D matrix with each geminin trace in a different matrix
    for c=1:size(signal,2)-size_of_window;
        moving_value(c,:)=signal(c:c+size_of_window-1);    %make a 3D matrix with a moving window of 10 frames for each timepoint. 
    end

signal_matrix_first_frame=repmat(signal_matrix_first_frame,1,size_of_window);    %creates a matrix of the first value of each window
signal_subtract_matrix=moving_value-signal_matrix_first_frame;     %subtract each value in window from the first value
signal_substract_max=max(signal_subtract_matrix,[],2);  %find the maximum value of each window
signal_subtract_matrix2=signal_subtract_matrix;         %re-define the matrix to a new variable
% 
signal_subtract_matrix2(:,1)=0.0001;         %first first value of the subtraction matrix is by definitin Zero. This intereferes with later operations, so this converts the first value to a very low, non-zero number
signal_subtract_matrix2(signal_subtract_matrix2<=0)=0;    %make any negative value, or zero, equal to zero.
signal_subtract_matrix2(isnan(signal_subtract_matrix2))=0;   %make any NaNs equal to zero
signal_subtract_matrix2(signal_subtract_matrix2>0)=1;        %make any positive number equal to 1. This makes the matrix binary
multiply_signal_subtract_matrix2=prod(signal_subtract_matrix2,2);   %multiply every value in every row of the binary matrix. Therefore, if any of the values of the subtraction matrix were negative, this value will be zero. If all values were positive, the value will be 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

ind=find(multiply_signal_subtract_matrix2(11:end,:)==1 & signal_substract_max(11:end,1)>cutoff_value_for_signal,1,'first'); %find time of geminin rise
    if isempty(ind)
        ind=NaN;
    end
    
risetime_steve(i)=ind+10+firstframe(i);

    if risetime(i)-risetime_steve(i)<15;
        risetime_steve_avg(i)=floor(((risetime_steve(i)*0.9)+risetime(i))/2);
    else risetime_steve_avg(i)=risetime_steve(i);
    end
    
if ~quiescentanalysis & ~Drugoption
    if risetime_steve_avg(i)-firstframe(i)<length(signal)-10;
        if signal_total(risetime_steve_avg(i)+10)<0.035;
            risetime_steve_avg(i)=risetime(i)-3;
        end
    end
    
    if risetime_steve_avg(i)-firstframe(i)<length(signal)-20;
        if signal_total(risetime_steve_avg(i)+20)<0.1;
            risetime_steve_avg(i)=risetime(i)-3;
        end
    end
    if risetime_steve_avg(i)-firstframe(i)<length(signal)-30;
        if signal_total(risetime_steve_avg(i)+30)<0.2;
            risetime_steve_avg(i)=risetime(i)-3;
        end
    end
    if risetime_steve_avg(i)-firstframe(i)<length(signal)-40;
        if signal_total(risetime_steve_avg(i)+40)<0.3;
            risetime_steve_avg(i)=risetime(i)-3;
        end
    end
    if risetime_steve_avg(i)-firstframe(i)<length(signal)-50;
        if signal_total(risetime_steve_avg(i)+50)<0.3;
            risetime_steve_avg(i)=risetime(i)-3;
        end
    end
    if risetime_steve_avg(i)-firstframe(i)<length(signal)-60;
        if signal_total(risetime_steve_avg(i)+60)<0.3;
            risetime_steve_avg(i)=risetime(i)-3;
        end
    end
end

if signal(5:end)<0.07;
    risetime_steve_avg(i)=NaN;
end

else
    risetime_steve_avg(i)=NaN;
end
    
end
risetime_new=risetime_steve_avg';
risetime=risetime_steve_avg';
% keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
for i=1:96
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    if ~isnan(frames)
       plot(1:length(frames),sigstore(i,frames));
       ylim([0 1]);
    end
    hold on;
    if ~isnan(risetime(i))
        plot(risetime(i)-firstframe(i),sigstore(i,frames(risetime(i)-firstframe(i))),'go','markerfacecolor','g','markersize',6);
    end
    if badtraces(i)==1 
        continue;
    end
    hold on;
    if ~isnan(frames)
        plot(1:length(frames),altstore(i,frames),'r');
    end
    title(num2str(i));
end
%}