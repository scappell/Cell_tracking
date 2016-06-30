function [risetime,badtraces]=getCdk2features_steve_drug(sampletraces,samplestats,minlength,Drugoption,quiescentanalysis)

cutoff_value_for_signal=0.085;  %0.025   %value cutoff for the maximum difference between points. 
size_of_window=10;  
totalframes=size(sampletraces,2);

slopewindow=min([minlength 10]); %default 40
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
%minval=ones(samplesize,1)*NaN;
risetime=ones(samplesize,1)*NaN;
risetimedb=risetime;
%riseslope=ones(samplesize,1)*NaN;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore1=ones(samplesize,tracelength)*NaN;
altvar1=ones(samplesize,1)*NaN;
% keyboard;
for i=1:samplesize
    signal_total=sampletraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    if firstframe(i)<totalframes-10 & lastframe(i)-firstframe(i)>10;
    %signal=signal(firstframe(i):lastframe(i));
    signal_total_smooth(i,:)=smooth(signal_total,5);
    signal=signal_total(firstframe(i):lastframe(i));
    signal_smooth=smooth(signal)'; %incoming signal always raw (unsmoothened)
    numframes=lastframe(i)-firstframe(i)+1;
    sigstore(i,firstframe(i):lastframe(i))=signal_smooth;
    %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     badtraces(i)=min(signal(1:slopewindow))>1; %will remove permanently high signals
    
    badtraces(i)=min(signal_total(firstframe(i):firstframe(i)+slopewindow))>1; %will remove permanently high signals
    
    
    %%%%%% general requirements
    gate=zeros(1,numframes);
    for j=2:numframes-slopewindow
        pastheight=max(signal_smooth(1:j))<1.5;
        futureheight=max(signal_smooth(j:end))>1;
        gate(j)=pastheight && futureheight;
    end
    %%%%%% filter for inflection points
    sig_fwdslope_short=getslope_forward_avg(signal_smooth,6:10);
    sig_fwdslope_long=getslope_forward_avg(signal_smooth,round(slopewindow/2+1):slopewindow);
    sig_time=(1:length(signal_smooth))/length(signal_smooth);
    filter=sig_fwdslope_short+sig_fwdslope_long-4*signal_smooth+sig_time+10;
    %%%%%% combine gate and filter
    filter=filter.*gate;
    %%% find cdk2start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filtermax=max(filter);
    if isnan(filtermax)
        filtermax=0;
    end
    if sum(~isnan(filter))>0
        risetime(i)=find(filter==filtermax,1,'first');
    else
        risetime(i)=NaN;
    end
    if filtermax<=0
        risetime(i)=NaN;
%         continue;
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
moving_value=NaN(size(signal_smooth,2),size_of_window); 
signal_matrix_first_frame=NaN(size(signal_smooth,2),1);    %make empty matrix the same size as moving_value. 

signal_matrix_first_frame=signal_smooth';    %make a 3D matrix with each geminin trace in a different matrix
    for c=1:size(signal_smooth,2)-size_of_window;
        moving_value(c,:)=signal_smooth(c:c+size_of_window-1);    %make a 3D matrix with a moving window of 10 frames for each timepoint. 
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

signal_substract_first_and_last=moving_value(:,end)-moving_value(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
risetime2_temp=find(multiply_signal_subtract_matrix2(1:end,:)==1 & signal_substract_max(1:end,1)>cutoff_value_for_signal & signal_substract_first_and_last(1:end,1)>0.01,1,'first'); %find time of geminin rise
if isempty(risetime2_temp)
    risetime2_temp=NaN;
end    
    risetime2db(i)=risetime2_temp;
    risetime2(i)=risetime2_temp+firstframe(i)-1; %return absolute POI rather than relative to mitosis
if ~Drugoption
    if ~isnan(risetime2(i))
        if risetime2(i)<=lastframe(i)-20;
            if signal_total_smooth(i,risetime2(i)+20)-signal_total_smooth(i,risetime2(i))<0.1 | signal_total_smooth(i,risetime2(i)+20)<0.75
                risetime2(i)=NaN;
            end
        end
    end
end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    %%% calc minval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     minval(i)=min(signal_smooth(1:risetime(i)));
    %%% collect additional info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     altvar1(i)=sig_fwdslope_long(risetime(i));
    risetimedb(i)=risetime(i);
    risetime(i)=risetime(i)+firstframe(i); %return absolute POI rather than relative to mitosis
    altstore1(i,firstframe(i):lastframe(i))=filter;
    
%%% Select which option works best %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
if~Drugoption
if ~isnan(risetime(i))
    if risetime(i)<=lastframe(i)-20;
        if signal_total_smooth(i,risetime(i)+20)-signal_total_smooth(i,risetime(i))<0.1 | signal_total_smooth(i,risetime(i)+20)<0.6
            risetime(i)=risetime2(i);
        end
    end
end
end

if isnan(risetime(i)) 
    risetime(i)=risetime2(i);
end

if ~isnan(risetime(i)) 
    if signal_total_smooth(i,risetime(i)+1)-signal_total_smooth(i,risetime(i))<0
        risetime(i)=risetime2(i);
    end
end

if ~isnan(risetime(i)) & ~isnan(risetime2(i))
    if signal_total_smooth(i,risetime2(i)) < signal_total_smooth(i,risetime(i));
        risetime(i)=risetime2(i);
    end
end
if quiescentanalysis
    if ~isnan(risetime(i)) & risetime(i)<lastframe(i)-10;
        min_signal_trace=min(signal_total_smooth(i,risetime(i)+20:lastframe(i)-10));
        if signal_total_smooth(i,risetime(i))>min_signal_trace;
            risetime(i)=NaN;
        end
    end
end

if isnan(risetime(i)) & signal_total_smooth(i,lastframe(i))>1.1
    badtraces(i)=1;
end

    
    else
       risetime(i)=NaN; 
    end

end

% keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
altstore1(altstore1==0)=NaN;
altstore_norm=normalizeMyTracesGeminin_alt2(altstore1,0.01);
for i=1:samplesize
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    plot(1:length(frames),sigstore(i,frames));
    axis([1 length(frames) 0.3 2.5]);

    if badtraces(i)==1 || isnan(risetimedb(i))
        continue;
    end
    hold on;
    plot(risetimedb(i),sigstore(i,frames(risetimedb(i))),'go','markerfacecolor','g','markersize',6);
    plot(1:length(frames),altstore_norm(i,frames)+0.5,'r');
    plot(risetimedb(i),altstore_norm(i,frames(risetimedb(i)))+0.5,'go','markerfacecolor','g','markersize',6);
end
%}

%{
altstore1(altstore1==0)=NaN;
altstore_norm=normalizeMyTracesGeminin_alt2(altstore1,0.01);
for i=1:96;%samplesize
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    numframes=size(altstore1,2);    
    frames=1:numframes;
    plot(1:length(frames),signal_total_smooth(i,frames));
    axis([1 length(frames) 0.3 2.5]);title(num2str(i));
    hold on;
    if ~isnan(risetime(i))
        plot(risetime(i),signal_total_smooth(i,frames(risetime(i))),'go','markerfacecolor','g','markersize',6);
    end
    if ~isnan(risetime2(i))
        plot(risetime2(i),signal_total_smooth(i,frames(risetime2(i))),'ko','markerfacecolor','k','markersize',6);
    else  plot(lastframe(i),signal_total_smooth(i,frames(lastframe(i))),'ko','markerfacecolor','k','markersize',15);
    end

    if badtraces(i)==1
        continue;
    end
    hold on;
    plot(1:length(frames),altstore_norm(i,frames)+0.5,'r');
end
%}






