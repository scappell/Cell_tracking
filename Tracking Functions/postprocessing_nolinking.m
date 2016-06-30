function [tracedata,genealogy,jitters]=postprocessing_nolinking(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr)
%%% interpolate merged traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxcellid=max(cellid);
for i=1:maxcellid
    [validframes,numsets]=bwlabel(~isnan(tracedata(i,:,1)));
    if numsets>1
        for j=1:(numsets-1)
            mergepreframe=find(validframes==j,1,'last');
            mergepostframe=find(validframes==j+1,1,'first');
            mergeframes=(mergepreframe+1:mergepostframe-1)';
            datainit=squeeze(tracedata(i,mergepreframe,:));
            datadelta=squeeze(tracedata(i,mergepostframe,:)-tracedata(i,mergepreframe,:));
            totalstep=mergepostframe-mergepreframe;
            eachstep=mergeframes-mergepreframe;
            tracedata(i,mergeframes,:)=ones(numel(mergeframes),1)*datainit'+(eachstep*datadelta')/totalstep;
        end
    end
end
%%% interpolate bad frames %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[tracedata,jitters]=interpolateframes(tracedata,jitters,badframes,tracking);
%%% remove excess tracks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracklengthcutoff=5; %default 5
excesstracks=(1:maxcellnum)'>maxcellid;
tracedata(excesstracks,:,:)=[];
tracking(excesstracks,:)=[];
%%% repair broken tracks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numcells=size(tracedata,1);
emptytracks=zeros(numcells,1);
% brokentracks=zeros(numcells,1);
tracestats=ones(numcells,2)*NaN;
for c=1:numcells
    if ~isempty(find(~isnan(tracedata(c,:,1)),1))
        tracestats(c,1)=find(~isnan(tracedata(c,:,1)),1,'first');
        tracestats(c,2)=find(~isnan(tracedata(c,:,1)),1,'last');        
    else
        emptytracks(c)=1;
    end
end
% totalframes=size(tracedata,2);
% fin=find(tracestats(:,2)<totalframes-1);
% nodaughters=zeros(numel(fin),1);
% for i=1:numel(fin)
%     nodaughters(i)=isempty(find(tracking(:,1)==fin(i),1));
% end
% fin=fin(logical(nodaughters));
% finend=tracestats(fin,2);
% [finend,tidx]=sort(finend,'descend');
% fin=fin(tidx);
% orphan=find(~isnan(tracestats(:,1)) & isnan(tracking(:,1)) & isnan(tracking(:,4)));
% ostart=tracestats(orphan,1);
% numorph=numel(orphan);
% omass=ones(numorph,1)*NaN; ox=omass; oy=omass;
% for i=1:numel(orphan)
%     omass(i)=tracedata(orphan(i),ostart(i),4); ox(i)=tracedata(orphan(i),ostart(i),1); oy(i)=tracedata(orphan(i),ostart(i),2);
% end
% winrad=4*nucr;
% for i=1:numel(fin)
%     tmass=tracedata(fin(i),finend(i),4); tx=tracedata(fin(i),finend(i),1); ty=tracedata(fin(i),finend(i),2);
%     c=find(ostart-finend(i)<4 & ostart-finend(i)>0 & abs(ox-tx)<winrad & abs(oy-ty)<winrad & abs((omass-tmass)/tmass)<0.3);
%     if numel(c)>0
%         c=c(find(ostart(c)==min(ostart(c)),1)); %multiple hits
%         concatframes=ostart(c):totalframes;
%         tracedata(fin(i),concatframes,:)=tracedata(orphan(c),concatframes,:);
%         %tracedata(term(i),tend(i)+1,:)=(tracedata(term(i),tend(i),:)+tracedata(term(i),ostart(c),:))/2;
%         framediff=ostart(c)-finend(i);
%         if framediff>1
%             datadiff=tracedata(fin(i),ostart(c),:) - tracedata(fin(i),finend(i),:);
%             datastep=datadiff/framediff;
%             for j=1:framediff-1
%                 tracedata(fin(i),finend(i)+j,:)=tracedata(fin(i),finend(i),:)+datastep*j;
%             end
%         end
%         tracking((tracking(:,1)==orphan(c)),1)=fin(i); %reassign daughters
%         ostart(c)=NaN; omass(c)=NaN; ox(c)=NaN; oy(c)=NaN;
%         brokentracks(orphan(c))=1;
%     end
% end
%%% remove bad tracks (empty, broken, short, merged) %%%%%%%%%%%%%%%%%%%%%%
tracklength=ones(maxcellid,1)*NaN;
for c=1:maxcellid
    lastframe=find(~isnan(tracedata(c,:,1)),1,'last');
    firstframe=find(~isnan(tracedata(c,:,1)),1,'first');
    if isempty(lastframe)
        tracklength(c)=0;
    else
        tracklength(c)=lastframe-firstframe;
    end
end
mothers=unique(tracking(~isnan(tracking(:,1)),1));
shorttracks=isnan(tracking(:,1)) & ~ismember((1:maxcellid)',mothers) & tracklength<tracklengthcutoff;
mergetracks=~isnan(tracking(:,4));
removetracks=find(emptytracks | shorttracks | mergetracks);
tracedata(removetracks,:,:)=[];
tracking(removetracks,:)=[];
%%% update genealogy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genealogy=tracking(:,1);
for m=mothers'
    numremovedbefore=sum(removetracks<m);
    daughtersidx=genealogy==m;
    newm=m-numremovedbefore;
    genealogy(daughtersidx)=newm;
end