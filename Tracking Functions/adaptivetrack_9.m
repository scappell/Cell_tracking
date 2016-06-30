function [tracedata,curdatatracked,tracking,nuc_label]=adaptivetrack_9(tf,lastgoodframe,curframe,tracedata,curdata,tracking,nuc_raw,nuc_label,jitter,trackparams,debugpackage)
nucr=trackparams{1};
maxjump=trackparams{2};
debrisarea=trackparams{3};
masschangethreshold=trackparams{4};
areachangethreshold=trackparams{5};
daughtervariance=trackparams{6};
daughtermin=-0.5-daughtervariance;
daughtermax=-0.5+daughtervariance;
%nuc_mask=bwmorph(nuc_label,'remove');
%nuc_label_mask=nuc_label.*nuc_mask;
absjitx=jitter(1); absjity=jitter(2);
%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%winrad=3*nucr; %prev:4*nucr
maxcellnum=size(tracedata,1);
numprevtotal=find(~isnan(tracedata(:,lastgoodframe,1)),1,'last');
xprev=tracedata(1:numprevtotal,lastgoodframe,1);
yprev=tracedata(1:numprevtotal,lastgoodframe,2);
areaprev=tracedata(1:numprevtotal,lastgoodframe,3);
massprev=tracedata(1:numprevtotal,lastgoodframe,4);
previd=find(~isnan(xprev));
numprevextant=numel(previd);
xcur=curdata(:,1); ycur=curdata(:,2); areacur=curdata(:,3); masscur=curdata(:,4);
rxcur=round(xcur); rycur=round(ycur);
numcurorg=numel(xcur);
ongoingmerge=~isnan(tracking(:,4)) & isnan(tracking(:,5));
%mergeduration=(curframe-tracking(:,4))>5;
shortmerge=(curframe-tracking(:,4))<=5;

%%% attempt to split current or previous merges %%%%%%%%%%%%%%%%%%%%%%%%%%%
bordermask=zeros(size(nuc_label));
borderflag=zeros(numcurorg,1);
nuc_mask=nuc_label>0;
%[B,denovo_label]=bwboundaries(nuc_mask,'noholes');
B=bwboundaries(nuc_mask,'noholes');
for i=1:numprevextant
    p=previd(i);
%     if abs(xprev(p)-1732)<maxjump && abs(yprev(p)-471)<maxjump && tf==2
%         keyboard;
%     end

    neighbors=find(abs(xcur-xprev(p))<maxjump & abs(ycur-yprev(p))<maxjump);
    if numel(neighbors)==0
        continue;
    end
    [~,cidx]=min(sqrt((xcur(neighbors)-xprev(p)).^2+(ycur(neighbors)-yprev(p)).^2));
    match=neighbors(cidx);
    massdiff=(masscur(match)-massprev(p))/massprev(p);
    if massdiff>masschangethreshold || ongoingmerge(p)
        %[bordermask,borderflag(match)]=attemptsplit_1(match,nuc_label,bordermask,nucr);
        %tempx=rxcur(match);
        %tempy=rycur(match);
        %denovo_id=denovo_label(tempy,tempx);
        %if denovo_id==0
        %    denovo_id=max(denovo_label(nuc_label==match));
        %end
        %[bordermask,borderflag(match)]=splitdeflections_1(B{denovo_id},bordermask,nucr);
        orderedset=B{match};
        orderedset=[orderedset(end:-1:1,2) orderedset(end:-1:1,1)];
        %bordermask=splitdeflections_4_bwboundaries(orderedset,bordermask,nucr);
        [bordermask,borderflag(match)]=splitdeflections_4_bwboundaries(orderedset,bordermask,nucr);
    end
end
%%% assign new IDs to un-merged cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
splitid=find(borderflag);
mergemaskorg=ismember(nuc_label,splitid);
mergemask=mergemaskorg & ~bordermask;
mergemask=~bwmorph(~mergemask,'diag');
[newlabels,numnew]=bwlabel(mergemask);
if numnew
    %%% extract features of new cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    new_info=struct2cell(regionprops(newlabels,nuc_raw,'Area','Centroid','MeanIntensity')');
    new_area=squeeze(cell2mat(new_info(1,1,:)));
    new_center=squeeze(cell2mat(new_info(2,1,:)))';
    if numnew==1
        new_center=[new_center(1,1) new_center(2,1)];
    end
    new_density=squeeze(cell2mat(new_info(3,1,:)));
    new_mass=new_density.*new_area;
    new_center(:,1)=new_center(:,1)+absjitx;
    new_center(:,2)=new_center(:,2)+absjity;
    %%% update with un-merged data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xcur=[xcur;new_center(:,1)]; xcur(splitid)=NaN;
    ycur=[ycur;new_center(:,2)]; ycur(splitid)=NaN;
    areacur=[areacur;new_area];  areacur(splitid)=NaN;
    masscur=[masscur;new_mass];  masscur(splitid)=NaN;
end
numcur=numel(xcur);
%%% update nuc_label with split cell IDs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_label(mergemaskorg>0)=0;
newlabels=newlabels+numcurorg;       %offset so that all label IDs are new
newlabels(newlabels==numcurorg)=0;
nuc_label=nuc_label+newlabels;

%%% match each cell from previous frame to a cell in the current frame %%%%
prevmatch=ones(maxcellnum,1)*NaN;
curmatch=zeros(numcur,1); %number of times matched
curtracking=ones(numcur,3)*NaN; %[mother, mergingcell1, mergingcell2]
newmerge=ones(maxcellnum,1)*NaN;
for i=1:numprevextant
    p=previd(i);
    neighbors=find(abs(xcur-xprev(p))<maxjump & abs(ycur-yprev(p))<maxjump);
    %%% in case of zero neighbors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(neighbors)
        continue;
    end
    %%% get distances and mass differences of neighbors %%%%%%%%%%%%%%%%%%%
%     if abs(xprev(p)-841)<maxjump && abs(yprev(p)-850)<maxjump && tf==9
%         keyboard;
%     elseif abs(xprev(p)-277)<maxjump && abs(yprev(p)-752)<maxjump && tf==6
%         keyboard;
%     end
    dist=sqrt((xcur(neighbors)-xprev(p)).^2+(ycur(neighbors)-yprev(p)).^2);
    massdiff=(masscur(neighbors)-massprev(p))/massprev(p);
    areadiff=areacur(neighbors)-areaprev(p);
    areadiffnorm=areadiff/areaprev(p);
    %areadiff=areacur(neighbors)-areaprev(p);
    %%% in case of only one neighbor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if numel(neighbors)==1
        if ~ongoingmerge(p) && massdiff>masschangethreshold
        %if ~ongoingmerge(p) && areadiff>debrisarea && areadiffnorm>areachangethreshold
            newmerge(p)=neighbors;
        elseif abs(massdiff)<masschangethreshold
            prevmatch(p)=neighbors;
            curmatch(neighbors)=curmatch(neighbors)+1;
        end
        continue;
    end
    %%% for multiple neighbors, find closest 2 neighbors %%%%%%%%%%%%%%%%%%
    [~,didx]=sort(dist);
    candidateidx=[didx(1);didx(2)];
    candidates=neighbors(candidateidx);
    dist=dist(candidateidx);
    massdiff=massdiff(candidateidx); areadiff=areadiff(candidateidx); areadiffnorm=areadiffnorm(candidateidx);
    %splitcheck=ongoingmerge(p) & ~mergeduration(p) & massdiff(1)<-masschangethreshold & massdiff(2)<-masschangethreshold & abs(sum(massdiff)+1)<masschangethreshold;
    %splitcheck=ongoingmerge(p) & shortmerge(p) & massdiff(1)<-masschangethreshold & massdiff(2)<-masschangethreshold & abs(sum(massdiff)+1)<masschangethreshold;
    splitcheck=ongoingmerge(p) & shortmerge(p) & abs(sum(massdiff)+1)<masschangethreshold;
    daughtercheck=~ongoingmerge(p) & massdiff>daughtermin & massdiff<daughtermax & dist<maxjump & areadiff<0; %H2B: -0.55 to -0.45; NLS: -0.70 to -0.30
    if splitcheck
        %tracking(p,5)=curframe-1; %finish tracking merge
        cell1=tracking(p,2); cell2=tracking(p,3);
        mergedframes=tracking(p,4):(curframe-1);
        premergeframe=find(~isnan(tracedata(cell1,:,1)),1,'last');
        mass1pre=tracedata(cell1,premergeframe,4); mass2pre=tracedata(cell2,premergeframe,4);
        massdiffpre=(mass2pre-mass1pre)/mass1pre;
        massdiffpost=(masscur(candidates(2))-masscur(candidates(1)))/masscur(candidates(1));
        if tracking(cell1,1)==tracking(cell2,1)  %sisters
            tracedata([cell1;cell2],mergedframes,:)=tracedata([p;p],mergedframes,:);
            tracedata([cell1;cell2],mergedframes,[3 4])=tracedata([cell1;cell2],mergedframes,[3 4])/2;
            prevmatch(cell1)=candidates(1);
            prevmatch(cell2)=candidates(2);
            curmatch(candidates)=curmatch(candidates)+1;
        elseif abs(massdiffpre)>masschangethreshold && abs(massdiffpost)>masschangethreshold
            splitorder=isequal(massdiffpre>0,massdiffpost>0);
            if splitorder==1
                cell1post=candidates(1); cell2post=candidates(2);
            else
                cell1post=candidates(2); cell2post=candidates(1);
            end
            %tracedata(cell1,mergedframes,:)=(tracedata(cell1post,curframe,:)-tracedata(cell1,premergeframe,:))*(curframe-mergedframes);
            %tracedata(cell2,mergedframes,:)=(tracedata(cell2post,curframe,:)-tracedata(cell2,premergeframe,:))*(curframe-mergedframes);
            prevmatch(cell1)=cell1post;
            prevmatch(cell2)=cell2post;
            curmatch(candidates)=curmatch(candidates)+1;
        else
            %don't track. will begin tracking as new traces
        end
    elseif sum(daughtercheck)==2
        curtracking(candidates,1)=p; %store mothercell
        curmatch(candidates)=curmatch(candidates)+1;
    elseif ~ongoingmerge(p) && massdiff(1)>masschangethreshold
    %elseif ~ongoingmerge(p) && areadiff(1)>debrisarea*4 && areadiffnorm(1)>areachangethreshold
        newmerge(p)=candidates(1);
    elseif abs(massdiff(1))<masschangethreshold
        prevmatch(p)=candidates(1);
        curmatch(candidates(1))=curmatch(candidates(1))+1;
    end
end
%%% confirm and track new merges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eachmerge=unique(newmerge);
eachmerge(isnan(eachmerge))=[];
for m=eachmerge'
    mergingcells=find(newmerge==m);
    if numel(mergingcells)==2
        massmerge=masscur(m);
        masspre=sum(massprev(mergingcells));
        massdiff=(massmerge-masspre)/masspre;
        if abs(massdiff)<masschangethreshold
            curtracking(m,[2 3])=[mergingcells(1) mergingcells(2)]; %store mergingcells
            curmatch(m)=curmatch(m)+1;
        end
    end
end
%%% remove conflicts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curconflict=find(curmatch>1);
for c=curconflict'
    tempprevs=find(prevmatch==c);
    tempmother=curtracking(c,1);
    if ~isnan(tempmother) %mitoses should be removed
        tempsisters=curtracking(:,1)==tempmother;
        curtracking(tempsisters,1)=NaN;
        curmatch(tempsisters)=0;
    end
    prevmatch(tempprevs)=NaN;
    curmatch(c)=0;
    curtracking(c,:)=NaN;
end
%%% stop tracking split or non-tracked merged cells %%%%%%%%%%%%%%%%%%%%%%%
tracking(ongoingmerge & isnan(prevmatch),5)=curframe-1;
%%% stop tracking chronic merges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chronicmerge=curframe-tracking(:,4)>5 & isnan(tracking(:,5));
% tracking(chronicmerge,5)=curframe;

%%% update tracked info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curdatatracked=ones(maxcellnum,4)*NaN;
relabelidx=ones(numcur,1)*NaN;
%%% add tracked cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matchprevidx=find(~isnan(prevmatch));
matchcuridx=prevmatch(matchprevidx);
curdatatracked(matchprevidx,:)=[xcur(matchcuridx) ycur(matchcuridx) areacur(matchcuridx) masscur(matchcuridx)];
relabelidx(matchcuridx)=matchprevidx;
%%% add daughters cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curtotal=numprevtotal;
daughtercx=find(~isnan(curtracking(:,1)));
numdaughter=numel(daughtercx);
if numdaughter>0
    daughterpx=(curtotal+1):(curtotal+numdaughter);
    curdatatracked(daughterpx,:)=[xcur(daughtercx) ycur(daughtercx) areacur(daughtercx) masscur(daughtercx)];
    tracking(daughterpx,1)=curtracking(daughtercx,1);
    relabelidx(daughtercx)=daughterpx;
    curtotal=curtotal+numdaughter;
end
%%% add newly merged cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mergecx=find(~isnan(curtracking(:,2)));
nummerge=numel(mergecx);
if nummerge>0
    mergepx=(curtotal+1):(curtotal+nummerge);
    curdatatracked(mergepx,:)=[xcur(mergecx) ycur(mergecx) areacur(mergecx) masscur(mergecx)];
    tracking(mergepx,[2 3])=curtracking(mergecx,[2 3]);
    tracking(mergepx,4)=curframe;
    relabelidx(mergecx)=mergepx;
    curtotal=curtotal+nummerge;
end
%%% add non-tracked cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nonmatchcx=find(curmatch==0);
numnonmatched=numel(nonmatchcx);
if numnonmatched>0
    nonmatchpx=(curtotal+1):(curtotal+numnonmatched);
    curdatatracked(nonmatchpx,:)=[xcur(nonmatchcx) ycur(nonmatchcx) areacur(nonmatchcx) masscur(nonmatchcx)];
    relabelidx(nonmatchcx)=nonmatchpx;
end
%%% re-label nuc_label %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_info=regionprops(nuc_label,'PixelIdxList');
for i=1:numcur
    nuc_label(nuc_info(i).PixelIdxList)=relabelidx(i);
end

%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%NOTE: must include extractmask, jitx, and jity in the arguments
%%%%%% view current cell and its prior neighbors %%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(nuc_label);
dxminprev=max([round(xprev(match)-winrad) 1]); dxmaxprev=min([round(xprev(match)+winrad) width]);
dyminprev=max([round(yprev(match)-winrad) 1]); dymaxprev=min([round(yprev(match)+winrad) height]);
dxmincur=round(dxminprev-reljitx); dxmindiff=double((1-dxmincur)*(dxmincur<1)); dxmincur=max([dxmincur 1]);
dxmaxcur=round(dxmaxprev-reljitx); dxmaxdiff=double((dxmaxcur-width)*(dxmaxcur>width)); dxmaxcur=min([dxmaxcur width]);
dymincur=round(dyminprev-reljity); dymindiff=double((1-dymincur)*(dymincur<1)); dymincur=max([dymincur 1]);
dymaxcur=round(dymaxprev-reljity); dymaxdiff=double((dymaxcur-height)*(dymaxcur>height)); dymaxcur=min([dymaxcur height]);
dbmaskcur=bwmorph(nuc_label,'remove');
dbmaskcur=dbmaskcur(dymincur:dymaxcur,dxmincur:dxmaxcur);
dbmaskcur=padarray(dbmaskcur,[dymindiff dxmindiff],'pre');
dbmaskcur=padarray(dbmaskcur,[dymaxdiff dxmaxdiff],'post');
dbmaskprev=extractmask(dyminprev:dymaxprev,dxminprev:dxmaxprev);
dbimage=mat2gray(dbmaskprev);
dbimage(:,:,2)=dbmaskcur;
dbimage(:,:,3)=0;
figure,imshow(dbimage);
%%%%%% view current cell w/ bridge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(nuc_label);
dxminprev=max([round(xprev(match)-winrad) 1]); dxmaxprev=min([round(xprev(match)+winrad) width]);
dyminprev=max([round(yprev(match)-winrad) 1]); dymaxprev=min([round(yprev(match)+winrad) height]);
dxmincur=round(dxminprev-reljitx); dxmindiff=double((1-dxmincur)*(dxmincur<1)); dxmincur=max([dxmincur 1]);
dxmaxcur=round(dxmaxprev-reljitx); dxmaxdiff=double((dxmaxcur-width)*(dxmaxcur>width)); dxmaxcur=min([dxmaxcur width]);
dymincur=round(dyminprev-reljity); dymindiff=double((1-dymincur)*(dymincur<1)); dymincur=max([dymincur 1]);
dymaxcur=round(dymaxprev-reljity); dymaxdiff=double((dymaxcur-height)*(dymaxcur>height)); dymaxcur=min([dymaxcur height]);
dbmaskcur=bwmorph(nuc_label,'remove');
dbmaskcur=dbmaskcur(dymincur:dymaxcur,dxmincur:dxmaxcur);
dbmaskcur=padarray(dbmaskcur,[dymindiff dxmindiff],'pre');
dbmaskcur=padarray(dbmaskcur,[dymaxdiff dxmaxdiff],'post');
dbbridgecur=bordermask(dymincur:dymaxcur,dxmincur:dxmaxcur);
dbbridgecur=padarray(dbbridgecur,[dymindiff dxmindiff],'pre');
dbbridgecur=padarray(dbbridgecur,[dymaxdiff dxmaxdiff],'post');
dbmaskprev=extractmask(dyminprev:dymaxprev,dxminprev:dxmaxprev);
dbimage=dbbridgecur;
dbimage(:,:,2)=dbmaskcur;
dbimage(:,:,3)=0;
figure,imshow(imresize(dbimage,5));
%%%%%% view previous cell and its future neighbors %%%%%%%%%%%%%%%%%%%%%%%%
p=3767;
prevcellandcurrentneighbors(debugpackage,nuc_label,winrad,xprev,yprev,p);
%}