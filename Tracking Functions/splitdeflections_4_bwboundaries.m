function [bordermask,bridgeflag]=splitdeflections_4_bwboundaries(orderedset,bordermask,nucr)
bridgeflag=0; %returned as 1 if deflections are bridged
nucr=round(nucr/4)*4; %make sure nucr is a multiple of 4
perilength=size(orderedset,1);
%%% detect deflection vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vIdx=getdeflections(orderedset,nucr); %returns boundary indices
%%% count vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vnum=length(vIdx);
if vnum<2
    return; %if less than two vertices are detected, exit function
end
%%% calculate perimeter distance between adjacent vertices %%%%%%%%%%%%%%%%
periIdx=vIdx;
periIdxadj1=[periIdx(2:end);perilength+periIdx(1)];
pairperi1=periIdxadj1-periIdx;
%%% pair and bridge vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while vnum>=2
    vpos=orderedset(vIdx,:);
    %%% Determine adjacent vertices that define the highest curvature %%%%%
    vposadj1=[vpos(2:end,:);vpos(1,:)];
    pair1=vposadj1-vpos;
    pairdist1=sqrt(sum(pair1.^2,2));
    curvature1=pairperi1./pairdist1;
    [bestcurve,curve1idx]=sort(curvature1);
    % If best curvature is too low, stop segmenting this object.
    if bestcurve(end)<2
        break
    end
    bestcurveidx=curve1idx(end);
    if bestcurveidx==vnum
        bestcurveidxadj=1;
    else
        bestcurveidxadj=bestcurveidx+1;
    end
    % If this point is reached, a split will be performed, so mark it.
    bridgeflag=1;
    %%% Bridge the vertices defining the best curvature %%%%%%%%%%%%%%%%%%%
    [bx,by]=bridge(vpos(bestcurveidx,:),vpos(bestcurveidxadj,:));
    % [bx,by] can equal NaN if the vertices were the same position (occurs
    % when two objects share a single coordinate. If this is the case,
    % skip mapping the bridge and proceed with perimeter update and vertex
    % removal.
    if ~isnan(bx)
        for bci=1:length(bx)
            %bridgeflag=1;
            bordermask(by(bci),bx(bci))=1;
        end
    end
    %%% assign new perimeter distances & remove old vertices %%%%%%%%%%%%%%
    previdx=bestcurveidx-1;
    if previdx==0
        previdx=vnum;
    end
    % Given vertex 3-4 gave best curvature and is now bridged, define the
    % perimeter from vertex 2 to 5: p(2-5)=p(2-3)+bridge+p(4-5).
    pairperi1(previdx)=pairperi1(previdx)+length(bx)-1+pairperi1(bestcurveidxadj);
    % Remove the vertices and perimeters of the vertices defining the best
    % curve.
    vIdx([bestcurveidx,bestcurveidxadj])=[];
    pairperi1([bestcurveidx,bestcurveidxadj])=[];
    vnum=length(vIdx);
end
%keyboard;
end
%%% debug: visualize deflections on boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% vpos=orderedset(vIdx,:);
% for vc=1:size(vpos,1)
%     bordermask(vpos(vc,2),vpos(vc,1))=1;
% end

tempbridgeimage=zeros(size(scmask));
for bci=1:length(bcx)
   tempbridgeimage(bcy(bci),bcx(bci))=1;
end
tempimage=mat2gray(scmask);
tempimage(:,:,2)=tempbridgeimage;
tempimage(:,:,3)=0;
imshow(tempimage);
%}
