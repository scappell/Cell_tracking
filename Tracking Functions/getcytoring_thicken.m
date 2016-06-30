function outerlabel=getcytoring_thicken(nuc_label,innerrad,outerrad,image)
innermask=bwmorph(nuc_label,'thicken',innerrad);
outerlabel=labelthicken(nuc_label,outerrad);
outerlabelorg=outerlabel;
outerlabel(innermask)=0;

%%% Remove cell-cell contact areas from ring calculations %%%%%%%%%%%%%%%%%
borders=bwmorph(outerlabel,'bothat');
borders=imdilate(borders,strel('disk',2,0));
outerlabel=outerlabel.*~borders;

%%% for absent rings, add a single pixel at centroid as a marker later %%%%
uniqueorg=unique(outerlabelorg);
uniquebordercleared=unique(outerlabel);
noring=uniqueorg(~ismember(uniqueorg,uniquebordercleared));
if ~isempty(noring)
    centers=round(squeeze(cell2mat(struct2cell(regionprops(outerlabelorg,'Centroid')')))');
    [height,width]=size(nuc_label);
    centers_idx=sub2ind([height,width],centers(:,2),centers(:,1));
    for i=noring'
        %ring_label(cytoring_rzb==i)=i;
        %outerlabel(cytoring==i)=i;
        outerlabel(centers_idx(i))=i;
    end
end
% keyboard;
%{
%%% visualization for debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempframe=imadjust(mat2gray(image));
innerdilate=imdilate(innermask,strel('disk',1));
outerlabel(innerdilate)=0;
outershell=bwmorph(outerlabel,'remove');
tempframe(:,:,2)=outershell;
%tempframe(:,:,2)=innershell | outershell;
%tempframe(:,:,3)=bwmorph(nuc_label,'remove');
%tempframe(:,:,3)=bwmorph(borders,'remove');
%tempframe(:,:,2)=ring_label>0;
tempframe(:,:,3)=0;
figure,
imshow(tempframe);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
end