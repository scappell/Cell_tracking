function thicklabelimg = labelthicken(labeledimg,thickenradius)
[height,width]=size(labeledimg);
thickenedmask=bwmorph(labeledimg,'thicken',thickenradius);
centers=round(squeeze(cell2mat(struct2cell(regionprops(thickenedmask,'Centroid')')))');
centers_idx=sub2ind([height,width],centers(:,2),centers(:,1));
center_label=labeledimg(centers_idx);
numobjects=size(centers,1);
pixels=regionprops(thickenedmask,'PixelIdxList');
thicklabelimg=zeros(height,width);
for i=1:numobjects
    thicklabelimg(pixels(i).PixelIdxList)=center_label(i);
end
end