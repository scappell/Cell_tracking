function newmask=secondthresh(image,blurradius,mask,boulderarea)
blur=imfilter(log(image),fspecial('disk',blurradius),'symmetric');
boulder_mask=bwareaopen(mask,boulderarea);
bouldervals=image(boulder_mask);
if isempty(bouldervals)
    newmask=mask;
    return;
end
normlogvals=mat2gray(log(bouldervals));
higherthresh=graythresh(normlogvals);
loghigherthresh=higherthresh*range(log(bouldervals))+min(log(bouldervals));
higher_mask=blur>loghigherthresh;
higher_mask(~boulder_mask)=0;
mask(boulder_mask)=0;
newmask=mask | higher_mask;
newmask=imfill(newmask,'holes');
end