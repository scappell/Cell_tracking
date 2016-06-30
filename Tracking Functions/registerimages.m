function [alignxrel,alignyrel]=registerimages(image1,image2)
[height,width]=size(image1);
crosscorrscore=abs(normxcorr2(image1,image2));
[~,idx]=max(crosscorrscore(:));
[y,x]=ind2sub(size(crosscorrscore),idx);
alignxrel=width-x;
alignyrel=height-y;